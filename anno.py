import os
import subprocess
import argparse
import yaml
import pandas as pd
import duckdb
import gzip
import shutil
from Bio import SeqIO
from multiprocessing import Pool

# 默认配置文件
DEFAULT_CONFIG = {
    "fasta_file": "input.fasta",
    "database": {
        "source_dir": "./databases",
        "uniprot_sprot": "uniprot_sprot.fasta.gz",
        "uniref90": "uniref90.fasta.gz",
        "idmapping": "idmapping_selected.tab.gz",
        "ncRNA": "rnacentral_active.fasta.gz",
        "cdd": "Cdd_LE.tar.gz"
    },
    "diamond": {
        "evalue": 1e-5,
        "outfmt": 6,
        "sensitive": True
    },
    "blastn": {
        "evalue": 1e-5,
        "outfmt": 6
    },
    "rpstblastn": {
        "evalue": 1e-5,
        "outfmt": 6
    },
    "cpat": {
        "hexamer": "human_hexamer.tsv",  # 需要用户提供 CPAT 训练文件
        "logit_model": "human_logit.RData"
    },
    "transdecoder": {
        "min_length": 100
    },
    "output": {
        "dir": "./output",
        "diamond_out": "diamond.out",
        "blastn_out": "blastn.out",
        "rpstblastn_out": "rpstblastn.out",
        "annotations": "annotations.tsv",
        "gff_dir": "./gff",
        "stats_dir": "./stats"
    }
}

def generate_config_template(user_dir):
    """生成配置文件模板"""
    config_path = os.path.join(user_dir, "annocript_config.yaml")
    if not os.path.exists(config_path):
        with open(config_path, "w") as f:
            yaml.dump(DEFAULT_CONFIG, f, default_flow_style=False)
        print(f"Generated config template at {config_path}. Please edit it and rerun.")
        exit(0)

def load_config(user_dir):
    """加载配置"""
    config_path = os.path.join(user_dir, "annocript_config.yaml")
    generate_config_template(user_dir)
    with open(config_path) as f:
        return yaml.safe_load(f)

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description="Annocript: Transcriptome Annotation Tool")
    parser.add_argument("--fasta", help="Input FASTA file", default=DEFAULT_CONFIG["fasta_file"])
    parser.add_argument("--db_dir", help="Database directory", default=DEFAULT_CONFIG["database"]["source_dir"])
    parser.add_argument("--threads", type=int, help="Number of CPU threads", default=8)
    parser.add_argument("--do_db_creation", action="store_true", help="Create database (skipped with pre-downloaded files)")
    parser.add_argument("--do_execute_programs", action="store_true", help="Execute alignment programs", default=True)
    parser.add_argument("--do_blastx", action="store_true", help="Run DIAMOND BLASTX", default=True)
    parser.add_argument("--do_blastn", action="store_true", help="Run BLASTN", default=True)
    parser.add_argument("--do_rpstblastn", action="store_true", help="Run RPSBLASTN", default=True)
    parser.add_argument("--do_lnc_prediction", action="store_true", help="Predict lncRNA with CPAT", default=True)
    parser.add_argument("--do_dna2pep", action="store_true", help="Search ORFs with TransDecoder", default=True)
    parser.add_argument("--do_build_output", action="store_true", help="Build final output", default=True)
    parser.add_argument("--extract_stats", action="store_true", help="Generate statistics", default=True)
    return parser.parse_args()

def check_files(config, args):
    """检查预下载文件"""
    db_dir = args.db_dir
    required_files = [
        config["database"]["uniprot_sprot"],
        config["database"]["uniref90"],
        config["database"]["idmapping"],
        config["database"]["ncRNA"],
        config["database"]["cdd"]
    ]
    for f in required_files:
        if not os.path.exists(os.path.join(db_dir, f)):
            raise FileNotFoundError(f"Missing {f} in {db_dir}")
    print("All required database files found.")

def build_index(db_file, output_db, tool, threads):
    """构建索引（DIAMOND 或 BLAST）"""
    if db_file.endswith(".gz"):
        uncompressed = output_db + ".fasta"
        if not os.path.exists(uncompressed):
            with gzip.open(db_file, "rb") as f_in, open(uncompressed, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        db_file = uncompressed
    
    if tool == "diamond" and not os.path.exists(f"{output_db}.dmnd"):
        cmd = ["diamond", "makedb", "--in", db_file, "--db", output_db, "--threads", str(threads)]
        subprocess.run(cmd, check=True)
        print(f"DIAMOND index built: {output_db}.dmnd")
    elif tool == "blast" and not os.path.exists(f"{output_db}.nhr"):
        cmd = ["makeblastdb", "-in", db_file, "-dbtype", "nucl", "-out", output_db]
        subprocess.run(cmd, check=True)
        print(f"BLAST index built: {output_db}")
    return output_db

def run_alignment(args_tuple):
    """并行运行比对任务"""
    tool, config, args, db_index, output_file = args_tuple
    if tool == "diamond":
        cmd = [
            "diamond", "blastx",
            "--query", args.fasta,
            "--db", db_index,
            "--out", output_file,
            "--evalue", str(config["diamond"]["evalue"]),
            "--threads", str(args.threads // 3),  # 分配线程
            "--outfmt", str(config["diamond"]["outfmt"])
        ]
        if config["diamond"]["sensitive"]:
            cmd.append("--sensitive")
    elif tool == "blastn":
        cmd = [
            "blastn",
            "-query", args.fasta,
            "-db", db_index,
            "-out", output_file,
            "-evalue", str(config["blastn"]["evalue"]),
            "-num_threads", str(args.threads // 3),
            "-outfmt", str(config["blastn"]["outfmt"])
        ]
    elif tool == "rpstblastn":
        cmd = [
            "rpstblastn",
            "-query", args.fasta,
            "-db", db_index,
            "-out", output_file,
            "-evalue", str(config["rpstblastn"]["evalue"]),
            "-num_threads", str(args.threads // 3),
            "-outfmt", str(config["rpstblastn"]["outfmt"])
        ]
    subprocess.run(cmd, check=True)
    print(f"{tool} alignment completed: {output_file}")
    return output_file

def predict_lncRNA(config, args, diamond_out, blastn_out):
    """使用 CPAT 预测 lncRNA"""
    output_dir = config["output"]["dir"]
    cpat_out = os.path.join(output_dir, "cpat_output")
    
    # 合并比对结果，过滤可能的编码序列
    diamond_hits = set(pd.read_csv(diamond_out, sep="\t", header=None)[0])
    blastn_hits = set(pd.read_csv(blastn_out, sep="\t", header=None)[0])
    lnc_candidates = os.path.join(output_dir, "lnc_candidates.fasta")
    
    with open(lnc_candidates, "w") as out:
        for record in SeqIO.parse(args.fasta, "fasta"):
            if record.id not in diamond_hits and record.id in blastn_hits:
                SeqIO.write(record, out, "fasta")
    
    # 运行 CPAT
    cmd = [
        "cpat.py",
        "-x", os.path.join(args.db_dir, config["cpat"]["hexamer"]),
        "-d", os.path.join(args.db_dir, config["cpat"]["logit_model"]),
        "-g", lnc_candidates,
        "-o", cpat_out
    ]
    subprocess.run(cmd, check=True)
    print(f"CPAT lncRNA prediction completed: {cpat_out}.ORF_seqs.fa")
    return f"{cpat_out}.ORF_seqs.fa"

def dna2pep(config, args):
    """使用 TransDecoder 搜索 ORF"""
    output_dir = config["output"]["dir"]
    transdecoder_dir = os.path.join(output_dir, "transdecoder")
    os.makedirs(transdecoder_dir, exist_ok=True)
    
    # 运行 TransDecoder.LongOrfs
    cmd = [
        "TransDecoder.LongOrfs",
        "-t", args.fasta,
        "-m", str(config["transdecoder"]["min_length"]),
        "--output_dir", transdecoder_dir
    ]
    subprocess.run(cmd, check=True)
    
    # 运行 TransDecoder.Predict
    cmd = [
        "TransDecoder.Predict",
        "-t", args.fasta,
        "--output_dir", transdecoder_dir
    ]
    subprocess.run(cmd, check=True)
    print(f"TransDecoder ORF search completed: {transdecoder_dir}/{args.fasta}.transdecoder.pep")
    return f"{transdecoder_dir}/{args.fasta}.transdecoder.pep"

def prepare_idmapping(config, args):
    """将 idmapping 转为 Parquet"""
    db_dir = args.db_dir
    idmapping_gz = os.path.join(db_dir, config["database"]["idmapping"])
    idmapping_parquet = os.path.join(db_dir, "idmapping.parquet")
    
    if not os.path.exists(idmapping_parquet):
        df = pd.read_csv(idmapping_gz, sep="\t", header=None, 
                         names=["uniprot_id", "type", "value"], compression="gzip")
        df.to_parquet(idmapping_parquet)
    return idmapping_parquet

def build_output(config, args, diamond_out, blastn_out, rpstblastn_out, idmapping_parquet):
    """生成最终输出（表格和 GFF3）"""
    output_dir = config["output"]["dir"]
    annotations_file = os.path.join(output_dir, config["output"]["annotations"])
    gff_dir = config["output"]["gff_dir"]
    os.makedirs(gff_dir, exist_ok=True)
    
    conn = duckdb.connect()
    conn.execute(f"""
        CREATE TABLE diamond AS 
        SELECT * FROM read_csv_auto('{diamond_out}', delim='\t', 
            names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                   'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    """)
    conn.execute(f"""
        CREATE TABLE blastn AS 
        SELECT * FROM read_csv_auto('{blastn_out}', delim='\t', 
            names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                   'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    """)
    conn.execute(f"""
        CREATE TABLE rpstblastn AS 
        SELECT * FROM read_csv_auto('{rpstblastn_out}', delim='\t', 
            names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                   'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    """)
    conn.execute(f"CREATE TABLE idmapping AS SELECT * FROM '{idmapping_parquet}'")
    
    # 生成表格
    query = """
    SELECT d.qseqid, d.sseqid AS protein, d.evalue, d.bitscore,
           b.sseqid AS ncRNA, b.evalue AS nc_evalue,
           r.sseqid AS domain, r.evalue AS domain_evalue,
           m.value AS go_term
    FROM diamond d
    LEFT JOIN blastn b ON d.qseqid = b.qseqid
    LEFT JOIN rpstblastn r ON d.qseqid = r.qseqid
    LEFT JOIN idmapping m ON d.sseqid = m.uniprot_id AND m.type = 'GO'
    """
    result = conn.execute(query).fetchdf()
    result.to_csv(annotations_file, sep="\t", index=False)
    
    # 生成 GFF3
    with open(os.path.join(gff_dir, "annotations.gff3"), "w") as gff:
        gff.write("##gff-version 3\n")
        for _, row in result.iterrows():
            if pd.notna(row["protein"]):
                gff.write(f"{row['qseqid']}\tDIAMOND\tmatch\t{row['qstart']}\t{row['qend']}\t"
                          f"{row['bitscore']}\t+\t.\tID={row['sseqid']};evalue={row['evalue']}\n")
    
    conn.close()
    print(f"Output generated: {annotations_file}, GFF3 in {gff_dir}")

def extract_statistics(config, annotations_file):
    """生成统计"""
    stats_dir = config["output"]["stats_dir"]
    os.makedirs(stats_dir, exist_ok=True)
    df = pd.read_csv(annotations_file, sep="\t")
    
    stats = {
        "total_sequences": len(df["qseqid"].unique()),
        "with_protein": len(df[df["protein"].notna()]["qseqid"].unique()),
        "with_ncRNA": len(df[df["ncRNA"].notna()]["qseqid"].unique()),
        "with_domain": len(df[df["domain"].notna()]["qseqid"].unique())
    }
    with open(os.path.join(stats_dir, "stats.txt"), "w") as f:
        yaml.dump(stats, f)
    print(f"Statistics saved to {stats_dir}/stats.txt")

def main():
    args = parse_arguments()
    config = load_config(os.getcwd())
    
    # 检查文件
    check_files(config, args)
    
    if args.do_execute_programs:
        diamond_out = blastn_out = rpstblastn_out = None
        output_dir = config["output"]["dir"]
        os.makedirs(output_dir, exist_ok=True)
        
        # 并行运行比对
        tasks = []
        if args.do_blastx:
            db_index = build_index(os.path.join(args.db_dir, config["database"]["uniprot_sprot"]), 
                                  os.path.join(args.db_dir, "uniprot_sprot"), "diamond", args.threads)
            tasks.append(("diamond", config, args, db_index, os.path.join(output_dir, config["output"]["diamond_out"])))
        if args.do_blastn:
            nc_index = build_index(os.path.join(args.db_dir, config["database"]["ncRNA"]), 
                                  os.path.join(args.db_dir, "ncRNA"), "blast", args.threads)
            tasks.append(("blastn", config, args, nc_index, os.path.join(output_dir, config["output"]["blastn_out"])))
        if args.do_rpstblastn:
            cdd_dir = os.path.join(args.db_dir, "cdd")
            if not os.path.exists(cdd_dir):
                os.makedirs(cdd_dir)
                subprocess.run(["tar", "-xzf", os.path.join(args.db_dir, config["database"]["cdd"]), "-C", cdd_dir], check=True)
            tasks.append(("rpstblastn", config, args, os.path.join(cdd_dir, "Cdd"), 
                         os.path.join(output_dir, config["output"]["rpstblastn_out"])))
        
        with Pool(processes=min(len(tasks), args.threads)) as pool:
            results = pool.map(run_alignment, tasks)
            diamond_out, blastn_out, rpstblastn_out = (results + [None, None, None])[:3]
        
        if args.do_lnc_prediction and diamond_out and blastn_out:
            predict_lncRNA(config, args, diamond_out, blastn_out)
        if args.do_dna2pep:
            dna2pep(config, args)
    
    if args.do_build_output:
        idmapping_parquet = prepare_idmapping(config, args)
        build_output(config, args, diamond_out, blastn_out, rpstblastn_out, idmapping_parquet)
    
    if args.extract_stats:
        extract_statistics(config, os.path.join(config["output"]["dir"], config["output"]["annotations"]))

if __name__ == "__main__":
    main()
