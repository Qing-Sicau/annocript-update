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

DEFAULT_CONFIG = {
    "fasta_file": "input.fasta",
    "database": {
        "source_dir": "./databases",
        "uniprot_sprot": "uniprot_sprot.fasta.gz",
        "uniref90": "uniref90.fasta.gz",
        "idmapping": "idmapping_selected.tab.gz",
        "ncRNA": "rnacentral_active.fasta.gz",
        "cdd": "Cdd_LE.tar.gz",
        "uniprot_dat": "uniprot_sprot.dat.gz" 
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
        "hexamer": "Arabidopsis_hexamer.tsv",
        "logit_model": "Arabidopsis_logit.RData"
    },
    "transdecoder": {
        "min_length": 50
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
    config_path = os.path.join(user_dir, "annocript_config.yaml")
    if not os.path.exists(config_path):
        with open(config_path, "w") as f:
            yaml.dump(DEFAULT_CONFIG, f, default_flow_style=False)
        print(f"Generated config template at {config_path}. Please edit it and rerun.")
        exit(0)

def load_config(user_dir):
    config_path = os.path.join(user_dir, "annocript_config.yaml")
    generate_config_template(user_dir)
    with open(config_path) as f:
        return yaml.safe_load(f)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Annocript: Plant Transcriptome Annotation Tool with Function Descriptions")
    parser.add_argument("--fasta", help="Input FASTA file", default=DEFAULT_CONFIG["fasta_file"])
    parser.add_argument("--db_dir", help="Database directory", default=DEFAULT_CONFIG["database"]["source_dir"])
    parser.add_argument("--threads", type=int, help="Number of CPU threads", default=8)
    parser.add_argument("--do_db_creation", action="store_true", help="Create database (skipped with pre-downloaded files)")
    parser.add_argument("--do_execute_programs", action="store_true", help="Execute alignment programs", default=True)
    parser.add_argument("--do_blastx", action="store_true", help="Run DIAMOND BLASTX", default=True)
    parser.add_argument("--do_blastn", action="store_true", help="Run BLASTN", default=True)
    parser.add_argument("--do_rpstblastn", action="store_true", help="Run RPSBLASTN", default=True)
    parser.add_argument("--do_lnc_prediction", action="store_true", help="Predict lncRNA with CPAT (plant)", default=True)
    parser.add_argument("--do_dna2pep", action="store_true", help="Search ORFs with TransDecoder (plant)", default=True)
    parser.add_argument("--do_build_output", action="store_true", help="Build final output", default=True)
    parser.add_argument("--extract_stats", action="store_true", help="Generate statistics", default=True)
    parser.add_argument("--do_kegg_annotation", action="store_true", help="Perform KEGG annotation", default=False)
    parser.add_argument("--do_extended_annotation", action="store_true", help="Perform Pfam, EC, TAIR, UniPathway annotations", default=False)
    parser.add_argument("--do_function_description", action="store_true", help="Add detailed gene function descriptions", default=False)
    return parser.parse_args()

def check_files(config, args):
    db_dir = args.db_dir
    required_files = [
        config["database"]["uniprot_sprot"],
        config["database"]["uniref90"],
        config["database"]["idmapping"],
        config["database"]["ncRNA"],
        config["database"]["cdd"],
        config["cpat"]["hexamer"],
        config["cpat"]["logit_model"]
    ]
    if args.do_function_description:
        required_files.append(config["database"]["uniprot_dat"])
    if args.do_kegg_annotation and os.path.exists(os.path.join(db_dir, config["kegg"]["pathway_file"])):
        required_files.append(config["kegg"]["pathway_file"])
    for f in required_files:
        if not os.path.exists(os.path.join(db_dir, f)):
            raise FileNotFoundError(f"Missing {f} in {db_dir}")
    print("All required database files found.")

def build_index(db_file, output_db, tool, threads):
    if db_file.endswith(".gz"):
        uncompressed = output_db + ".fasta"
        if not os.path.exists(uncompressed):
            with gzip.open(db_file, "rb") as f_in, open(uncompressed, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        db_file = uncompressed
    
    if tool == "diamond" and not os.path.exists(f"{output_db}.dmnd"):
        cmd = ["diamond", "makedb", "--in", db_file, "--db", output_db, "--threads", str(threads)]
        subprocess.run(cmd, check=True)
    elif tool == "blast" and not os.path.exists(f"{output_db}.nhr"):
        cmd = ["makeblastdb", "-in", db_file, "-dbtype", "nucl", "-out", output_db]
        subprocess.run(cmd, check=True)
    return output_db

def run_alignment(args_tuple):
    tool, config, args, db_index, output_file = args_tuple
    threads_per_task = max(1, args.threads // 3)
    if tool == "diamond":
        cmd = [
            "diamond", "blastx",
            "--query", args.fasta,
            "--db", db_index,
            "--out", output_file,
            "--evalue", str(config["diamond"]["evalue"]),
            "--threads", str(threads_per_task),
            "--outfmt", str(config["diamond"]["outfmt"]),
            "--sensitive" if config["diamond"]["sensitive"] else ""
        ]
    elif tool == "blastn":
        cmd = [
            "blastn",
            "-query", args.fasta,
            "-db", db_index,
            "-out", output_file,
            "-evalue", str(config["blastn"]["evalue"]),
            "-num_threads", str(threads_per_task),
            "-outfmt", str(config["blastn"]["outfmt"])
        ]
    elif tool == "rpstblastn":
        cmd = [
            "rpstblastn",
            "-query", args.fasta,
            "-db", db_index,
            "-out", output_file,
            "-evalue", str(config["rpstblastn"]["evalue"]),
            "-num_threads", str(threads_per_task),
            "-outfmt", str(config["rpstblastn"]["outfmt"])
        ]
    subprocess.run([arg for arg in cmd if arg], check=True)
    print(f"{tool} alignment completed: {output_file}")
    return output_file

def predict_lncRNA(config, args, diamond_out, blastn_out):
    output_dir = config["output"]["dir"]
    cpat_out = os.path.join(output_dir, "cpat_output")
    lnc_candidates = os.path.join(output_dir, "lnc_candidates.fasta")
    
    diamond_hits = set(pd.read_csv(diamond_out, sep="\t", header=None)[0])
    blastn_hits = set(pd.read_csv(blastn_out, sep="\t", header=None)[0])
    with open(lnc_candidates, "w") as out:
        for record in SeqIO.parse(args.fasta, "fasta"):
            if record.id not in diamond_hits and record.id in blastn_hits:
                SeqIO.write(record, out, "fasta")
    
    cmd = [
        "cpat.py",
        "-x", os.path.join(args.db_dir, config["cpat"]["hexamer"]),
        "-d", os.path.join(args.db_dir, config["cpat"]["logit_model"]),
        "-g", lnc_candidates,
        "-o", cpat_out
    ]
    subprocess.run(cmd, check=True)
    
    lnc_file = os.path.join(output_dir, "lncRNA.fasta")
    df = pd.read_csv(f"{cpat_out}.txt", sep="\t")
    lnc_ids = df[df["coding_prob"] < 0.5]["mRNA"].tolist()
    with open(lnc_file, "w") as out:
        for record in SeqIO.parse(lnc_candidates, "fasta"):
            if record.id in lnc_ids:
                SeqIO.write(record, out, "fasta")
    print(f"CPAT lncRNA prediction completed: {lnc_file}")
    return lnc_file

def dna2pep(config, args):
    output_dir = config["output"]["dir"]
    transdecoder_dir = os.path.join(output_dir, "transdecoder")
    os.makedirs(transdecoder_dir, exist_ok=True)
    
    cmd = [
        "TransDecoder.LongOrfs",
        "-t", args.fasta,
        "-m", str(config["transdecoder"]["min_length"]),
        "--output_dir", transdecoder_dir
    ]
    subprocess.run(cmd, check=True)
    
    cmd = [
        "TransDecoder.Predict",
        "-t", args.fasta,
        "--output_dir", transdecoder_dir,
        "--single_best_only"
    ]
    subprocess.run(cmd, check=True)
    print(f"TransDecoder ORF search completed: {transdecoder_dir}/{args.fasta}.transdecoder.pep")
    return f"{transdecoder_dir}/{args.fasta}.transdecoder.pep"

def prepare_idmapping(config, args):
    db_dir = args.db_dir
    idmapping_gz = os.path.join(db_dir, config["database"]["idmapping"])
    idmapping_parquet = os.path.join(db_dir, "idmapping.parquet")
    
    if not os.path.exists(idmapping_parquet):
        df = pd.read_csv(idmapping_gz, sep="\t", header=None, 
                         names=["uniprot_id", "type", "value"], compression="gzip")
        df.to_parquet(idmapping_parquet)
    return idmapping_parquet

def parse_uniprot_dat(dat_file):
    """解析 uniprot_sprot.dat.gz 
    descriptions = {}
    with gzip.open(dat_file, "rt") as f:
        uniprot_id = None
        desc = []
        for line in f:
            if line.startswith("ID"):
                uniprot_id = line.split()[1]
            elif line.startswith("DE   RecName: Full="):
                desc.append(line.split("=", 1)[1].strip(";"))
            elif line.startswith("CC   -!- FUNCTION:"):
                desc.append(line.split(":", 1)[1].strip())
            elif line.startswith("//"):
                if uniprot_id and desc:
                    descriptions[uniprot_id] = " ".join(desc)
                uniprot_id = None
                desc = []
    return descriptions

def build_output(config, args, diamond_out, blastn_out, rpstblastn_out, idmapping_parquet):
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
    
    query = """
    SELECT d.qseqid, d.sseqid AS protein, d.evalue, d.bitscore,
           b.sseqid AS ncRNA, b.evalue AS nc_evalue,
           r.sseqid AS domain, r.evalue AS domain_evalue
    """
    if args.do_kegg_annotation or args.do_extended_annotation:
        query += """,
           GROUP_CONCAT(DISTINCT gm.value) AS go_terms,
           GROUP_CONCAT(DISTINCT km.value) AS kegg_ids
        """
    if args.do_extended_annotation:
        query += """,
           GROUP_CONCAT(DISTINCT pm.value) AS pfam_ids,
           GROUP_CONCAT(DISTINCT em.value) AS ec_numbers,
           GROUP_CONCAT(DISTINCT tm.value) AS tair_ids,
           GROUP_CONCAT(DISTINCT um.value) AS unipathway_ids
        """
    query += """
    FROM diamond d
    LEFT JOIN blastn b ON d.qseqid = b.qseqid
    LEFT JOIN rpstblastn r ON d.qseqid = r.qseqid
    """
    if args.do_kegg_annotation or args.do_extended_annotation:
        query += """
        LEFT JOIN idmapping gm ON d.sseqid = gm.uniprot_id AND gm.type = 'GO'
        LEFT JOIN idmapping km ON d.sseqid = km.uniprot_id AND km.type = 'KEGG'
        """
    if args.do_extended_annotation:
        query += """
        LEFT JOIN idmapping pm ON d.sseqid = pm.uniprot_id AND pm.type = 'Pfam'
        LEFT JOIN idmapping em ON d.sseqid = em.uniprot_id AND em.type = 'EC'
        LEFT JOIN idmapping tm ON d.sseqid = tm.uniprot_id AND tm.type = 'TAIR'
        LEFT JOIN idmapping um ON d.sseqid = um.uniprot_id AND um.type = 'UniPathway'
        """
    query += "GROUP BY ALL"
    
    result = conn.execute(query).fetchdf()
    
    if args.do_kegg_annotation and os.path.exists(os.path.join(args.db_dir, config["kegg"]["pathway_file"])):
        pathway_file = os.path.join(args.db_dir, config["kegg"]["pathway_file"])
        pathways = {}
        with open(pathway_file) as f:
            for line in f:
                if line.startswith("PATH:"):
                    parts = line.strip().split("\t")
                    if len(parts) > 1:
                        kegg_id = parts[0].replace("PATH:", "")
                        pathways[kegg_id] = parts[1]
        
        def map_pathways(kegg_ids):
            if pd.isna(kegg_ids):
                return None
            return ";".join([f"{kid} ({pathways.get(kid.split(':')[-1], 'Unknown')})" for kid in kegg_ids.split(",")])
        
        result["kegg_pathways"] = result["kegg_ids"].apply(map_pathways)
    else:
        result["kegg_pathways"] = result.get("kegg_ids", None)
    
    # 添加详细功能描述
    if args.do_function_description:
        dat_file = os.path.join(args.db_dir, config["database"]["uniprot_dat"])
        descriptions = parse_uniprot_dat(dat_file)
        result["function_description"] = result["protein"].map(lambda x: descriptions.get(x, "No description available"))
    
    result.to_csv(annotations_file, sep="\t", index=False)
    
    with open(os.path.join(gff_dir, "annotations.gff3"), "w") as gff:
        gff.write("##gff-version 3\n")
        for _, row in result.iterrows():
            if pd.notna(row["protein"]):
                gff.write(f"{row['qseqid']}\tDIAMOND\tmatch\t{row['qstart']}\t{row['qend']}\t"
                          f"{row['bitscore']}\t+\t.\tID={row['sseqid']};evalue={row['evalue']}\n")
    
    conn.close()
    print(f"Output with function descriptions generated: {annotations_file}, GFF3 in {gff_dir}")
    return annotations_file

def extract_statistics(config, annotations_file):
    stats_dir = config["output"]["stats_dir"]
    os.makedirs(stats_dir, exist_ok=True)
    df = pd.read_csv(annotations_file, sep="\t")
    
    stats = {
        "total_sequences": len(df["qseqid"].unique()),
        "with_protein": len(df[df["protein"].notna()]["qseqid"].unique()),
        "with_ncRNA": len(df[df["ncRNA"].notna()]["qseqid"].unique()),
        "with_domain": len(df[df["domain"].notna()]["qseqid"].unique()),
        "with_kegg": len(df[df["kegg_ids"].notna()]["qseqid"].unique()) if "kegg_ids" in df else 0,
        "with_pfam": len(df[df["pfam_ids"].notna()]["qseqid"].unique()) if "pfam_ids" in df else 0,
        "with_ec": len(df[df["ec_numbers"].notna()]["qseqid"].unique()) if "ec_numbers" in df else 0,
        "with_tair": len(df[df["tair_ids"].notna()]["qseqid"].unique()) if "tair_ids" in df else 0,
        "with_unipathway": len(df[df["unipathway_ids"].notna()]["qseqid"].unique()) if "unipathway_ids" in df else 0,
        "with_function": len(df[df["function_description"].notna()]["qseqid"].unique()) if "function_description" in df else 0
    }
    with open(os.path.join(stats_dir, "stats.txt"), "w") as f:
        yaml.dump(stats, f)
    print(f"Statistics saved to {stats_dir}/stats.txt")

def main():
    args = parse_arguments()
    config = load_config(os.getcwd())
    
    check_files(config, args)
    
    if args.do_execute_programs:
        diamond_out = blastn_out = rpstblastn_out = None
        output_dir = config["output"]["dir"]
        os.makedirs(output_dir, exist_ok=True)
        
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
    
    annotations_file = None
    if args.do_build_output:
        idmapping_parquet = prepare_idmapping(config, args)
        annotations_file = build_output(config, args, diamond_out, blastn_out, rpstblastn_out, idmapping_parquet)
    
    if args.extract_stats and annotations_file:
        extract_statistics(config, annotations_file)

if __name__ == "__main__":
    main()
