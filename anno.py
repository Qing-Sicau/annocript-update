import os
import subprocess
import argparse
import yaml
import pandas as pd
import duckdb
import gzip
import shutil
import tarfile
from Bio import SeqIO
from multiprocessing import Pool
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Default configuration
DEFAULT_CONFIG = {
    "fasta_file": "input.fasta",
    "database": {
        "source_dir": "./databases",
        "uniprot_sprot": "uniprot_sprot.fasta.gz",
        "uniref90": "uniref90.fasta.gz",
        "uniprot_trembl": "uniprot_trembl.fasta.gz",
        "ncRNA": "rnacentral_active.fasta.gz",
        "idmapping": "idmapping_selected.tab.gz",
        "cdd": "Cdd_LE.tar.gz",
        "uniprot_dat": "uniprot_sprot.dat.gz",
        "pathway_file": "pathway.txt",
        "go_obo": "go.obo",
        "enzyme_dat": "enzyme.dat",
        "pfam2go": "pfam2go.txt"
    },
    "diamond": {
        "evalue": 1e-5,
        "outfmt": 6,
        "sensitive": False
    },
    "blastn": {
        "evalue": 1e-5,
        "outfmt": 6
    },
    "rpstblastn": {
        "evalue": 1e-5,
        "outfmt": 6
    },
    "lncRNA": {
        "min_length": 200
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
        logging.info(f"Generated config template at {config_path}. Please edit it and rerun.")
        exit(0)

def load_config(user_dir):
    config_path = os.path.join(user_dir, "annocript_config.yaml")
    generate_config_template(user_dir)
    with open(config_path) as f:
        return yaml.safe_load(f)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Annocript: Plant Transcriptome Annotation Tool with Checkpoints")
    parser.add_argument("--fasta", help="Input FASTA file", default=DEFAULT_CONFIG["fasta_file"])
    parser.add_argument("--db_dir", help="Database directory", default=DEFAULT_CONFIG["database"]["source_dir"])
    parser.add_argument("--threads", type=int, help="Number of CPU threads", default=8)
    parser.add_argument("--db_option", choices=["sprot", "trembl", "uniref90"], nargs="+", default=["sprot"],
                        help="DIAMOND databases: 'sprot' (Swiss-Prot), 'trembl' (TrEMBL), 'uniref90' (UniRef90).")
    parser.add_argument("--do_db_creation", action="store_true", help="Create database (skipped with pre-downloaded files)")
    parser.add_argument("--do_execute_programs", action="store_true", help="Execute alignment programs", default=True)
    parser.add_argument("--do_blastx", action="store_true", help="Run DIAMOND BLASTX", default=True)
    parser.add_argument("--do_blastn", action="store_true", help="Run BLASTN", default=True)
    parser.add_argument("--do_rpstblastn", action="store_true", help="Run RPSBLASTN", default=True)
    parser.add_argument("--do_lnc_prediction", action="store_true", help="Predict lncRNA with simple logic", default=True)
    parser.add_argument("--do_dna2pep", action="store_true", help="Search ORFs with TransDecoder", default=True)
    parser.add_argument("--do_build_output", action="store_true", help="Build final output", default=True)
    parser.add_argument("--extract_stats", action="store_true", help="Generate statistics", default=True)
    parser.add_argument("--do_kegg_annotation", action="store_true", help="Perform KEGG annotation", default=False)
    parser.add_argument("--do_extended_annotation", action="store_true", help="Perform Pfam, EC, TAIR, UniPathway annotations", default=False)
    parser.add_argument("--do_function_description", action="store_true", help="Add detailed gene function descriptions", default=False)
    return parser.parse_args()

def check_files(config, args):
    db_dir = os.path.abspath(args.db_dir)
    if not os.path.exists(args.fasta):
        raise FileNotFoundError(f"Input FASTA file not found: {args.fasta}")
    required_files = [
        config["database"]["uniprot_sprot"],
        config["database"]["uniref90"],
        config["database"]["uniprot_trembl"],
        config["database"]["ncRNA"],
        config["database"]["cdd"],
        config["database"]["go_obo"],
        config["database"]["enzyme_dat"],
        config["database"]["pfam2go"]
    ]
    if args.do_function_description:
        required_files.append(config["database"]["uniprot_dat"])
    if args.do_kegg_annotation and "pathway_file" in config["database"]:
        required_files.append(config["database"]["pathway_file"])
    for f in required_files:
        if not os.path.exists(os.path.join(db_dir, f)):
            raise FileNotFoundError(f"Missing {f} in {db_dir}")
    logging.info("All required database files found.")

def build_index(db_file, output_db, tool, threads):
    db_file = os.path.abspath(db_file)
    if not os.path.exists(db_file):
        raise FileNotFoundError(f"Database file not found: {db_file}")

    if db_file.endswith(".gz"):
        uncompressed = output_db + ".fasta"
        if not os.path.exists(uncompressed):
            with gzip.open(db_file, "rb") as f_in, open(uncompressed, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            logging.info(f"Uncompressed {db_file} to {uncompressed}")
        db_file = uncompressed

    index_exists = False
    if tool == "diamond":
        index_file = f"{output_db}.dmnd"
        index_exists = os.path.exists(index_file)
    elif tool == "blast":
        index_file = f"{output_db}.nhr"
        index_exists = os.path.exists(index_file) and os.path.exists(f"{output_db}.nin") and os.path.exists(f"{output_db}.nsq")

    if index_exists and os.path.getmtime(index_file) >= os.path.getmtime(db_file):
        logging.info(f"Index already exists and is up-to-date: {index_file}")
        return output_db

    if tool == "diamond":
        cmd = ["diamond", "makedb", "--in", db_file, "--db", output_db, "--threads", str(threads)]
        subprocess.run(cmd, check=True)
        logging.info(f"DIAMOND index built: {output_db}.dmnd")
    elif tool == "blast":
        cmd = ["makeblastdb", "-in", db_file, "-dbtype", "nucl", "-out", output_db]
        subprocess.run(cmd, check=True)
        logging.info(f"BLAST index built: {output_db}")
    return output_db

def run_alignment(args_tuple):
    tool, config, args, db_index, output_file, checkpoint_file, threads_per_task = args_tuple
    if os.path.exists(checkpoint_file):
        logging.info(f"{tool} alignment already completed: {output_file}")
        return output_file

    try:
        if tool == "diamond":
            cmd = [
                "diamond", "blastx",
                "--query", args.fasta,
                "--db", db_index,
                "--out", output_file,
                "--evalue", str(config["diamond"]["evalue"]),
                "--threads", str(threads_per_task),
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
        subprocess.run(cmd, check=True)
        with open(checkpoint_file, "w") as f:
            f.write("done")
        logging.info(f"{tool} alignment completed: {output_file}")
        return output_file
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"{tool} failed with error: {e}")

def predict_lncRNA(config, args, diamond_out_files, blastn_out):
    output_dir = config["output"]["dir"]
    lnc_file = os.path.join(output_dir, "lncRNA.fasta")
    checkpoint_file = os.path.join(output_dir, "lncRNA.done")

    if os.path.exists(checkpoint_file):
        logging.info(f"lncRNA prediction already completed: {lnc_file}")
        return lnc_file

    diamond_hits = set()
    for diamond_file in diamond_out_files:
        diamond_hits.update(pd.read_csv(diamond_file, sep="\t", header=None)[0])
    blastn_hits = set(pd.read_csv(blastn_out, sep="\t", header=None)[0])

    with open(lnc_file, "w") as out:
        for record in SeqIO.parse(args.fasta, "fasta"):
            if (record.id not in diamond_hits and
                record.id in blastn_hits and
                len(record.seq) >= config["lncRNA"]["min_length"]):
                SeqIO.write(record, out, "fasta")

    with open(checkpoint_file, "w") as f:
        f.write("done")
    logging.info(f"Simple lncRNA prediction completed: {lnc_file}")
    return lnc_file

def dna2pep(config, args):
    output_dir = config["output"]["dir"]
    transdecoder_dir = os.path.join(output_dir, "transdecoder")
    checkpoint_file = os.path.join(output_dir, "dna2pep.done")
    orf_file = f"{transdecoder_dir}/{os.path.basename(args.fasta)}.transdecoder.pep"

    if os.path.exists(checkpoint_file):
        logging.info(f"ORF search already completed: {orf_file}")
        return orf_file

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

    with open(checkpoint_file, "w") as f:
        f.write("done")
    logging.info(f"TransDecoder ORF search completed: {orf_file}")
    return orf_file

def prepare_idmapping(config, args):
    db_dir = args.db_dir
    idmapping_gz = os.path.join(db_dir, config["database"]["idmapping"])
    idmapping_parquet = os.path.join(db_dir, "idmapping.parquet")

    if os.path.exists(idmapping_parquet):
        logging.info(f"idmapping already prepared: {idmapping_parquet}")
        return idmapping_parquet

    df = pd.read_csv(idmapping_gz, sep="\t", header=None,
                     names=["uniprot_id", "type", "value"], compression="gzip")
    df.to_parquet(idmapping_parquet)
    logging.info(f"idmapping prepared: {idmapping_parquet}")
    return idmapping_parquet

def parse_uniprot_dat(dat_file):
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

def parse_go_obo(obo_file):
    go_terms = {}
    current_term = {}
    
    with open(obo_file, "r") as f:
        for line in f:
            line = line.strip()
            if line == "[Term]":
                if current_term.get("id") and current_term.get("name"):
                    go_terms[current_term["id"]] = current_term["name"]
                current_term = {}
            elif line.startswith("id:"):
                current_term["id"] = line.split("id: ")[1]
            elif line.startswith("name:"):
                current_term["name"] = line.split("name: ")[1]
    
    if current_term.get("id") and current_term.get("name"):
        go_terms[current_term["id"]] = current_term["name"]
    
    logging.info(f"Parsed {len(go_terms)} GO terms from {obo_file}")
    return go_terms

def parse_enzyme_dat(enzyme_file):
    enzyme_data = {}
    current_ec = None
    description = []
    
    with open(enzyme_file, "r") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("ID"):
                if current_ec and description:
                    enzyme_data[current_ec] = " ".join(description)
                current_ec = line.split()[1]
                description = []
            elif line.startswith("DE"):
                description.append(line[5:].strip("."))
            elif line.startswith("//"):
                if current_ec and description:
                    enzyme_data[current_ec] = " ".join(description)
                current_ec = None
                description = []
    
    logging.info(f"Parsed {len(enzyme_data)} enzyme entries from {enzyme_file}")
    return enzyme_data

def parse_cdd_metadata(cdd_tar_file, extract_dir):
    cdd_data = {}
    if not os.path.exists(extract_dir):
        with tarfile.open(cdd_tar_file, "r:gz") as tar:
            tar.extractall(extract_dir)
        logging.info(f"Extracted {cdd_tar_file} to {extract_dir}")
    
    for smp_file in os.listdir(extract_dir):
        if smp_file.endswith(".smp"):
            with open(os.path.join(extract_dir, smp_file), "r") as f:
                for line in f:
                    if line.startswith(">"):
                        parts = line[1:].strip().split()
                        domain_id = parts[0]  # e.g., cd00001
                        desc = " ".join(parts[1:])
                        cdd_data[domain_id] = desc
    logging.info(f"Parsed {len(cdd_data)} CDD entries from {extract_dir}")
    return cdd_data

def parse_pfam2go(pfam2go_file):
    pfam_to_go = {}
    with open(pfam2go_file, "r") as f:
        for line in f:
            if line.startswith("Pfam:"):
                parts = line.strip().split(" > ")
                pfam_part = parts[0].split()
                pfam_id = pfam_part[0].replace("Pfam:", "")  # e.g., PF00001
                go_terms = parts[1].split(" ; ")
                pfam_to_go[pfam_id] = [go.split()[-1] for go in go_terms]  # Extract GO IDs
    logging.info(f"Parsed {len(pfam_to_go)} Pfam-to-GO mappings from {pfam2go_file}")
    return pfam_to_go

def build_output(config, args, diamond_out_files, blastn_out, rpstblastn_out, idmapping_parquet):
    output_dir = config["output"]["dir"]
    annotations_file = os.path.join(output_dir, config["output"]["annotations"])
    gff_dir = config["output"]["gff_dir"]
    checkpoint_file = os.path.join(output_dir, "output.done")

    if os.path.exists(checkpoint_file):
        logging.info(f"Output already built: {annotations_file}")
        return annotations_file

    os.makedirs(gff_dir, exist_ok=True)

    conn = duckdb.connect()

    diamond_tables = []
    for i, diamond_file in enumerate(diamond_out_files):
        conn.execute(f"""
            CREATE TABLE diamond_{i} AS 
            SELECT * FROM read_csv_auto('{diamond_file}', delim='\t', 
                names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                       'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
        """)
        diamond_tables.append(f"diamond_{i}")

    diamond_query = " UNION ALL ".join(f"SELECT * FROM {table}" for table in diamond_tables)
    conn.execute(f"""
        CREATE TABLE diamond AS 
        SELECT qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
        FROM ({diamond_query})
        QUALIFY ROW_NUMBER() OVER (PARTITION BY qseqid ORDER BY bitscore DESC) = 1
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

    # Parse additional data sources
    go_obo_file = os.path.join(args.db_dir, config["database"]["go_obo"])
    go_terms = parse_go_obo(go_obo_file)
    
    enzyme_file = os.path.join(args.db_dir, config["database"]["enzyme_dat"])
    enzyme_data = parse_enzyme_dat(enzyme_file)
    
    cdd_dir = os.path.join(args.db_dir, "cdd")
    cdd_data = parse_cdd_metadata(os.path.join(args.db_dir, config["database"]["cdd"]), cdd_dir)
    
    pfam2go_file = os.path.join(args.db_dir, config["database"]["pfam2go"])
    pfam_to_go = parse_pfam2go(pfam2go_file)

    query = """
    SELECT d.qseqid, d.sseqid AS protein, d.evalue, d.bitscore,
           b.sseqid AS ncRNA, b.evalue AS nc_evalue,
           r.sseqid AS domain, r.evalue AS domain_evalue,
           GROUP_CONCAT(DISTINCT gm.value) AS go_ids,
           GROUP_CONCAT(DISTINCT em.value) AS ec_numbers,
           GROUP_CONCAT(DISTINCT pm.value) AS pfam_ids
    """
    if args.do_kegg_annotation or args.do_extended_annotation:
        query += """,
           GROUP_CONCAT(DISTINCT km.value) AS kegg_ids
        """
    if args.do_extended_annotation:
        query += """,
           GROUP_CONCAT(DISTINCT tm.value) AS tair_ids,
           GROUP_CONCAT(DISTINCT um.value) AS unipathway_ids
        """
    query += """
    FROM diamond d
    LEFT JOIN blastn b ON d.qseqid = b.qseqid
    LEFT JOIN rpstblastn r ON d.qseqid = r.qseqid
    LEFT JOIN idmapping gm ON d.sseqid = gm.uniprot_id AND gm.type = 'GO'
    LEFT JOIN idmapping em ON d.sseqid = em.uniprot_id AND em.type = 'EC'
    LEFT JOIN idmapping pm ON d.sseqid = pm.uniprot_id AND pm.type = 'Pfam'
    """
    if args.do_kegg_annotation or args.do_extended_annotation:
        query += """
        LEFT JOIN idmapping km ON d.sseqid = km.uniprot_id AND km.type = 'KEGG'
        """
    if args.do_extended_annotation:
        query += """
        LEFT JOIN idmapping tm ON d.sseqid = tm.uniprot_id AND tm.type = 'TAIR'
        LEFT JOIN idmapping um ON d.sseqid = um.uniprot_id AND um.type = 'UniPathway'
        """
    query += "GROUP BY ALL"

    result = conn.execute(query).fetchdf()

    # Map GO terms
    def map_go_terms(go_ids):
        if pd.isna(go_ids):
            return None
        return ";".join([f"{gid} ({go_terms.get(gid, 'Unknown')})" for gid in go_ids.split(",")])
    result["go_descriptions"] = result["go_ids"].apply(map_go_terms)

    # Map EC numbers to enzyme descriptions
    def map_enzyme(ec_numbers):
        if pd.isna(ec_numbers):
            return None
        return ";".join([f"{ec} ({enzyme_data.get(ec, 'Unknown')})" for ec in ec_numbers.split(",")])
    result["enzyme_descriptions"] = result["ec_numbers"].apply(map_enzyme)

    # Map CDD domains
    def map_cdd_domains(domain):
        if pd.isna(domain):
            return None
        return f"{domain} ({cdd_data.get(domain.split('.')[0], 'Unknown')})"
    result["domain_descriptions"] = result["domain"].apply(map_cdd_domains)

    # Add Pfam-to-GO inferred GO terms
    def infer_go_from_pfam(pfam_ids):
        if pd.isna(pfam_ids):
            return None
        inferred_go = set()
        for pfam in pfam_ids.split(","):
            inferred_go.update(pfam_to_go.get(pfam, []))
        return ";".join([f"{gid} ({go_terms.get(gid, 'Unknown')})" for gid in inferred_go]) if inferred_go else None
    result["pfam_inferred_go"] = result["pfam_ids"].apply(infer_go_from_pfam)

    if args.do_kegg_annotation and "pathway_file" in config["database"] and os.path.exists(os.path.join(args.db_dir, config["database"]["pathway_file"])):
        pathway_file = os.path.join(args.db_dir, config["database"]["pathway_file"])
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

    if args.do_function_description:
        dat_file = os.path.join(args.db_dir, config["database"]["uniprot_dat"])
        descriptions = parse_uniprot_dat(dat_file)
        result["function_description"] = result["protein"].map(
            lambda x: descriptions.get(x.split("|")[1] if "|" in x else x.split("_")[-1] if "UniRef90_" in x else x, "No description available")
        )

    result.to_csv(annotations_file, sep="\t", index=False)

    with open(os.path.join(gff_dir, "annotations.gff3"), "w") as gff:
        gff.write("##gff-version 3\n")
        for _, row in result.iterrows():
            if pd.notna(row["protein"]):
                gff.write(f"{row['qseqid']}\tDIAMOND\tmatch\t{row['qstart']}\t{row['qend']}\t"
                          f"{row['bitscore']}\t+\t.\tID={row['sseqid']};evalue={row['evalue']}\n")

    with open(checkpoint_file, "w") as f:
        f.write("done")
    conn.close()
    logging.info(f"Output generated: {annotations_file}, GFF3 in {gff_dir}")
    return annotations_file

def extract_statistics(config, args, annotations_file):
    stats_dir = config["output"]["stats_dir"]
    os.makedirs(stats_dir, exist_ok=True)
    checkpoint_file = os.path.join(stats_dir, "stats.done")

    if os.path.exists(checkpoint_file):
        logging.info(f"Statistics already generated: {stats_dir}/stats.txt")
        return

    df = pd.read_csv(annotations_file, sep="\t")

    stats = {
        "total_sequences": len(df["qseqid"].unique()),
        "with_protein": len(df[df["protein"].notna()]["qseqid"].unique()),
        "with_ncRNA": len(df[df["ncRNA"].notna()]["qseqid"].unique()),
        "with_domain": len(df[df["domain"].notna()]["qseqid"].unique()),
        "with_go": len(df[df["go_ids"].notna()]["qseqid"].unique()),
        "with_enzyme": len(df[df["ec_numbers"].notna()]["qseqid"].unique()),
        "with_pfam": len(df[df["pfam_ids"].notna()]["qseqid"].unique()),
        "with_pfam_inferred_go": len(df[df["pfam_inferred_go"].notna()]["qseqid"].unique()),
        "with_kegg": len(df[df["kegg_ids"].notna()]["qseqid"].unique()) if "kegg_ids" in df else 0,
        "with_tair": len(df[df["tair_ids"].notna()]["qseqid"].unique()) if "tair_ids" in df else 0,
        "with_unipathway": len(df[df["unipathway_ids"].notna()]["qseqid"].unique()) if "unipathway_ids" in df else 0,
        "with_function": len(df[df["function_description"].notna()]["qseqid"].unique()) if "function_description" in df else 0
    }
    with open(os.path.join(stats_dir, "stats.txt"), "w") as f:
        yaml.dump(stats, f)

    with open(checkpoint_file, "w") as f:
        f.write("done")
    logging.info(f"Statistics saved to {stats_dir}/stats.txt")

def main():
    args = parse_arguments()
    config = load_config(os.getcwd())

    check_files(config, args)

    if args.do_execute_programs:
        diamond_out = blastn_out = rpstblastn_out = None
        output_dir = config["output"]["dir"]
        os.makedirs(output_dir, exist_ok=True)

        tasks = []
        db_map = {"sprot": "uniprot_sprot", "trembl": "uniprot_trembl", "uniref90": "uniref90"}
        diamond_out_files = []

        total_tasks = sum([args.do_blastx * len(args.db_option), args.do_blastn, args.do_rpstblastn])
        threads_per_task = max(1, args.threads // max(1, total_tasks))

        if args.do_blastx:
            for db_option in args.db_option:
                db_name = db_map[db_option]
                db_index = build_index(
                    os.path.join(args.db_dir, config["database"][db_name]),
                    os.path.join(args.db_dir, db_name),
                    "diamond",
                    threads_per_task
                )
                diamond_out_file = os.path.join(output_dir, f"{db_option}_{config['output']['diamond_out']}")
                tasks.append(("diamond", config, args, db_index, diamond_out_file,
                              os.path.join(output_dir, f"{db_option}_diamond.done"), threads_per_task))
                diamond_out_files.append(diamond_out_file)

        if args.do_blastn:
            nc_index = build_index(os.path.join(args.db_dir, config["database"]["ncRNA"]),
                                   os.path.join(args.db_dir, "ncRNA"), "blast", threads_per_task)
            blastn_out_file = os.path.join(output_dir, config["output"]["blastn_out"])
            tasks.append(("blastn", config, args, nc_index, blastn_out_file,
                          os.path.join(output_dir, "blastn.done"), threads_per_task))
            blastn_out = blastn_out_file

        if args.do_rpstblastn:
            cdd_dir = os.path.join(args.db_dir, "cdd")
            if not os.path.exists(cdd_dir):
                os.makedirs(cdd_dir)
                subprocess.run(["tar", "-xzf", os.path.join(args.db_dir, config["database"]["cdd"]), "-C", cdd_dir], check=True)
            rpstblastn_out_file = os.path.join(output_dir, config["output"]["rpstblastn_out"])
            tasks.append(("rpstblastn", config, args, os.path.join(cdd_dir, "Cdd"), rpstblastn_out_file,
                          os.path.join(output_dir, "rpstblastn.done"), threads_per_task))
            rpstblastn_out = rpstblastn_out_file

        with Pool(processes=min(len(tasks), args.threads)) as pool:
            results = pool.map(run_alignment, tasks)

        diamond_out = [r for r in results if "diamond" in r] or diamond_out_files
        blastn_out = blastn_out or next((r for r in results if "blastn" in r), None)
        rpstblastn_out = rpstblastn_out or next((r for r in results if "rpstblastn" in r), None)

        if args.do_lnc_prediction and diamond_out and blastn_out:
            predict_lncRNA(config, args, diamond_out, blastn_out)
        if args.do_dna2pep:
            dna2pep(config, args)

    annotations_file = None
    if args.do_build_output:
        idmapping_parquet = prepare_idmapping(config, args)
        annotations_file = build_output(config, args, diamond_out, blastn_out, rpstblastn_out, idmapping_parquet)

    if args.extract_stats and annotations_file:
        extract_statistics(config, args, annotations_file)

if __name__ == "__main__":
    main()
