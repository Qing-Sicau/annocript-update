import os
import subprocess
import argparse
import yaml
import pandas as pd
import dask.dataframe as dd
import gzip
import shutil
import tarfile
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import logging
import glob
import duckdb
import pickle
import hashlib
from pathlib import Path
import pyfaidx
import time
import pyarrow as pa
import pyarrow.parquet as pq

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
    "duckdb": {
        "path": "annocript.duckdb"
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
        "stats_dir": "./stats",
        "cache_dir": "./cache"
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
    parser = argparse.ArgumentParser(description="Optimized Annocript: Plant Transcriptome Annotation Tool")
    parser.add_argument("--fasta", help="Input FASTA file", default=DEFAULT_CONFIG["fasta_file"])
    parser.add_argument("--db_dir", help="Database directory", default=DEFAULT_CONFIG["database"]["source_dir"])
    parser.add_argument("--threads", type=int, help="Number of CPU threads", default=8)
    parser.add_argument("--db_option", choices=["sprot", "trembl", "uniref90"], nargs="+", default=["sprot"])
    parser.add_argument("--do_db_creation", action="store_true")
    parser.add_argument("--do_execute_programs", action="store_true", default=True)
    parser.add_argument("--do_blastx", action="store_true", default=True)
    parser.add_argument("--do_blastn", action="store_true", default=True)
    parser.add_argument("--do_rpstblastn", action="store_true", default=True)
    parser.add_argument("--do_lnc_prediction", action="store_true", default=True)
    parser.add_argument("--do_dna2pep", action="store_true", default=True)
    parser.add_argument("--do_build_output", action="store_true", default=True)
    parser.add_argument("--extract_stats", action="store_true", default=True)
    parser.add_argument("--do_kegg_annotation", action="store_true", default=False)
    parser.add_argument("--do_extended_annotation", action="store_true", default=False)
    parser.add_argument("--do_function_description", action="store_true", default=False)
    return parser.parse_args()

def pre_uncompress_if_needed(gz_file, output_file):
    if not os.path.exists(output_file) or os.path.getmtime(gz_file) > os.path.getmtime(output_file):
        with gzip.open(gz_file, "rb") as f_in, open(output_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        logging.info(f"Pre-uncompressed {gz_file} to {output_file}")
    return output_file

def check_files(config, args):
    db_dir = os.path.abspath(args.db_dir)
    if not os.path.exists(args.fasta):
        raise FileNotFoundError(f"Input FASTA file not found: {args.fasta}")
    required_files = ["uniprot_sprot", "uniref90", "uniprot_trembl", "ncRNA", "cdd", "go_obo", "enzyme_dat", "pfam2go"]
    if args.do_function_description:
        required_files.append("uniprot_dat")
    if args.do_kegg_annotation and "pathway_file" in config["database"]:
        required_files.append("pathway_file")
    missing_files = [f for f in required_files if not os.path.exists(os.path.join(db_dir, config["database"][f]))]
    if missing_files:
        raise FileNotFoundError(f"Missing files in {db_dir}: {', '.join(missing_files)}")
    logging.info("All required database files found.")

def build_index(db_file, output_db, tool, threads):
    """Build index for alignment tools."""
    db_file = os.path.abspath(db_file)
    output_db = os.path.abspath(output_db)

    if not os.path.exists(db_file):
        raise FileNotFoundError(f"Database file not found: {db_file}")

    if db_file.endswith(".gz"):
        uncompressed = output_db + ".fasta"
        if not os.path.exists(uncompressed):
            with gzip.open(db_file, "rb") as f_in, open(uncompressed, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            logging.info(f"Uncompressed {db_file} to {uncompressed}")
        db_file = uncompressed

    if tool == "diamond":
        index_file = f"{output_db}.dmnd"
    elif tool == "blast":
        index_file = f"{output_db}.00.nhr"
    else:
        raise ValueError(f"Unsupported tool: {tool}")

    logging.info(f"Checking index file: {index_file}, exists: {os.path.exists(index_file)}")

    if os.path.exists(index_file):
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
    tool, config, args, db_index, output_file, checkpoint_file, threads = args_tuple
    if os.path.exists(checkpoint_file):
        logging.info(f"{tool} alignment already completed: {output_file}")
        return output_file

    cmd_base = {
        "diamond": ["diamond", "blastx", "--query", args.fasta, "--db", db_index, "--out", output_file,
                    "--evalue", str(config["diamond"]["evalue"]), "--threads", str(threads),
                    "--outfmt", str(config["diamond"]["outfmt"])] + (["--sensitive"] if config["diamond"]["sensitive"] else []),
        "blastn": ["blastn", "-query", args.fasta, "-db", db_index, "-out", output_file,
                   "-evalue", str(config["blastn"]["evalue"]), "-num_threads", str(threads),
                   "-outfmt", str(config["blastn"]["outfmt"])],
        "rpstblastn": ["rpstblastn", "-query", args.fasta, "-db", db_index, "-out", output_file,
                       "-evalue", str(config["rpstblastn"]["evalue"]), "-num_threads", str(threads),
                       "-outfmt", str(config["rpstblastn"]["outfmt"])]
    }
    subprocess.run(cmd_base[tool], check=True)
    with open(checkpoint_file, "w") as f:
        f.write("done")
    logging.info(f"{tool} alignment completed: {output_file}")
    return output_file

def predict_lncRNA(config, args, diamond_out_files, blastn_out):
    output_dir = config["output"]["dir"]
    lnc_file = os.path.join(output_dir, "lncRNA.fasta")
    checkpoint_file = os.path.join(output_dir, "lncRNA.done")
    if os.path.exists(checkpoint_file):
        logging.info(f"lncRNA prediction already completed: {lnc_file}")
        return lnc_file

    diamond_hits = set()
    for diamond_file in diamond_out_files:
        df = dd.read_csv(diamond_file, sep="\t", header=None, usecols=[0])
        diamond_hits.update(df[0].compute())
    blastn_df = dd.read_csv(blastn_out, sep="\t", header=None, usecols=[0])
    blastn_hits = set(blastn_df[0].compute())

    fasta = pyfaidx.Fasta(args.fasta)
    with open(lnc_file, "w") as out:
        for seq_id in fasta.keys():
            if (seq_id not in diamond_hits and seq_id in blastn_hits and
                len(fasta[seq_id]) >= config["lncRNA"]["min_length"]):
                out.write(f">{seq_id}\n{str(fasta[seq_id])}\n")
    with open(checkpoint_file, "w") as f:
        f.write("done")
    logging.info(f"lncRNA prediction completed: {lnc_file}")
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
    subprocess.run(["TransDecoder.LongOrfs", "-t", args.fasta, "-m", str(config["transdecoder"]["min_length"]),
                    "--output_dir", transdecoder_dir], check=True)
    subprocess.run(["TransDecoder.Predict", "-t", args.fasta, "--output_dir", transdecoder_dir, "--single_best_only"],
                   check=True)
    with open(checkpoint_file, "w") as f:
        f.write("done")
    logging.info(f"TransDecoder ORF search completed: {orf_file}")
    return orf_file

def prepare_idmapping(config, args):
    db_dir = args.db_dir
    idmapping_gz = os.path.join(db_dir, config["database"]["idmapping"])
    idmapping_parquet = os.path.join(db_dir, "idmapping.parquet")
    if os.path.exists(idmapping_parquet):
        logging.info(f"ID mapping already prepared: {idmapping_parquet}")
        return idmapping_parquet

    chunk_size = 10_000_000
    schema = None
    writer = None

    # Define all 22 columns from UniProt idmapping_selected.tab.gz
    column_names = [
        "UniProtKB_AC", "UniProtKB_ID", "GeneID", "RefSeq", "GI", "PDB", "GO", "UniRef100",
        "UniRef90", "UniRef50", "UniParc", "PIR", "NCBI_taxon", "MIM", "UniGene", "PubMed",
        "EMBL", "EMBL_CDS", "Ensembl", "Ensembl_TRS", "Ensembl_PRO", "Additional_PubMed"
    ]

    try:
        for i, chunk in enumerate(pd.read_csv(
            idmapping_gz,
            sep="\t",
            header=None,
            names=column_names,
            dtype=str,
            compression="gzip",
            chunksize=chunk_size
        )):
            table = pa.Table.from_pandas(chunk)
            if i == 0:
                schema = table.schema
                writer = pq.ParquetWriter(idmapping_parquet, schema, compression="snappy")
                writer.write_table(table)
            else:
                table = table.cast(schema)
                writer.write_table(table)
            logging.info(f"Processed chunk {i + 1} with {len(chunk)} rows")
    finally:
        if writer:
            writer.close()
    logging.info(f"ID mapping prepared: {idmapping_parquet}")
    return idmapping_parquet

def cached_parse(file_path, parse_func, cache_file):
    if os.path.exists(cache_file) and os.path.getmtime(cache_file) > os.path.getmtime(file_path):
        with open(cache_file, "rb") as f:
            return pickle.load(f)
    result = parse_func(file_path)
    with open(cache_file, "wb") as f:
        pickle.dump(result, f)
    return result

def parse_uniprot_dat(dat_file):
    def parse(dat_file):
        descriptions = {}
        ec_numbers = {}
        with gzip.open(dat_file, "rt") as f:
            uniprot_id, desc, ec = None, [], []
            for line in f:
                if line.startswith("ID"):
                    uniprot_id = line.split()[1]
                elif line.startswith("DE   RecName: Full="):
                    desc.append(line.split("=", 1)[1].strip(";"))
                elif line.startswith("CC   -!- FUNCTION:"):
                    desc.append(line.split(":", 1)[1].strip())
                elif line.startswith("DR   EC;"):
                    ec.append(line.split(";")[1].strip())
                elif line.startswith("//"):
                    if uniprot_id:
                        if desc:
                            descriptions[uniprot_id] = " ".join(desc)
                        if ec:
                            ec_numbers[uniprot_id] = ";".join(ec)
                    uniprot_id, desc, ec = None, [], []
        return descriptions, ec_numbers
    cache_file = dat_file.replace(".dat.gz", ".pickle")
    if os.path.exists(cache_file) and os.path.getmtime(cache_file) > os.path.getmtime(dat_file):
        with open(cache_file, "rb") as f:
            return pickle.load(f)
    result = parse(dat_file)
    with open(cache_file, "wb") as f:
        pickle.dump(result, f)
    return result

def parse_go_obo(obo_file):
    def parse(obo_file):
        go_terms = {}
        with open(obo_file) as f:
            current_term = {}
            for line in f:
                line = line.strip()
                if line == "[Term]":
                    if "id" in current_term and "name" in current_term:
                        go_terms[current_term["id"]] = current_term["name"]
                    current_term = {}
                elif line.startswith("id:"):
                    current_term["id"] = line.split("id: ")[1]
                elif line.startswith("name:"):
                    current_term["name"] = line.split("name: ")[1]
        if "id" in current_term and "name" in current_term:
            go_terms[current_term["id"]] = current_term["name"]
        return go_terms
    return cached_parse(obo_file, parse, obo_file + ".pickle")

def parse_enzyme_dat(enzyme_file):
    def parse(enzyme_file):
        enzyme_data = {}
        with open(enzyme_file) as f:
            current_ec, description = None, []
            for line in f:
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
                    current_ec, description = None, []
        return enzyme_data
    return cached_parse(enzyme_file, parse, enzyme_file + ".pickle")

def parse_cdd_metadata(cdd_tar_file, extract_dir):
    def parse(extract_dir):
        cdd_data = {}
        for smp_file in glob.glob(os.path.join(extract_dir, "*.smp")):
            with open(smp_file) as f:
                for line in f:
                    if line.startswith(">"):
                        parts = line[1:].strip().split()
                        cdd_data[parts[0]] = " ".join(parts[1:])
        return cdd_data
    cache_file = os.path.join(extract_dir, "cdd.pickle")
    if not os.path.exists(extract_dir):
        with tarfile.open(cdd_tar_file, "r:gz") as tar:
            tar.extractall(extract_dir)
        logging.info(f"Extracted {cdd_tar_file} to {extract_dir}")
    return cached_parse(extract_dir, parse, cache_file)

def parse_pfam2go(pfam2go_file):
    def parse(pfam2go_file):
        pfam_to_go = {}
        with open(pfam2go_file) as f:
            for line in f:
                if line.startswith("Pfam:"):
                    parts = line.strip().split(" > ")
                    pfam_id = parts[0].split()[0].replace("Pfam:", "")
                    go_terms = parts[1].split(" ; ")
                    pfam_to_go[pfam_id] = [go.split()[-1] for go in go_terms]
        return pfam_to_go
    return cached_parse(pfam2go_file, parse, pfam2go_file + ".pickle")

def build_output(config, args, diamond_out_files, blastn_out, rpstblastn_out, idmapping_parquet):
    output_dir = config["output"]["dir"]
    annotations_file = os.path.join(output_dir, config["output"]["annotations"])
    gff_dir = config["output"]["gff_dir"]
    checkpoint_file = os.path.join(output_dir, "output.done")
    if os.path.exists(checkpoint_file):
        logging.info(f"Output already built: {annotations_file}")
        return annotations_file

    os.makedirs(gff_dir, exist_ok=True)
    conn = duckdb.connect(config["duckdb"]["path"])
    conn.execute("CREATE TABLE IF NOT EXISTS diamond (qseqid VARCHAR, sseqid VARCHAR, pident FLOAT, length INT, mismatch INT, gapopen INT, qstart INT, qend INT, sstart INT, send INT, evalue DOUBLE, bitscore FLOAT)")
    conn.execute("CREATE TABLE IF NOT EXISTS blastn (qseqid VARCHAR, sseqid VARCHAR, pident FLOAT, length INT, mismatch INT, gapopen INT, qstart INT, qend INT, sstart INT, send INT, evalue DOUBLE, bitscore FLOAT)")
    conn.execute("CREATE TABLE IF NOT EXISTS rpstblastn (qseqid VARCHAR, sseqid VARCHAR, pident FLOAT, length INT, mismatch INT, gapopen INT, qstart INT, qend INT, sstart INT, send INT, evalue DOUBLE, bitscore FLOAT)")
    
    # Create idmapping table with all 22 columns
    conn.execute("""
        CREATE TABLE IF NOT EXISTS idmapping (
            UniProtKB_AC VARCHAR, UniProtKB_ID VARCHAR, GeneID VARCHAR, RefSeq VARCHAR, GI VARCHAR,
            PDB VARCHAR, GO VARCHAR, UniRef100 VARCHAR, UniRef90 VARCHAR, UniRef50 VARCHAR,
            UniParc VARCHAR, PIR VARCHAR, NCBI_taxon VARCHAR, MIM VARCHAR, UniGene VARCHAR,
            PubMed VARCHAR, EMBL VARCHAR, EMBL_CDS VARCHAR, Ensembl VARCHAR, Ensembl_TRS VARCHAR,
            Ensembl_PRO VARCHAR, Additional_PubMed VARCHAR
        )
    """)

    for diamond_file in diamond_out_files:
        conn.execute(f"INSERT INTO diamond SELECT * FROM read_csv_auto('{diamond_file}', delim='\t', header=false, columns={{'qseqid': 'VARCHAR', 'sseqid': 'VARCHAR', 'pident': 'FLOAT', 'length': 'INT', 'mismatch': 'INT', 'gapopen': 'INT', 'qstart': 'INT', 'qend': 'INT', 'sstart': 'INT', 'send': 'INT', 'evalue': 'DOUBLE', 'bitscore': 'FLOAT'}})")
    if blastn_out:
        conn.execute(f"INSERT INTO blastn SELECT * FROM read_csv_auto('{blastn_out}', delim='\t', header=false, columns={{'qseqid': 'VARCHAR', 'sseqid': 'VARCHAR', 'pident': 'FLOAT', 'length': 'INT', 'mismatch': 'INT', 'gapopen': 'INT', 'qstart': 'INT', 'qend': 'INT', 'sstart': 'INT', 'send': 'INT', 'evalue': 'DOUBLE', 'bitscore': 'FLOAT'}})")
    if rpstblastn_out:
        conn.execute(f"INSERT INTO rpstblastn SELECT * FROM read_csv_auto('{rpstblastn_out}', delim='\t', header=false, columns={{'qseqid': 'VARCHAR', 'sseqid': 'VARCHAR', 'pident': 'FLOAT', 'length': 'INT', 'mismatch': 'INT', 'gapopen': 'INT', 'qstart': 'INT', 'qend': 'INT', 'sstart': 'INT', 'send': 'INT', 'evalue': 'DOUBLE', 'bitscore': 'FLOAT'}})")
    conn.execute(f"INSERT INTO idmapping SELECT * FROM parquet_scan('{idmapping_parquet}')")

    conn.execute("""
        CREATE TABLE diamond_best AS
        SELECT d.* FROM diamond d
        INNER JOIN (SELECT qseqid, MAX(bitscore) AS max_bitscore FROM diamond GROUP BY qseqid) db
        ON d.qseqid = db.qseqid AND d.bitscore = db.max_bitscore
    """)

    query = """
    SELECT d.qseqid, d.sseqid AS protein, d.evalue, d.bitscore, b.sseqid AS ncRNA, b.evalue AS nc_evalue,
           r.sseqid AS domain, r.evalue AS domain_evalue,
           im.GO AS go_ids,
           im.Pfam AS pfam_ids
    """
    if args.do_kegg_annotation or args.do_extended_annotation:
        query += ", im.KEGG AS kegg_ids"
    if args.do_extended_annotation:
        query += ", im.TAIR AS tair_ids, im.UniPathway AS unipathway_ids"
    query += """
    FROM diamond_best d
    LEFT JOIN blastn b ON d.qseqid = b.qseqid
    LEFT JOIN rpstblastn r ON d.qseqid = r.qseqid
    LEFT JOIN idmapping im ON d.sseqid = im.UniProtKB_AC
    """

    result = conn.execute(query).fetchdf()

    with ThreadPoolExecutor() as executor:
        go_terms = executor.submit(parse_go_obo, os.path.join(args.db_dir, config["database"]["go_obo"])).result()
        enzyme_data = executor.submit(parse_enzyme_dat, os.path.join(args.db_dir, config["database"]["enzyme_dat"])).result()
        cdd_data = executor.submit(parse_cdd_metadata, os.path.join(args.db_dir, config["database"]["cdd"]), os.path.join(args.db_dir, "cdd")).result()
        pfam_to_go = executor.submit(parse_pfam2go, os.path.join(args.db_dir, config["database"]["pfam2go"])).result()
        uniprot_data = executor.submit(parse_uniprot_dat, os.path.join(args.db_dir, config["database"]["uniprot_dat"])).result()
        descriptions, ec_numbers = uniprot_data

    # Add EC numbers and descriptions
    result["ec_numbers"] = result["protein"].apply(lambda x: ec_numbers.get(x.split("|")[1] if "|" in x else x, None))
    result["go_descriptions"] = result["go_ids"].apply(lambda x: map_go_terms(x, go_terms) if pd.notna(x) else None)
    result["enzyme_descriptions"] = result["ec_numbers"].apply(lambda x: map_enzyme(x, enzyme_data) if pd.notna(x) else None)
    result["domain_descriptions"] = result["domain"].apply(lambda x: map_cdd_domains(x, cdd_data) if pd.notna(x) else None)
    result["pfam_inferred_go"] = result["pfam_ids"].apply(lambda x: infer_go_from_pfam(x, pfam_to_go, go_terms) if pd.notna(x) else None)

    if args.do_kegg_annotation and "pathway_file" in config["database"]:
        pathways = {line.split("\t")[0].replace("PATH:", ""): line.split("\t")[1].strip()
                    for line in open(os.path.join(args.db_dir, config["database"]["pathway_file"])) if line.startswith("PATH:")}
        result["kegg_pathways"] = result["kegg_ids"].apply(lambda x: map_pathways(x, pathways) if pd.notna(x) else None)

    if args.do_function_description:
        result["function_description"] = result["protein"].apply(lambda x: descriptions.get(x.split("|")[1] if "|" in x else x, "No description"))

    result.to_csv(annotations_file, sep="\t", index=False)

    with open(os.path.join(gff_dir, "annotations.gff3"), "w") as gff:
        gff.write("##gff-version 3\n")
        for row in conn.execute("SELECT qseqid, sseqid, qstart, qend, bitscore, evalue FROM diamond_best WHERE sseqid IS NOT NULL").fetchall():
            gff.write(f"{row[0]}\tDIAMOND\tmatch\t{row[2]}\t{row[3]}\t{row[4]}\t+\t.\tID={row[1]};evalue={row[5]}\n")

    with open(checkpoint_file, "w") as f:
        f.write("done")
    conn.close()
    logging.info(f"Output generated: {annotations_file}, GFF3 in {gff_dir}")
    return annotations_file

def map_go_terms(go_ids, go_terms):
    return ";".join([f"{gid} ({go_terms.get(gid, 'Unknown')})" for gid in go_ids.split(";")])

def map_enzyme(ec_numbers, enzyme_data):
    if not ec_numbers:
        return None
    return ";".join([f"{ec} ({enzyme_data.get(ec, 'Unknown')})" for ec in ec_numbers.split(";")])

def map_cdd_domains(domain, cdd_data):
    return f"{domain} ({cdd_data.get(domain.split('.')[0], 'Unknown')})"

def infer_go_from_pfam(pfam_ids, pfam_to_go, go_terms):
    inferred_go = set()
    for pfam in pfam_ids.split(";"):
        inferred_go.update(pfam_to_go.get(pfam, []))
    return ";".join([f"{gid} ({go_terms.get(gid, 'Unknown')})" for gid in inferred_go]) if inferred_go else None

def map_pathways(kegg_ids, pathways):
    return ";".join([f"{kid} ({pathways.get(kid.split(':')[-1], 'Unknown')})" for kid in kegg_ids.split(";")])

def extract_statistics(config, args, annotations_file):
    stats_dir = config["output"]["stats_dir"]
    os.makedirs(stats_dir, exist_ok=True)
    checkpoint_file = os.path.join(stats_dir, "stats.done")
    if os.path.exists(checkpoint_file):
        logging.info(f"Statistics already generated: {stats_dir}/stats.txt")
        return

    df = dd.read_csv(annotations_file, sep="\t")
    stats = {
        "total_sequences": len(df["qseqid"].unique().compute()),
        "with_protein": len(df[df["protein"].notna()]["qseqid"].unique().compute()),
        "with_ncRNA": len(df[df["ncRNA"].notna()]["qseqid"].unique().compute()),
        "with_domain": len(df[df["domain"].notna()]["qseqid"].unique().compute()),
        "with_go": len(df[df["go_ids"].notna()]["qseqid"].unique().compute()),
        "with_enzyme": len(df[df["ec_numbers"].notna()]["qseqid"].unique().compute()),
        "with_pfam": len(df[df["pfam_ids"].notna()]["qseqid"].unique().compute()),
        "with_pfam_inferred_go": len(df[df["pfam_inferred_go"].notna()]["qseqid"].unique().compute()),
        "with_kegg": len(df[df["kegg_ids"].notna()]["qseqid"].unique().compute()) if "kegg_ids" in df else 0,
        "with_tair": len(df[df["tair_ids"].notna()]["qseqid"].unique().compute()) if "tair_ids" in df else 0,
        "with_unipathway": len(df[df["unipathway_ids"].notna()]["qseqid"].unique().compute()) if "unipathway_ids" in df else 0,
        "with_function": len(df[df["function_description"].notna()]["qseqid"].unique().compute()) if "function_description" in df else 0
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

    output_dir = config["output"]["dir"]
    os.makedirs(output_dir, exist_ok=True)
    diamond_out_files, blastn_out, rpstblastn_out = [], None, None

    if args.do_execute_programs:
        tasks = []
        db_map = {"sprot": "uniprot_sprot", "trembl": "uniprot_trembl", "uniref90": "uniref90"}
        total_tasks = sum([args.do_blastx * len(args.db_option), args.do_blastn, args.do_rpstblastn])
        threads_per_task = max(1, args.threads // max(1, total_tasks))

        if args.do_blastx:
            for db_option in args.db_option:
                db_index = build_index(os.path.join(args.db_dir, config["database"][db_map[db_option]]),
                                      os.path.join(args.db_dir, db_option), "diamond", threads_per_task)
                diamond_out_file = os.path.join(output_dir, f"{db_option}_{config['output']['diamond_out']}")
                tasks.append(("diamond", config, args, db_index, diamond_out_file,
                              os.path.join(output_dir, f"{db_option}_diamond.done"), threads_per_task))
                diamond_out_files.append(diamond_out_file)

        if args.do_blastn:
            nc_index = build_index(os.path.join(args.db_dir, config["database"]["ncRNA"]),
                                   os.path.join(args.db_dir, "ncRNA"), "blast", threads_per_task)
            blastn_out = os.path.join(output_dir, config["output"]["blastn_out"])
            tasks.append(("blastn", config, args, nc_index, blastn_out,
                          os.path.join(output_dir, "blastn.done"), threads_per_task))

        if args.do_rpstblastn:
            cdd_dir = os.path.join(args.db_dir, "cdd")
            if not os.path.exists(cdd_dir):
                os.makedirs(cdd_dir)
                subprocess.run(["tar", "-xzf", os.path.join(args.db_dir, config["database"]["cdd"]), "-C", cdd_dir], check=True)
            rpstblastn_out = os.path.join(output_dir, config["output"]["rpstblastn_out"])
            tasks.append(("rpstblastn", config, args, os.path.join(cdd_dir, "Cdd"), rpstblastn_out,
                          os.path.join(output_dir, "rpstblastn.done"), threads_per_task))

        if tasks:
            with ProcessPoolExecutor(max_workers=min(len(tasks), args.threads)) as pool:
                results = pool.map(run_alignment, tasks)
            diamond_out_files = [r for r in results if "diamond" in r] or diamond_out_files
            blastn_out = next((r for r in results if "blastn" in r), blastn_out)
            rpstblastn_out = next((r for r in results if "rpstblastn" in r), rpstblastn_out)

        if args.do_lnc_prediction and diamond_out_files and blastn_out:
            predict_lncRNA(config, args, diamond_out_files, blastn_out)
        if args.do_dna2pep:
            dna2pep(config, args)

    if args.do_build_output or args.extract_stats:
        if not diamond_out_files:
            diamond_out_files = glob.glob(os.path.join(output_dir, "*_diamond.out"))
            if not diamond_out_files:
                raise FileNotFoundError(f"No DIAMOND output files found in {output_dir}")
            logging.info(f"Detected DIAMOND output files: {diamond_out_files}")

        if not blastn_out and os.path.exists(os.path.join(output_dir, config["output"]["blastn_out"])):
            blastn_out = os.path.join(output_dir, config["output"]["blastn_out"])
            logging.info(f"Detected BLASTN output file: {blastn_out}")

        if not rpstblastn_out and os.path.exists(os.path.join(output_dir, config["output"]["rpstblastn_out"])):
            rpstblastn_out = os.path.join(output_dir, config["output"]["rpstblastn_out"])
            logging.info(f"Detected RPS-BLASTN output file: {rpstblastn_out}")

        if args.do_build_output:
            idmapping_parquet = prepare_idmapping(config, args)
            annotations_file = build_output(config, args, diamond_out_files, blastn_out, rpstblastn_out, idmapping_parquet)
        else:
            annotations_file = os.path.join(output_dir, config["output"]["annotations"])

        if args.extract_stats:
            if not os.path.exists(annotations_file):
                raise FileNotFoundError(f"Annotations file not found: {annotations_file}")
            extract_statistics(config, args, annotations_file)

if __name__ == "__main__":
    main()
