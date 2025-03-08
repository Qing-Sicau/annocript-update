import os
import subprocess
import argparse
import yaml
import pandas as pd
import gzip
import shutil
import tarfile
from Bio import SeqIO
from multiprocessing import Pool
import logging
import glob
import mysql.connector

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Default configuration with MySQL section added
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
    "mysql": {
        "host": "localhost",
        "port": 3306,
        "user": "pasa_write",
        "password": "pasa_write_pwd",
        "database": "annocript_db"
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
    """Generate a default config file if it doesn't exist."""
    config_path = os.path.join(user_dir, "annocript_config.yaml")
    if not os.path.exists(config_path):
        with open(config_path, "w") as f:
            yaml.dump(DEFAULT_CONFIG, f, default_flow_style=False)
        logging.info(f"Generated config template at {config_path}. Please edit it and rerun.")
        exit(0)

def load_config(user_dir):
    """Load the configuration file, generating it if missing."""
    config_path = os.path.join(user_dir, "annocript_config.yaml")
    generate_config_template(user_dir)
    with open(config_path) as f:
        return yaml.safe_load(f)

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Annocript: Plant Transcriptome Annotation Tool with Checkpoints")
    parser.add_argument("--fasta", help="Input FASTA file", default=DEFAULT_CONFIG["fasta_file"])
    parser.add_argument("--db_dir", help="Database directory", default=DEFAULT_CONFIG["database"]["source_dir"])
    parser.add_argument("--threads", type=int, help="Number of CPU threads", default=8)
    parser.add_argument("--db_option", choices=["sprot", "trembl", "uniref90"], nargs="+", default=["sprot"],
                        help="DIAMOND databases: 'sprot' (Swiss-Prot), 'trembl' (TrEMBL), 'uniref90' (UniRef90).")
    parser.add_argument("--do_db_creation", action="store_true", help="Create database (skipped with pre-downloaded files)")
    parser.add_argument("--do_execute_programs", action="store_true", help="Execute alignment programs", default=True)
    parser.add_argument("--do_blastx", action="store_true", help="Run DIAMOND BLASTX", default=True)
    parser.add_argument("--do_blastn", action="store_true", help="Run BLASTN against ncRNA database", default=True)
    parser.add_argument("--do_rpstblastn", action="store_true", help="Run RPSBLASTN against CDD domains", default=True)
    parser.add_argument("--do_lnc_prediction", action="store_true", help="Predict lncRNA", default=True)
    parser.add_argument("--do_dna2pep", action="store_true", help="Search ORFs with TransDecoder", default=True)
    parser.add_argument("--do_build_output", action="store_true", help="Build final output", default=True)
    parser.add_argument("--extract_stats", action="store_true", help="Generate statistics", default=True)
    parser.add_argument("--do_kegg_annotation", action="store_true", help="Perform KEGG annotation", default=False)
    parser.add_argument("--do_extended_annotation", action="store_true", help="Perform Pfam, EC, TAIR, UniPathway annotations", default=False)
    parser.add_argument("--do_function_description", action="store_true", help="Add detailed gene function descriptions", default=False)
    return parser.parse_args()

def check_files(config, args):
    """Check if required files exist."""
    db_dir = os.path.abspath(args.db_dir)
    if not os.path.exists(args.fasta):
        raise FileNotFoundError(f"Input FASTA file not found: {args.fasta}")

    if args.do_execute_programs:
        required_files = [
            "uniprot_sprot", "uniref90", "uniprot_trembl", "ncRNA", "cdd",
            "go_obo", "enzyme_dat", "pfam2go"
        ]
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
    """Run alignment tools with checkpointing."""
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

    try:
        subprocess.run(cmd_base[tool], check=True)
        with open(checkpoint_file, "w") as f:
            f.write("done")
        logging.info(f"{tool} alignment completed: {output_file}")
        return output_file
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"{tool} failed with error: {e}")

def predict_lncRNA(config, args, diamond_out_files, blastn_out):
    """Predict long non-coding RNAs."""
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
            if (record.id not in diamond_hits and record.id in blastn_hits and
                len(record.seq) >= config["lncRNA"]["min_length"]):
                SeqIO.write(record, out, "fasta")

    with open(checkpoint_file, "w") as f:
        f.write("done")
    logging.info(f"lncRNA prediction completed: {lnc_file}")
    return lnc_file

def dna2pep(config, args):
    """Run TransDecoder to predict ORFs."""
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
    """Prepare ID mapping file in Parquet format by processing it in chunks."""
    db_dir = args.db_dir
    idmapping_gz = os.path.join(db_dir, config["database"]["idmapping"])
    idmapping_parquet = os.path.join(db_dir, "idmapping.parquet")

    if os.path.exists(idmapping_parquet):
        logging.info(f"ID mapping already prepared: {idmapping_parquet}")
        return idmapping_parquet

    # Define chunk size (adjust based on your memory constraints)
    chunk_size = 1_000_000
    schema = None  # To store the schema for consistency across chunks
    writer = None  # Parquet writer object

    try:
        # Read the gzipped file in chunks
        for i, chunk in enumerate(pd.read_csv(
            idmapping_gz,
            sep="\t",
            header=None,
            names=["uniprot_id", "type", "value"],
            dtype=str,  # Force all columns to be strings to avoid mixed types
            compression="gzip",
            chunksize=chunk_size
        )):
            # Convert chunk to pyarrow Table
            table = pa.Table.from_pandas(chunk)

            if i == 0:
                # Initialize the Parquet writer with the schema from the first chunk
                schema = table.schema
                writer = pq.ParquetWriter(idmapping_parquet, schema, compression="snappy")
                writer.write_table(table)
            else:
                # Ensure subsequent chunks match the schema
                table = table.cast(schema)
                writer.write_table(table)

            logging.info(f"Processed chunk {i + 1} with {len(chunk)} rows")

    finally:
        # Close the writer to ensure the file is properly finalized
        if writer:
            writer.close()

    logging.info(f"ID mapping prepared: {idmapping_parquet}")
    return idmapping_parquet

def parse_uniprot_dat(dat_file):
    """Parse UniProt .dat file for functional descriptions."""
    descriptions = {}
    with gzip.open(dat_file, "rt") as f:
        uniprot_id, desc = None, []
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
                uniprot_id, desc = None, []
    return descriptions

def parse_go_obo(obo_file):
    """Parse GO ontology file."""
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
    logging.info(f"Parsed {len(go_terms)} GO terms from {obo_file}")
    return go_terms

def parse_enzyme_dat(enzyme_file):
    """Parse enzyme.dat file."""
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
    logging.info(f"Parsed {len(enzyme_data)} enzyme entries from {enzyme_file}")
    return enzyme_data

def parse_cdd_metadata(cdd_tar_file, extract_dir):
    """Parse CDD metadata from tarball."""
    cdd_data = {}
    if not os.path.exists(extract_dir):
        with tarfile.open(cdd_tar_file, "r:gz") as tar:
            tar.extractall(extract_dir)
        logging.info(f"Extracted {cdd_tar_file} to {extract_dir}")

    for smp_file in glob.glob(os.path.join(extract_dir, "*.smp")):
        with open(smp_file) as f:
            for line in f:
                if line.startswith(">"):
                    parts = line[1:].strip().split()
                    cdd_data[parts[0]] = " ".join(parts[1:])
    logging.info(f"Parsed {len(cdd_data)} CDD entries from {extract_dir}")
    return cdd_data

def parse_pfam2go(pfam2go_file):
    """Parse Pfam-to-GO mappings."""
    pfam_to_go = {}
    with open(pfam2go_file) as f:
        for line in f:
            if line.startswith("Pfam:"):
                parts = line.strip().split(" > ")
                pfam_id = parts[0].split()[0].replace("Pfam:", "")
                go_terms = parts[1].split(" ; ")
                pfam_to_go[pfam_id] = [go.split()[-1] for go in go_terms]
    logging.info(f"Parsed {len(pfam_to_go)} Pfam-to-GO mappings from {pfam2go_file}")
    return pfam_to_go

def build_output(config, args, diamond_out_files, blastn_out, rpstblastn_out, idmapping_parquet):
    """Build final annotations using MySQL."""
    output_dir = config["output"]["dir"]
    annotations_file = os.path.join(output_dir, config["output"]["annotations"])
    gff_dir = config["output"]["gff_dir"]
    checkpoint_file = os.path.join(output_dir, "output.done")

    if os.path.exists(checkpoint_file):
        logging.info(f"Output already built: {annotations_file}")
        return annotations_file

    os.makedirs(gff_dir, exist_ok=True)
    conn = mysql.connector.connect(**config["mysql"])
    cursor = conn.cursor()
    logging.info("Connected to MySQL database")

    cursor.execute(f"CREATE DATABASE IF NOT EXISTS {config['mysql']['database']}")
    cursor.execute(f"USE {config['mysql']['database']}")

    # Define table creation queries
    tables = {
        "diamond": """
            CREATE TABLE IF NOT EXISTS diamond (
                qseqid VARCHAR(255), sseqid VARCHAR(255), pident FLOAT, length INT, mismatch INT, gapopen INT,
                qstart INT, qend INT, sstart INT, send INT, evalue DOUBLE, bitscore FLOAT,
                PRIMARY KEY (qseqid, sseqid)
            )""",
        "blastn": """
            CREATE TABLE IF NOT EXISTS blastn (
                qseqid VARCHAR(255), sseqid VARCHAR(255), pident FLOAT, length INT, mismatch INT, gapopen INT,
                qstart INT, qend INT, sstart INT, send INT, evalue DOUBLE, bitscore FLOAT,
                PRIMARY KEY (qseqid, sseqid)
            )""",
        "rpstblastn": """
            CREATE TABLE IF NOT EXISTS rpstblastn (
                qseqid VARCHAR(255), sseqid VARCHAR(255), pident FLOAT, length INT, mismatch INT, gapopen INT,
                qstart INT, qend INT, sstart INT, send INT, evalue DOUBLE, bitscore FLOAT,
                PRIMARY KEY (qseqid, sseqid)
            )""",
        "idmapping": """
            CREATE TABLE IF NOT EXISTS idmapping (
                uniprot_id VARCHAR(255), type VARCHAR(50), value VARCHAR(255), INDEX (uniprot_id, type)
            )"""
    }
    for table, query in tables.items():
        cursor.execute(query)

    # Load data into tables
    for diamond_file in diamond_out_files:
        cursor.execute(f"LOAD DATA LOCAL INFILE '{diamond_file}' INTO TABLE diamond FIELDS TERMINATED BY '\t' "
                       "(qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)")
        logging.info(f"Loaded DIAMOND data from {diamond_file}")

    if blastn_out:
        cursor.execute(f"LOAD DATA LOCAL INFILE '{blastn_out}' INTO TABLE blastn FIELDS TERMINATED BY '\t' "
                       "(qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)")
        logging.info(f"Loaded BLASTN data from {blastn_out}")

    if rpstblastn_out:
        cursor.execute(f"LOAD DATA LOCAL INFILE '{rpstblastn_out}' INTO TABLE rpstblastn FIELDS TERMINATED BY '\t' "
                       "(qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)")
        logging.info(f"Loaded RPSBLASTN data from {rpstblastn_out}")

    idmapping_tsv = os.path.splitext(idmapping_parquet)[0] + ".tsv"
    if not os.path.exists(idmapping_tsv):
        pd.read_parquet(idmapping_parquet).to_csv(idmapping_tsv, sep="\t", index=False)
    cursor.execute(f"LOAD DATA LOCAL INFILE '{idmapping_tsv}' INTO TABLE idmapping FIELDS TERMINATED BY '\t' "
                   "(uniprot_id, type, value)")
    logging.info(f"Loaded IDmapping data from {idmapping_tsv}")

    # Create best-match view
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS diamond_best AS
        SELECT d.* FROM diamond d
        INNER JOIN (SELECT qseqid, MAX(bitscore) AS max_bitscore FROM diamond GROUP BY qseqid) db
        ON d.qseqid = db.qseqid AND d.bitscore = db.max_bitscore
    """)

    # Construct main query
    query = """
    SELECT d.qseqid, d.sseqid AS protein, d.evalue, d.bitscore, b.sseqid AS ncRNA, b.evalue AS nc_evalue,
           r.sseqid AS domain, r.evalue AS domain_evalue, GROUP_CONCAT(DISTINCT gm.value) AS go_ids,
           GROUP_CONCAT(DISTINCT em.value) AS ec_numbers, GROUP_CONCAT(DISTINCT pm.value) AS pfam_ids
    """
    if args.do_kegg_annotation or args.do_extended_annotation:
        query += ", GROUP_CONCAT(DISTINCT km.value) AS kegg_ids"
    if args.do_extended_annotation:
        query += ", GROUP_CONCAT(DISTINCT tm.value) AS tair_ids, GROUP_CONCAT(DISTINCT um.value) AS unipathway_ids"
    query += """
    FROM diamond_best d
    LEFT JOIN blastn b ON d.qseqid = b.qseqid
    LEFT JOIN rpstblastn r ON d.qseqid = r.qseqid
    LEFT JOIN idmapping gm ON d.sseqid = gm.uniprot_id AND gm.type = 'GO'
    LEFT JOIN idmapping em ON d.sseqid = em.uniprot_id AND em.type = 'EC'
    LEFT JOIN idmapping pm ON d.sseqid = pm.uniprot_id AND pm.type = 'Pfam'
    """
    if args.do_kegg_annotation or args.do_extended_annotation:
        query += "LEFT JOIN idmapping km ON d.sseqid = km.uniprot_id AND km.type = 'KEGG'"
    if args.do_extended_annotation:
        query += """
        LEFT JOIN idmapping tm ON d.sseqid = tm.uniprot_id AND tm.type = 'TAIR'
        LEFT JOIN idmapping um ON d.sseqid = um.uniprot_id AND um.type = 'UniPathway'
        """
    query += "GROUP BY d.qseqid, d.sseqid, d.evalue, d.bitscore, b.sseqid, b.evalue, r.sseqid, r.evalue"

    cursor.execute(query)
    columns = ["qseqid", "protein", "evalue", "bitscore", "ncRNA", "nc_evalue", "domain", "domain_evalue",
               "go_ids", "ec_numbers", "pfam_ids"]
    if args.do_kegg_annotation or args.do_extended_annotation:
        columns.append("kegg_ids")
    if args.do_extended_annotation:
        columns.extend(["tair_ids", "unipathway_ids"])
    result = pd.DataFrame(cursor.fetchall(), columns=columns)

    # Parse additional annotations
    go_terms = parse_go_obo(os.path.join(args.db_dir, config["database"]["go_obo"]))
    enzyme_data = parse_enzyme_dat(os.path.join(args.db_dir, config["database"]["enzyme_dat"]))
    cdd_data = parse_cdd_metadata(os.path.join(args.db_dir, config["database"]["cdd"]), os.path.join(args.db_dir, "cdd"))
    pfam_to_go = parse_pfam2go(os.path.join(args.db_dir, config["database"]["pfam2go"]))

    result["go_descriptions"] = result["go_ids"].apply(lambda x: map_go_terms(x, go_terms) if pd.notna(x) else None)
    result["enzyme_descriptions"] = result["ec_numbers"].apply(lambda x: map_enzyme(x, enzyme_data) if pd.notna(x) else None)
    result["domain_descriptions"] = result["domain"].apply(lambda x: map_cdd_domains(x, cdd_data) if pd.notna(x) else None)
    result["pfam_inferred_go"] = result["pfam_ids"].apply(lambda x: infer_go_from_pfam(x, pfam_to_go, go_terms) if pd.notna(x) else None)

    if args.do_kegg_annotation and "pathway_file" in config["database"]:
        pathways = {line.split("\t")[0].replace("PATH:", ""): line.split("\t")[1].strip()
                    for line in open(os.path.join(args.db_dir, config["database"]["pathway_file"])) if line.startswith("PATH:")}
        result["kegg_pathways"] = result["kegg_ids"].apply(lambda x: map_pathways(x, pathways) if pd.notna(x) else None)

    if args.do_function_description:
        descriptions = parse_uniprot_dat(os.path.join(args.db_dir, config["database"]["uniprot_dat"]))
        result["function_description"] = result["protein"].map(lambda x: descriptions.get(x.split("|")[1] if "|" in x else x, "No description"))

    result.to_csv(annotations_file, sep="\t", index=False)

    # Generate GFF3
    with open(os.path.join(gff_dir, "annotations.gff3"), "w") as gff:
        gff.write("##gff-version 3\n")
        cursor.execute("SELECT qseqid, sseqid, qstart, qend, bitscore, evalue FROM diamond_best WHERE sseqid IS NOT NULL")
        for row in cursor.fetchall():
            gff.write(f"{row[0]}\tDIAMOND\tmatch\t{row[2]}\t{row[3]}\t{row[4]}\t+\t.\tID={row[1]};evalue={row[5]}\n")

    with open(checkpoint_file, "w") as f:
        f.write("done")
    conn.commit()
    cursor.close()
    conn.close()
    logging.info(f"Output generated: {annotations_file}, GFF3 in {gff_dir}")
    return annotations_file

# Helper functions for mappings
def map_go_terms(go_ids, go_terms):
    return ";".join([f"{gid} ({go_terms.get(gid, 'Unknown')})" for gid in go_ids.split(",")])

def map_enzyme(ec_numbers, enzyme_data):
    return ";".join([f"{ec} ({enzyme_data.get(ec, 'Unknown')})" for ec in ec_numbers.split(",")])

def map_cdd_domains(domain, cdd_data):
    return f"{domain} ({cdd_data.get(domain.split('.')[0], 'Unknown')})"

def infer_go_from_pfam(pfam_ids, pfam_to_go, go_terms):
    inferred_go = set()
    for pfam in pfam_ids.split(","):
        inferred_go.update(pfam_to_go.get(pfam, []))
    return ";".join([f"{gid} ({go_terms.get(gid, 'Unknown')})" for gid in inferred_go]) if inferred_go else None

def map_pathways(kegg_ids, pathways):
    return ";".join([f"{kid} ({pathways.get(kid.split(':')[-1], 'Unknown')})" for kid in kegg_ids.split(",")])

def extract_statistics(config, args, annotations_file):
    """Extract and save statistics."""
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
    """Main execution logic."""
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
            with Pool(processes=min(len(tasks), args.threads)) as pool:
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
