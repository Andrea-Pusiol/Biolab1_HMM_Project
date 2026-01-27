import os
import gzip
import requests
import subprocess
from Bio import SeqIO

# --- CONFIGURATION ---
TRAINING_SET_FILE = "kunitz_dataset_clustered.fasta" 
SWISSPROT_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
PFAM_ID = "PF00014" # Pfam ID for Kunitz domain
BLAST_IDENTITY_THRESHOLD = 95.0 # Identity threshold to avoid overfitting

def download_file(url, output_path):
    """Downloads a file if it does not exist locally."""
    if os.path.exists(output_path):
        print(f"File {output_path} already exists. Skipping download.")
        return
    print(f"Downloading {url}...")
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(output_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    print("Download complete.")

def get_kunitz_ids():
    """Retrieves list of UniProt IDs associated with the Kunitz domain (PF00014)."""
    print("Fetching Kunitz IDs from UniProt API...")
    # Query: database: (type:swiss) AND xref:pfam-PF00014
    url = "https://rest.uniprot.org/uniprotkb/stream?format=list&query=%28reviewed%3Atrue%29+AND+%28xref%3Apfam-PF00014%29"
    
    r = requests.get(url)
    ids = set(r.text.strip().split("\n"))
    print(f"Found {len(ids)} Kunitz entries in Swiss-Prot.")
    return ids

def run_blast_filter(positives_file, training_file, final_output):
    """
    Filters out positive sequences that are too similar to the training set 
    (Identity >= 95%) using BLASTP to ensure fair validation.
    """
    print("\n--- BLAST FILTERING (Removing Training Redundancy) ---")
    
    # 1. Create BLAST database from training set
    print("Creating local BLAST database...")
    cmd_db = f"makeblastdb -in {training_file} -dbtype prot -out temp_training_db"
    subprocess.run(cmd_db, shell=True, check=True)
    
    # 2. Run Blastp (Query: Positives, DB: Training)
    print("Running BLASTP...")
    blast_out = "blast_results.txt"
    # Output fmt 6: qseqid sseqid pident ...
    cmd_blast = f"blastp -query {positives_file} -db temp_training_db -out {blast_out} -outfmt 6 -evalue 1e-5"
    subprocess.run(cmd_blast, shell=True, check=True)
    
    # 3. Parse results
    to_exclude = set()
    with open(blast_out, "r") as f:
        for line in f:
            cols = line.split("\t")
            q_id = cols[0] # Query ID (Positive)
            pident = float(cols[2]) # % Identity
            
            if pident >= BLAST_IDENTITY_THRESHOLD:
                to_exclude.add(q_id)
    
    print(f"Sequences excluded due to similarity (>={BLAST_IDENTITY_THRESHOLD}%): {len(to_exclude)}")
    
    # 4. Write final file
    count = 0
    with open(final_output, "w") as out_handle:
        for record in SeqIO.parse(positives_file, "fasta"):
            # Clean ID extraction (handling formats like sp|ID|Name)
            clean_id = record.id.split("|")[1] if "|" in record.id else record.id
            
            if clean_id not in to_exclude and record.id not in to_exclude:
                SeqIO.write(record, out_handle, "fasta")
                count += 1
                
    print(f"Final Positive Dataset created: {final_output} ({count} sequences)")
    
    # Cleanup temporary BLAST files
    for ext in ["phr", "pin", "psq", "pdb", "pot", "ptf", "pto"]:
        try: os.remove(f"temp_training_db.{ext}") 
        except: pass
    if os.path.exists(blast_out): os.remove(blast_out)

def main():
    # 1. Download Swiss-Prot
    sp_file_gz = "uniprot_sprot.fasta.gz"
    download_file(SWISSPROT_URL, sp_file_gz)
    
    # 2. Get True Positive IDs
    kunitz_ids = get_kunitz_ids()
    
    # 3. Split Swiss-Prot into Raw Positives and Negatives
    print("\nSplitting Swiss-Prot into Positives (Raw) and Negatives...")
    pos_raw = "positives_raw.fasta"
    neg_file = "dataset_negatives.fasta"
    
    count_pos = 0
    count_neg = 0
    
    with gzip.open(sp_file_gz, "rt") as handle:
        with open(pos_raw, "w") as f_pos, open(neg_file, "w") as f_neg:
            for record in SeqIO.parse(handle, "fasta"):
                try:
                    # Extract Unique ID (e.g., Q9UII4 from sp|Q9UII4|...)
                    uid = record.id.split("|")[1]
                except IndexError:
                    uid = record.id
                
                if uid in kunitz_ids:
                    SeqIO.write(record, f_pos, "fasta")
                    count_pos += 1
                else:
                    SeqIO.write(record, f_neg, "fasta")
                    count_neg += 1
                    
    print(f"Splitting complete: {count_pos} Positives (raw), {count_neg} Negatives.")

    # 4. Filter Positives (Redundancy check)
    final_pos = "dataset_positives.fasta"
    
    if os.path.exists(TRAINING_SET_FILE):
        run_blast_filter(pos_raw, TRAINING_SET_FILE, final_pos)
    else:
        print(f"WARNING: Training file {TRAINING_SET_FILE} not found. Skipping BLAST filter.")
        os.rename(pos_raw, final_pos)

    # Cleanup raw positive file
    if os.path.exists(pos_raw): os.remove(pos_raw)
    
    print("\n--- PROCESS COMPLETED ---")
    print(f"1. Positive Dataset: {final_pos}")
    print(f"2. Negative Dataset: {neg_file}")

if __name__ == "__main__":
    main()
