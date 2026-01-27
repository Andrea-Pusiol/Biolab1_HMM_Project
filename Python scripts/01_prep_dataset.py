import os
import subprocess
from Bio import SeqIO

# --- CONFIGURATION ---
INPUT_PDB_FILE = "pdb_results.fasta" 
OUTPUT_BASE = "kunitz_dataset"
CDHIT_THRESHOLD = 0.9  # 90% sequence identity
MIN_LEN = 50
MAX_LEN = 90

def run_cdhit_pipeline():
    """
    Filters sequences by length and runs CD-HIT to remove redundancy.
    """
    cleaned_seqs = []
    original_count = 0
    
    # 1. Length Filtering
    print(f"--- 1. Sequence filtering (Length: {MIN_LEN}-{MAX_LEN} aa) ---")
    
    if not os.path.exists(INPUT_PDB_FILE):
        print(f"Error: File '{INPUT_PDB_FILE}' not found.")
        return

    for record in SeqIO.parse(INPUT_PDB_FILE, "fasta"):
        original_count += 1
        # Standardize sequence to uppercase
        seq_str = str(record.seq).upper()
        
        if MIN_LEN <= len(seq_str) <= MAX_LEN:
            cleaned_seqs.append(record)
        else:
            print(f"Discarded: {record.id} (Length: {len(seq_str)})")
            
    print(f"Original sequences: {original_count}")
    print(f"Valid sequences: {len(cleaned_seqs)}")
    
    # Save temporary file for CD-HIT input
    temp_fasta = "temp_cleaned.fasta"
    SeqIO.write(cleaned_seqs, temp_fasta, "fasta")
    
    # 2. CD-HIT Execution
    print(f"\n--- 2. Running CD-HIT (Threshold: {CDHIT_THRESHOLD}) ---")
    clustered_fasta = f"{OUTPUT_BASE}_clustered.fasta"
    
    cmd = [
        "cd-hit",
        "-i", temp_fasta,
        "-o", clustered_fasta,
        "-c", str(CDHIT_THRESHOLD),
        "-n", "5",  # Recommended word size for 0.9 threshold
        "-d", "0",  # Use full header description
        "-M", "16000"
    ]
    
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
        print("ERROR: 'cd-hit' tool not found. Please ensure it is installed.")
        return
    except subprocess.CalledProcessError as e:
        print(f"ERROR: CD-HIT execution failed: {e}")
        return

    # 3. Output Analysis
    final_seqs = list(SeqIO.parse(clustered_fasta, "fasta"))
    print(f"\n--- Final Output ---")
    print(f"Representative sequences (Training Set): {len(final_seqs)}")
    print(f"Output file: {clustered_fasta}")
    
    print("\n--- Manual Inspection Required ---")
    print("Please check the following sequences (potential outliers):")
    for rec in final_seqs:
        # Visual alert for Cysteine composition
        C_counter = 0
        for residue in rec.seq:
            if residue == "C":
                C_counter += 1
                
        if C_counter < 6:
            print(f" > {rec.id} | C composition low")

    # Cleanup
    if os.path.exists(temp_fasta):
        os.remove(temp_fasta)

if __name__ == "__main__":
    run_cdhit_pipeline()
