import os
import subprocess
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

# --- CONFIGURATION ---
HMM_MODEL = "kunitz_model.hmm"
POS_FILE = "dataset_positives.fasta"
NEG_FILE = "dataset_negatives.fasta"
CPUS = 4
SEED = 42 # Ensures reproducibility of the random split

def count_sequences(fasta_file):
    """Counts sequences in a FASTA file using grep for speed."""
    try:
        result = subprocess.check_output(f"grep -c '^>' {fasta_file}", shell=True)
        return int(result.strip())
    except:
        return 0

def run_hmmsearch(input_fasta, output_tbl, total_db_size):
    """Executes hmmsearch with --max to disable heuristics."""
    print(f"Running hmmsearch on {input_fasta}...")
    cmd = [
        "hmmsearch", "--max", "-Z", str(total_db_size), "--cpu", str(CPUS),
        "--tblout", output_tbl, "-E", "1000",
        HMM_MODEL, input_fasta
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    print("Done.")

def load_data(fasta_file, tbl_file, label_value):
    """
    Parses FASTA to get all IDs and maps them to HMMER E-values.
    Assigns a dummy E-value (999) to sequences not detected by HMMER.
    """
    # 1. Load HMMER results
    hits = {}
    if os.path.exists(tbl_file):
        try:
            # Use sep='\s+' to handle variable whitespace
            df_hits = pd.read_csv(tbl_file, sep='\s+', comment='#', header=None, usecols=[0, 4], names=['id', 'evalue'])
            for _, row in df_hits.iterrows():
                hits[row['id']] = row['evalue']
        except Exception as e:
            print(f"Warning reading {tbl_file}: {e}")

    # 2. Map results to full sequence list
    data = []
    print(f"Indexing sequences from {fasta_file}...")
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id
        short_id = seq_id.split()[0] # Truncate ID at first space if necessary
        
        e_val = 999.0
        
        if seq_id in hits:
            e_val = hits[seq_id]
        elif short_id in hits:
            e_val = hits[short_id]
        
        data.append({
            'id': seq_id,
            'label': label_value, # 1 = Positive, 0 = Negative
            'evalue': e_val
        })
        
    return pd.DataFrame(data)

def calculate_mcc(tp, tn, fp, fn):
    """Calculates Matthews Correlation Coefficient."""
    numerator = (tp * tn) - (fp * fn)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denominator == 0: return 0
    return numerator / denominator

def test_threshold(df, threshold):
    """Computes metrics (MCC, ACC, Confusion Matrix) for a specific threshold."""
    pred = df['evalue'] <= threshold
    tp = ((df['label'] == 1) & (pred == True)).sum()
    tn = ((df['label'] == 0) & (pred == False)).sum()
    fp = ((df['label'] == 0) & (pred == True)).sum()
    fn = ((df['label'] == 1) & (pred == False)).sum()
    
    mcc = calculate_mcc(tp, tn, fp, fn)
    acc = (tp + tn) / len(df) if len(df) > 0 else 0
    return mcc, acc, (tp, tn, fp, fn)

def main():
    random.seed(SEED)
    
    # 1. Database Sizing (for E-value calculation)
    n_pos = count_sequences(POS_FILE)
    n_neg = count_sequences(NEG_FILE)
    total_z = n_pos + n_neg
    print(f"Total Database Size (Z): {total_z}")

    # 2. Run HMMSEARCH
    pos_tbl = "results_positives.tbl"
    neg_tbl = "results_negatives.tbl"
    
    if not os.path.exists(pos_tbl): run_hmmsearch(POS_FILE, pos_tbl, total_z)
    if not os.path.exists(neg_tbl): run_hmmsearch(NEG_FILE, neg_tbl, total_z)

    # 3. Load Data into Memory
    print("\n--- Loading and Merging Data ---")
    df_pos = load_data(POS_FILE, pos_tbl, 1)
    df_neg = load_data(NEG_FILE, neg_tbl, 0)
    
    # 4. Data Splitting (50/50 Shuffle)
    print("\n--- Creating Optimization Set (Subset 1) and Verification Set (Subset 2) ---")
    
    # Shuffle Positives
    pos_ids = df_pos.index.tolist()
    random.shuffle(pos_ids)
    split_p = int(len(pos_ids) * 0.5)
    pos_subset1 = df_pos.loc[pos_ids[:split_p]]
    pos_subset2 = df_pos.loc[pos_ids[split_p:]]
    
    # Shuffle Negatives
    neg_ids = df_neg.index.tolist()
    random.shuffle(neg_ids)
    split_n = int(len(neg_ids) * 0.5)
    neg_subset1 = df_neg.loc[neg_ids[:split_n]]
    neg_subset2 = df_neg.loc[neg_ids[split_n:]]
    
    # Combine
    df_set1 = pd.concat([pos_subset1, neg_subset1]) # Optimization
    df_set2 = pd.concat([pos_subset2, neg_subset2]) # Blind Test
    
    print(f"Subset 1: {len(df_set1)} sequences")
    print(f"Subset 2: {len(df_set2)} sequences")

    # 5. OPTIMIZATION PHASE (Subset 1)
    print("\n--- PHASE A: Optimization on Subset 1 ---")
    thresholds = [1e-10, 1e-5, 1e-4, 1e-3, 0.01, 0.1, 1, 10]
    results_set1 = []
    
    best_mcc = -1
    best_th = None
    
    print(f"{'Thresh':<10} | {'MCC (Set1)':<10} | {'ACC':<8}")
    for th in thresholds:
        mcc, acc, _ = test_threshold(df_set1, th)
        results_set1.append(mcc)
        print(f"{th:<10.1e} | {mcc:.4f}     | {acc:.4f}")
        
        if mcc > best_mcc:
            best_mcc = mcc
            best_th = th
            
    print(f"\n>>> Best Threshold found on Subset 1: {best_th} (MCC: {best_mcc:.4f})")

    # 6. VERIFICATION PHASE (Subset 2)
    print("\n--- PHASE B: Verification (Blind Test) on Subset 2 ---")
    mcc_ver, acc_ver, cm_ver = test_threshold(df_set2, best_th)
    tp, tn, fp, fn = cm_ver
    
    print(f"Applying threshold {best_th} to Subset 2:")
    print(f"Verification MCC: {mcc_ver:.4f}")
    print(f"Accuracy:         {acc_ver:.4f}")
    print(f"Confusion Matrix (Set 2): TP={tp}, FP={fp}, FN={fn}, TN={tn}")
    
    # Check Consistency
    diff = abs(best_mcc - mcc_ver)
    if diff < 0.02:
        print(">>> RESULT: Model is ROBUST (Similar performance across subsets).")
    else:
        print(f">>> RESULT: Warning, discrepancy found ({diff:.4f}). Possible overfitting.")

    # 7. FINAL RESULTS (Full Dataset)
    print("\n--- PHASE C: Final Results (Full Dataset) ---")
    df_full = pd.concat([df_pos, df_neg])
    mcc_fin, acc_fin, cm_fin = test_threshold(df_full, best_th)
    tp, tn, fp, fn = cm_fin
    
    print(f"Final Threshold: {best_th}")
    print(f"Final MCC: {mcc_fin:.4f}")
    print(f"Final Confusion Matrix:\nTP: {tp}\tFP: {fp}\nFN: {fn}\tTN: {tn}")

    # Plot generation
    plt.figure(figsize=(8, 5))
    plt.plot([str(t) for t in thresholds], results_set1, marker='o', label='Optimization (Set 1)')
    plt.axhline(y=mcc_ver, color='r', linestyle='--', label=f'Verification MCC (Set 2) at {best_th}')
    plt.title(f'Threshold Optimization (Best: {best_th})')
    plt.xlabel('E-value')
    plt.ylabel('MCC')
    plt.legend()
    plt.grid(True)
    plt.savefig("split_validation_plot.png")
    print("\nPlot saved: split_validation_plot.png")

if __name__ == "__main__":
    main()
