import pandas as pd
import os

# --- CONFIGURATION ---
# Use the threshold obtained from script 04 (e.g., 1e-5 or 0.001)
THRESHOLD = 1e-05 
POS_TBL = "results_positives.tbl"
NEG_TBL = "results_negatives.tbl"

def load_hmm_results(tbl_file):
    """Parses HMMER table output into a DataFrame."""
    if not os.path.exists(tbl_file):
        print(f"Error: File '{tbl_file}' not found.")
        return None
    
    try:
        # Use sep='\s+' to handle variable whitespace in HMMER output
        df = pd.read_csv(tbl_file, sep='\s+', comment='#', header=None, 
                         usecols=[0, 4], names=['id', 'evalue'])
        return df
    except Exception as e:
        print(f"Error reading {tbl_file}: {e}")
        return None

def main():
    print(f"--- HMMER Result Investigation (Threshold: {THRESHOLD}) ---")

    # 1. Load Data
    df_pos = load_hmm_results(POS_TBL)
    df_neg = load_hmm_results(NEG_TBL)

    if df_pos is None or df_neg is None:
        return

    # 2. Analyze False Negatives (FN)
    # Definition: Positive sequences with E-value > Threshold (or not found)
    print("\n>>> ANALYSIS: False Negatives (FN)")
    print("Criteria: Known Kunitz proteins missed by the model (E-value > Threshold).")
    
    fns = df_pos[df_pos['evalue'] > THRESHOLD].sort_values(by='evalue')
    
    if not fns.empty:
        print("-" * 65)
        print(f"{'UniProt ID':<30} | {'E-value':<15} | {'Interpretation'}")
        print("-" * 65)
        
        for _, row in fns.iterrows():
            e_val = row['evalue']
            # Automatic interpretation based on score
            if e_val > 10:
                interp = "Non-functional / Outlier"
            elif e_val > 0.001:
                interp = "Weak Match / Divergent"
            else:
                interp = "Borderline (Strict Threshold)"
                
            print(f"{row['id']:<30} | {e_val:<15.2e} | {interp}")
        print("-" * 65)
    else:
        print("No False Negatives found (Perfect Recall).")

    # 3. Analyze False Positives (FP)
    # Definition: Negative sequences with E-value <= Threshold
    print("\n>>> ANALYSIS: False Positives (FP)")
    print("Criteria: Non-Kunitz proteins (background) detected as Kunitz.")
    
    fps = df_neg[df_neg['evalue'] <= THRESHOLD].sort_values(by='evalue')
    
    if not fps.empty:
        print("-" * 65)
        print(f"{'UniProt ID':<30} | {'E-value':<15} | {'Interpretation'}")
        print("-" * 65)
        
        for _, row in fps.iterrows():
            print(f"{row['id']:<30} | {row['evalue']:<15.2e} | Potential New Annotation?")
        print("-" * 65)
        print("Note: Check UniProt for these IDs. They might be unannotated Kunitz domains.")
    else:
        print("No False Positives found (Perfect Precision).")

if __name__ == "__main__":
    main()
