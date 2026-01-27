import sys

# --- CONFIGURATION ---
INPUT_FILE = "promals_output.txt"   # Paste the web content into this file
OUTPUT_FILE = "final_alignment.fasta" # Output for HMMER

def parse_promals_to_fasta():
    """
    Parses PROMALS3D text output into a clean FASTA file.
    Merges sequence blocks and ignores secondary structure (ss) lines.
    """
    sequences = {}
    order = [] # To maintain the original order

    print(f"--- Parsing {INPUT_FILE} ---")
    with open(INPUT_FILE,"r") as lines:
        for line in lines:
            line = line.strip()
            if not line: continue
            
            # Skip metadata headers and secondary structure lines
            if (line.startswith("blocksize:") or 
                line.startswith("alilen:") or 
                line.startswith("nal:") or 
                line.startswith("nblocks:") or 
                line.startswith("ss:") or 
                line.startswith("Final") or 
                line.startswith("Consensus_ss:")):
                continue
            
            # Parse sequence line (Format: Identifier Sequence)
            parts = line.split()
            if len(parts) < 2: continue 
            
            # The first part is the ID, the last part is the sequence segment
            seq_id = parts[0]
            seq_segment = parts[-1]
            
            # Initialize if new ID
            if seq_id not in sequences:
                sequences[seq_id] = ""
                order.append(seq_id)
            
            # Append sequence block
            sequences[seq_id] += seq_segment

    # Write to FASTA
    with open(OUTPUT_FILE, 'w') as out:
        count = 0
        for seq_id in order:
            # Write header and sequence (removing potential internal spaces)
            out.write(f">{seq_id}\n")
            out.write(f"{sequences[seq_id]}\n")
            count += 1
            
    print(f"--- Done. Wrote {count} sequences to {OUTPUT_FILE} ---")

if __name__ == "__main__":
    parse_promals_to_fasta()
