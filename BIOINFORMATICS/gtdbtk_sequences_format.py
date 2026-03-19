#!/usr/bin/env python3

# ----------------------------
# 1. Input and output files
# ----------------------------
input_file = "../pp_db/tree_sequences.txt"        # file from previous script
output_file = "./MAG/phylophlan/db/accession/formatted_sequences.txt"  # new formatted file

# ----------------------------
# 2. Read sequences and format
# ----------------------------
formatted_sequences = []

with open(input_file, "r") as f:
    for line in f:
        seq = line.strip()
        # Skip bin entries
        if seq.startswith("bin_"):
            continue
        # Remove GB_ or RS_ prefix
        if seq.startswith("GB_GCA_"):
            seq_formatted = seq.replace("GB_", "")
        elif seq.startswith("GB_GCF_"):
            seq_formatted = seq.replace("GB_", "")
        elif seq.startswith("RS_GCA_"):
            seq_formatted = seq.replace("RS_", "")
        elif seq.startswith("RS_GCF_"):
            seq_formatted = seq.replace("RS_", "")
        else:
            seq_formatted = seq  # leave as is if no known prefix
        formatted_sequences.append(seq_formatted)

# ----------------------------
# 3. Write formatted sequences
# ----------------------------
with open(output_file, "w") as f:
    for seq in formatted_sequences:
        f.write(seq + "\n")

print(f"Formatted {len(formatted_sequences)} sequences and saved to {output_file}")
