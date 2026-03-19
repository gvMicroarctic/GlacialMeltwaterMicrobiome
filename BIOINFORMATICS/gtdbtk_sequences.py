#!/usr/bin/env python3

from Bio import Phylo

# ----------------------------
# 1. Input and output files
# ----------------------------
tree_file = "./MAG/gtdbtk_drep/classify/gtdbtk.backbone.bac120.classify.tree"  # your Newick tree
output_file = "../pp_db/tree_sequences.txt"

# ----------------------------
# 2. Load tree
# ----------------------------
tree = Phylo.read(tree_file, "newick")

# ----------------------------
# 3. Extract all leaf names
# ----------------------------
sequences = []
for leaf in tree.get_terminals():
    if leaf.name:
        sequences.append(leaf.name.strip("'"))  # remove quotes if any

# ----------------------------
# 4. Write to output file
# ----------------------------
with open(output_file, "w") as f:
    for seq in sequences:
        f.write(seq + "\n")

print(f"Extracted {len(sequences)} sequences to {output_file}")
