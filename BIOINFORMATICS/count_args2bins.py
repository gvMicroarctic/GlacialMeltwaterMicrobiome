#!/usr/bin/env python3

import sys
import re
import pandas as pd

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <input_table.tsv> <output_prefix>")
    sys.exit(1)

input_file = sys.argv[1]
prefix = sys.argv[2]

# Load table
df = pd.read_csv(input_file, sep="\t")

# Extract bin name (e.g., bin_65)
def extract_bin(name):
    m = re.search(r"(bin_\d+)", name)
    return m.group(1) if m else "unknown_bin"

df["bin"] = df["new_protein_name"].apply(extract_bin)

# ---------- AMR_category ----------
cat = (
    df.groupby(["bin", "AMR_category"])
    .size()
    .reset_index(name="count")
    .pivot(index="bin", columns="AMR_category", values="count")
    .fillna(0)
    .astype(int)
)

cat.to_csv(f"{prefix}_bin_AMR_category.tsv", sep="\t")

# ---------- AMR_sub_class ----------
sub = (
    df.groupby(["bin", "AMR_sub_class"])
    .size()
    .reset_index(name="count")
    .pivot(index="bin", columns="AMR_sub_class", values="count")
    .fillna(0)
    .astype(int)
)

sub.to_csv(f"{prefix}_bin_AMR_sub_class.tsv", sep="\t")

# ---------- Resistance_mechanism ----------
mech = (
    df.groupby(["bin", "Resistance_mechanism"])
    .size()
    .reset_index(name="count")
    .pivot(index="bin", columns="Resistance_mechanism", values="count")
    .fillna(0)
    .astype(int)
)

mech.to_csv(f"{prefix}_bin_Resistance_mechanism.tsv", sep="\t")

# ---------- TOTAL COUNTS PER BIN ----------
totals = (
    df.groupby("bin")
    .size()
    .reset_index(name="total_AMR_genes")
)

totals.to_csv(f"{prefix}_bin_total_counts.tsv", sep="\t", index=False)

print("Done.")
