#!/usr/bin/env python3

import sys
import pandas as pd

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <input_table.tsv> <output_prefix>")
    sys.exit(1)

input_file = sys.argv[1]
prefix = sys.argv[2]

# Read table
df = pd.read_csv(input_file, sep="\t")

# Count AMR_category
cat_counts = (
    df['AMR_category']
    .value_counts()
    .reset_index()
)
cat_counts.columns = ['AMR_category', 'count']
cat_counts.to_csv(f"{prefix}_AMR_category_counts.tsv", sep="\t", index=False)

# Count AMR_sub_class
sub_counts = (
    df['AMR_sub_class']
    .value_counts()
    .reset_index()
)
sub_counts.columns = ['AMR_sub_class', 'count']
sub_counts.to_csv(f"{prefix}_AMR_sub_class_counts.tsv", sep="\t", index=False)

# Count Resistance_mechanism
mech_counts = (
    df['Resistance_mechanism']
    .value_counts()
    .reset_index()
)
mech_counts.columns = ['Resistance_mechanism', 'count']
mech_counts.to_csv(f"{prefix}_Resistance_mechanism_counts.tsv", sep="\t", index=False)

print("Done.")
