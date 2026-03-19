#!/usr/bin/env python3

import pandas as pd
import sys

# ------------------------------------------------------------
# Command line arguments
# ------------------------------------------------------------
if len(sys.argv) != 3:
    print("\nUsage:")
    print("  python count_proteins_intervals.py pathogenic_intervals.tsv interval_counts.tsv\n")
    sys.exit(1)

INPUT_FILE = sys.argv[1]
OUTPUT_FILE = sys.argv[2]

# ------------------------------------------------------------
# 1. Read pathogenic intervals file
# ------------------------------------------------------------
df = pd.read_csv(INPUT_FILE, sep="\t")

# Identify interval columns (skip first 4: cluster_id, contig, cluster_start, cluster_end)
interval_cols = df.columns[4:]

# ------------------------------------------------------------
# 2. Count proteins per interval
# ------------------------------------------------------------
counts = {}
for col in interval_cols:
    # split by ',' if multiple proteins, ignore '-' entries
    counts[col] = df[col].apply(lambda x: 0 if x == "-" else len(str(x).split(","))).sum()

# ------------------------------------------------------------
# 3. Write output
# ------------------------------------------------------------
counts_df = pd.DataFrame(list(counts.items()), columns=["interval", "protein_count"])
counts_df.to_csv(OUTPUT_FILE, sep="\t", index=False)

print(f"\n✔ Interval protein counts written to {OUTPUT_FILE}\n")
