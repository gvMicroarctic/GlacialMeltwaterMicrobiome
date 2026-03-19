#!/usr/bin/env python3

import pandas as pd
import sys
import re

# ------------------------------------------------------------
# Command line arguments
# ------------------------------------------------------------
if len(sys.argv) != 4:
    print("\nUsage:")
    print("  python count_transposase_intervals.py transposase_clusters_flanks.txt pathogenic_intervals.tsv output_counts.tsv\n")
    sys.exit(1)

TRANS_FILE = sys.argv[1]
INTERVALS_FILE = sys.argv[2]
OUTPUT_FILE = sys.argv[3]

# ------------------------------------------------------------
# Helper: normalize contig names
# ------------------------------------------------------------
def normalize_contig(name):
    return re.sub(r'^bin_\d+_(.*)', r'\1', str(name))

# ------------------------------------------------------------
# 1. Load transposase clusters
# ------------------------------------------------------------
trans_df = pd.read_csv(TRANS_FILE, sep="\t")
trans_df["contig"] = trans_df["contig"].apply(normalize_contig)

# ------------------------------------------------------------
# 2. Load intervals (to filter clusters)
# ------------------------------------------------------------
intervals_df = pd.read_csv(INTERVALS_FILE, sep="\t")
intervals_df["contig"] = intervals_df["contig"].apply(normalize_contig)

# Keep only clusters present in intervals (match by contig and cluster_start/cluster_end)
merged_df = pd.merge(
    trans_df,
    intervals_df[["contig","cluster_start","cluster_end"]],
    on=["contig","cluster_start","cluster_end"],
    how="inner"
)

# ------------------------------------------------------------
# 3. Identify interval columns (skip first 6 columns)
# ------------------------------------------------------------
interval_cols = merged_df.columns[6:]

# ------------------------------------------------------------
# 4. Count proteins per interval
# ------------------------------------------------------------
counts = {}
for col in interval_cols:
    counts[col] = merged_df[col].apply(lambda x: 0 if x == "-" else len(str(x).split(","))).sum()

# ------------------------------------------------------------
# 5. Write output file
# ------------------------------------------------------------
counts_df = pd.DataFrame(list(counts.items()), columns=["interval", "protein_count"])
counts_df.to_csv(OUTPUT_FILE, sep="\t", index=False)

print(f"\n✔ Total proteins per interval written to {OUTPUT_FILE}\n")
