#!/usr/bin/env python3

import pandas as pd
import sys

# -----------------------------
# Command line input
# -----------------------------
if len(sys.argv) != 2:
    print("\nUsage:")
    print("  python count_unique_proteins.py flank_virulence_interval_counts.tsv\n")
    sys.exit(1)

INPUT = sys.argv[1]

# -----------------------------
# Load table
# -----------------------------
df = pd.read_csv(INPUT, sep="\t")

# -----------------------------
# Identify interval columns -10000 .. +10000
# -----------------------------
valid_intervals = []

for col in df.columns:

    if "_" not in col:
        continue

    # interval like "-9000_10000"
    body = col.replace("+", "").replace("-", "")
    try:
        a, b = map(int, body.split("_"))
    except:
        continue

    # keep columns where both ends lie between -10000 and +10000
    # original sign determines direction
    if col.startswith("-"):
        start = -b
        end = -a
    elif col.startswith("+"):
        start = a
        end = b
    else:
        continue
    
    if start >= -10000 and end <= 10000:
        valid_intervals.append(col)

# sort intervals human-friendly
def interval_sort_key(x):
    sgn = -1 if x.startswith("-") else 1
    vals = x.replace("+", "").replace("-", "")
    a, b = map(int, vals.split("_"))
    start = a * sgn
    return start

valid_intervals = sorted(valid_intervals, key=interval_sort_key)

# -----------------------------
# Count unique proteins
# -----------------------------
all_proteins = set()

for _, row in df.iterrows():
    for col in valid_intervals:
        val = row[col]
        if isinstance(val, str) and val != "-":
            for p in val.split(","):
                p = p.strip()
                if p:
                    all_proteins.add(p)

print(f"\n✔ Found {len(all_proteins)} unique proteins within -10kb..+10kb")
