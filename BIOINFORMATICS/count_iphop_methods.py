#!/usr/bin/env python3

import csv
import sys
import re
from collections import Counter

# -------------------------
# Check command-line input
# -------------------------
if len(sys.argv) != 2:
    print("Usage: python3 count_methods.py <input_csv>")
    sys.exit(1)

input_csv = sys.argv[1]

method_counts = Counter()

# -------------------------
# Read and parse the CSV
# -------------------------
with open(input_csv, newline="") as f:
    reader = csv.DictReader(f)

    for row in reader:
        methods_field = row["List of methods"].strip()
        if not methods_field:
            continue

        # Split by whitespace
        parts = re.split(r"\s+", methods_field)

        for part in parts:
            if ";" not in part:
                continue

            method, score_str = part.split(";")

            try:
                score = float(score_str)
            except ValueError:
                continue

            # Count only if score >= 90
            if score >= 90:
                method_counts[method] += 1

# -------------------------
# Output results
# -------------------------
print("Counts of methods with score >= 90:")
for method, count in method_counts.items():
    print(f"{method}: {count}")
