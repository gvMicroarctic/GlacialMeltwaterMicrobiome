#!/usr/bin/env python3
import csv
import argparse

parser = argparse.ArgumentParser(description="Match COG IDs between two files.")
parser.add_argument("transposase_file", help="Input file with transposase IDs (first column)")
parser.add_argument("eggnog_file", help="EggNOG emapper annotation file")
parser.add_argument("output_file", help="Output file for matched rows")
args = parser.parse_args()

# ---------------------------------------------------------
# Read first file: store first-column IDs
# ---------------------------------------------------------
transposase_ids = set()

with open(args.transposase_file, newline='') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        if not row:
            continue
        # Your transposase ID is in column 1 (second column)
        transposase_ids.add(row[1])

print(f"Loaded {len(transposase_ids)} transposase codes.")

# ---------------------------------------------------------
# Read EggNOG file and match ANY ID in column 5
# ---------------------------------------------------------
matches = []

with open(args.eggnog_file, newline='') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:

        # Skip comments or malformed lines
        if not row or row[0].startswith("#") or len(row) < 5:
            continue

        col5 = row[4]

        # Split by commas → ["COG0373@1|root", "COG0373@2|Bacteria", ...]
        entries = col5.split(",")

        # Extract all codes before "@"
        # e.g. "COG0373@1|root" → "COG0373"
        codes = [e.split("@")[0] for e in entries if "@" in e]

        # If ANY code matches transposase list
        if any(code in transposase_ids for code in codes):
            matches.append(row)

# ---------------------------------------------------------
# Write output
# ---------------------------------------------------------
with open(args.output_file, "w", newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerows(matches)

print(f"Saved {len(matches)} matched rows to {args.output_file}")
