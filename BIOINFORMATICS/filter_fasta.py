#!/usr/bin/env python3

# Filter FASTA file into two outputs:
#   1) sequences whose IDs are listed in the CSV file
#   2) sequences NOT listed in the CSV file

import argparse
from Bio import SeqIO

# --- Command-line arguments ---
parser = argparse.ArgumentParser(
    description="Split a FASTA file into two based on contig names from a CSV (first column)."
)
parser.add_argument("csv_file", help="CSV file with contig names in the first column")
parser.add_argument("fasta_file", help="Input FASTA file")
parser.add_argument("included_fasta", help="Output FASTA for contigs listed in CSV")
parser.add_argument("excluded_fasta", help="Output FASTA for contigs NOT listed in CSV")
args = parser.parse_args()

# --- Step 1: Read contig names from CSV ---
contigs_to_keep = set()
with open(args.csv_file, "r") as f:
    for line in f:
        contig_name = line.strip().split(",")[0]
        if contig_name:
            contigs_to_keep.add(contig_name)

# --- Step 2: Filter FASTA file ---
included_count = 0
excluded_count = 0

with open(args.included_fasta, "w") as keep_f, open(args.excluded_fasta, "w") as drop_f:
    for record in SeqIO.parse(args.fasta_file, "fasta"):
        if record.id in contigs_to_keep:
            SeqIO.write(record, keep_f, "fasta")
            included_count += 1
        else:
            SeqIO.write(record, drop_f, "fasta")
            excluded_count += 1

print(f"Eukaryotic {included_count} sequences → {args.included_fasta}")
print(f"Prokaryotic {excluded_count} sequences → {args.excluded_fasta}")
