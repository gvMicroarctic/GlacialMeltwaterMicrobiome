#!/usr/bin/env python3
import argparse
from Bio import SeqIO

# 1. Parse command-line arguments
parser = argparse.ArgumentParser(description="Filter FASTA sequences based on headers in a TSV file.")
parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
parser.add_argument("-a", "--annotations", required=True, help="Annotations TSV file (headers in first column)")
parser.add_argument("-o", "--output", required=True, help="Output filtered FASTA file")
args = parser.parse_args()

# 2. Read headers to remove from the annotations file
headers_to_remove = set()
with open(args.annotations, "r") as f:
    next(f)  # skip header line
    for line in f:
        parts = line.strip().split("\t")
        if parts:
            headers_to_remove.add(parts[0])

# 3. Filter FASTA sequences
with open(args.output, "w") as out_f:
    for record in SeqIO.parse(args.fasta, "fasta"):
        # Use record.id (first word of FASTA header)
        if record.id not in headers_to_remove:
            SeqIO.write(record, out_f, "fasta")

print(f"Filtered FASTA saved to {args.output}")
