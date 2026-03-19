#!/usr/bin/env python3
import argparse
from collections import defaultdict
from Bio import SeqIO

# -----------------------------
# Arguments
# -----------------------------
parser = argparse.ArgumentParser(description="Extract contig sequences based on GFF CDS coordinates.")
parser.add_argument("gff", help="Input GFF file")
parser.add_argument("fasta", help="Input FASTA file")
parser.add_argument("output", help="Output FASTA file with extracted contigs")
args = parser.parse_args()

# -----------------------------
# Load GFF and get contig ranges
# -----------------------------
contig_ranges = defaultdict(lambda: [float('inf'), 0])  # contig -> [min_start, max_end]

with open(args.gff) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) != 9:
            continue
        contig, source, feature, start, end, score, strand, phase, attrs = parts
        if feature != "CDS":
            continue

        # Adjust contig name to match FASTA: remove bin_XX_ prefix if present
        if contig.startswith("bin_"):
            contig_fasta = "_".join(contig.split("_")[2:])  # NODE_51_bin_74
        else:
            contig_fasta = contig

        start = int(start)
        end = int(end)

        # Update min/max coordinates
        contig_ranges[contig_fasta][0] = min(contig_ranges[contig_fasta][0], start)
        contig_ranges[contig_fasta][1] = max(contig_ranges[contig_fasta][1], end)
        
# -----------------------------
# Read FASTA and extract sequences
# -----------------------------
output_records = []

for record in SeqIO.parse(args.fasta, "fasta"):
    contig_id = record.id
    if contig_id in contig_ranges:
        min_start, max_end = contig_ranges[contig_id]
        # GFF coordinates are 1-based inclusive
        seq_slice = record.seq[min_start-1:max_end]
        record.seq = seq_slice

        # Reformat header: >bin_52_k143_4527602_bin_52
        parts = contig_id.split("_")
        bin_number = parts[-1]
        prefix = "bin_" + bin_number
        rest = "_".join(parts[:-1])
        new_id = f"{prefix}_{rest}_{bin_number}"

        record.id = new_id
        record.description = ""  # remove everything else
        output_records.append(record)

# -----------------------------
# Write output FASTA
# -----------------------------
with open(args.output, "w") as out_f:
    SeqIO.write(output_records, out_f, "fasta")

print(f"Extracted {len(output_records)} contigs to {args.output}")
