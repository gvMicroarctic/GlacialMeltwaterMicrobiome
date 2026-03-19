#!/usr/bin/env python3

import sys
from Bio import SeqIO
from collections import defaultdict

if len(sys.argv) != 3:
    print("Usage: python fasta_length_windows.py input.fasta output_prefix")
    sys.exit(1)

input_fasta = sys.argv[1]
prefix = sys.argv[2]

window_size = 1000

# Dictionaries to store counts
seq_counts = defaultdict(int)
bp_counts = defaultdict(int)

# Process each sequence in the FASTA
for record in SeqIO.parse(input_fasta, "fasta"):
    length = len(record.seq)
    window = (length // window_size) * window_size
    seq_counts[window] += 1
    bp_counts[window] += length

# Output 1: number of sequences per window
with open(f"{prefix}_sequence_counts.txt", "w") as f:
    f.write("Window_start\tNumber_of_sequences\n")
    for window in sorted(seq_counts):
        f.write(f"{window}\t{seq_counts[window]}\n")

# Output 2: number of base pairs per window
with open(f"{prefix}_bp_counts.txt", "w") as f:
    f.write("Window_start\tNumber_of_bp\n")
    for window in sorted(bp_counts):
        f.write(f"{window}\t{bp_counts[window]}\n")
