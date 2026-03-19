#!/usr/bin/env python3

import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(
        description="Remove sequences from a FASTA file based on header names"
    )
    parser.add_argument(
        "--sequence", required=True,
        help="Comma-separated list of sequence headers to remove"
    )
    parser.add_argument(
        "--fasta", required=True,
        help="Input FASTA file"
    )
    parser.add_argument(
        "--output", required=True,
        help="Output FASTA file after removing sequences"
    )
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Create a set of sequence headers to remove
    remove_set = set(args.sequence.split(","))
    
    with open(args.fasta) as infile, open(args.output, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id not in remove_set:
                SeqIO.write(record, outfile, "fasta")
    
    print(f"Removed {len(remove_set)} sequences. Output saved to {args.output}")

if __name__ == "__main__":
    main()
