#!/usr/bin/env python3

import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Parse Bacphlip output and classify viral lifestyle with length and taxonomy filtering"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Bacphlip output file (.bacphlip)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file with lifestyle classification"
    )
    parser.add_argument(
        "-l", "--lengths",
        required=True,
        help="File with contig lengths (contig_id <tab> length)"
    )
    parser.add_argument(
        "-x", "--taxonomy",
        required=True,
        help="Taxonomy file (contig_id d__ p__ c__ o__ f__)"
    )
    parser.add_argument(
        "-m", "--min_length",
        type=int,
        required=True,
        help="Minimum vOTU length to consider"
    )
    parser.add_argument(
        "-t", "--threshold",
        type=float,
        default=0.80,
        help="Probability threshold for classification (default: 0.80)"
    )

    args = parser.parse_args()

    # -----------------------------
    # Read contig lengths
    # -----------------------------
    contig_lengths = {}
    with open(args.lengths) as lf:
        for line in lf:
            contig, length = line.strip().split()
            contig_lengths[contig] = int(length)

    # -----------------------------
    # Read taxonomy and keep only p__ != Unclassified
    # -----------------------------
    valid_taxa = set()
    with open(args.taxonomy) as tf:
        for line in tf:
            fields = line.strip().split("\t")
            contig_id = fields[0]

            # Find phylum field (starts with p__)
            phylum = next(f for f in fields if f.startswith("p__"))
            if phylum != "p__Unclassified":
                valid_taxa.add(contig_id)

    count_virulent = 0
    count_temperate = 0
    count_none = 0
    count_filtered = 0

    with open(args.output, "w") as out_file, open(args.input, "r") as in_file:
        in_file.readline()  # skip header

        for line in in_file:
            info = line.strip().split("\t")
            contig_id = info[0]

            # Apply filters
            if contig_id not in contig_lengths:
                continue
            if contig_lengths[contig_id] < args.min_length:
                count_filtered += 1
                continue
            if contig_id not in valid_taxa:
                continue

            # Classification
            if float(info[1]) >= args.threshold:
                out_file.write(f"{contig_id}\tvirulent\n")
                count_virulent += 1
            elif float(info[2]) >= args.threshold:
                out_file.write(f"{contig_id}\ttemperate\n")
                count_temperate += 1
            else:
                out_file.write(f"{contig_id}\tnone\n")
                count_none += 1

    print(f"Virulent viruses are {count_virulent}.")
    print(f"Temperate viruses are {count_temperate}.")
    print(f"Unclassified viruses are {count_none}.")
    print(f"Filtered out due to length < {args.min_length}: {count_filtered}.")


if __name__ == "__main__":
    main()
