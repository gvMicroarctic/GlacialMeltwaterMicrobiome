#!/usr/bin/env python3

import argparse

def load_high_quality_contigs(quality_file, completeness_threshold):
    """Read quality_summary.tsv and return contigs passing completeness cutoff"""
    good_contigs = set()

    with open(quality_file, "r") as f:
        header = f.readline().strip().split("\t")
        contig_idx = header.index("contig_id")
        completeness_idx = header.index("completeness")

        for line in f:
            fields = line.strip().split("\t")
            completeness = fields[completeness_idx]

            if completeness == "NA":
                continue

            if float(completeness) >= completeness_threshold:
                good_contigs.add(fields[contig_idx])

    return good_contigs


def main():
    parser = argparse.ArgumentParser(
        description="Parse Bacphlip output, filter by CheckV completeness, and classify viral lifestyle"
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
        "-q", "--quality",
        required=True,
        help="CheckV quality_summary.tsv file"
    )
    parser.add_argument(
        "-t", "--threshold",
        type=float,
        default=0.80,
        help="Probability threshold for lifestyle classification (default: 0.80)"
    )
    parser.add_argument(
        "-c", "--completeness_threshold",
        type=float,
        default=70.0,
        help="Minimum completeness percentage to keep contig (default: 70)"
    )

    args = parser.parse_args()

    # Load high-quality contigs
    good_contigs = load_high_quality_contigs(
        args.quality, args.completeness_threshold
    )

    count_virulent = 0
    count_temperate = 0
    count_none = 0
    count_filtered = 0

    with open(args.output, "w") as out_file, open(args.input, "r") as in_file:
        in_file.readline()  # skip header

        for line in in_file:
            info = line.strip().split("\t")
            contig_id = info[0]

            # Skip low-quality contigs
            if contig_id not in good_contigs:
                count_filtered += 1
                continue

            if float(info[1]) >= args.threshold:
                out_file.write(f"{contig_id}\tvirulent\n")
                count_virulent += 1
            elif float(info[2]) >= args.threshold:
                out_file.write(f"{contig_id}\ttemperate\n")
                count_temperate += 1
            else:
                out_file.write(f"{contig_id}\tnone\n")
                count_none += 1

    print(f"Virulent viruses: {count_virulent}")
    print(f"Temperate viruses: {count_temperate}")
    print(f"Unclassified viruses: {count_none}")
    print(f"Filtered out due to low completeness: {count_filtered}")


if __name__ == "__main__":
    main()
