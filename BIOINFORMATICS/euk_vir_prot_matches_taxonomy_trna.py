#!/usr/bin/env python3

import argparse
import csv
from collections import Counter


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize taxonomy of m8 hits by phylum/class/family"
    )
    parser.add_argument(
        "-m", "--m8",
        required=True,
        help="Input m8 file"
    )
    parser.add_argument(
        "-t", "--taxonomy",
        required=True,
        help="Taxonomy CSV file (contig in first column)"
    )
    parser.add_argument(
        "-x", "--exclude_m8",
        help="m8 file containing contigs to exclude (first column)"
    )
    parser.add_argument(
        "-o", "--out_prefix",
        required=True,
        help="Output prefix"
    )
    return parser.parse_args()


def load_taxonomy(tax_file):
    taxonomy = {}
    with open(tax_file) as f:
        reader = csv.reader(f)
        for row in reader:
            contig = row[0]
            taxonomy[contig] = {
                "phylum": row[2],
                "class": row[3],
                "family": row[5],
            }
    return taxonomy


def extract_contig(query_id):
    # Example:
    # tRNA-Glu(ctc)_c[239847,239931]_k143_374875
    parts = query_id.split("]_")
    for p in parts:
        if p.startswith("k143_"):
            return p
    return None


def load_excluded_contigs(m8_file):
    excluded = set()
    with open(m8_file) as f:
        for line in f:
            if not line.strip():
                continue
            query = line.split("\t")[0]
            contig = extract_contig(query)
            if contig:
                excluded.add(contig)
    return excluded


def summarize_m8(m8_file, taxonomy, excluded_contigs=None):
    if excluded_contigs is None:
        excluded_contigs = set()

    counts = {
        "phylum": Counter(),
        "class": Counter(),
        "family": Counter(),
    }

    with open(m8_file) as f:
        for line in f:
            if not line.strip():
                continue

            query = line.split("\t")[0]
            contig = extract_contig(query)

            if contig is None:
                continue

            if contig in excluded_contigs:
                continue

            tax = taxonomy.get(contig)
            if tax is None:
                continue

            for level in counts:
                counts[level][tax[level]] += 1

    return counts


def write_counts(counter, outfile):
    with open(outfile, "w") as f:
        f.write("taxonomy,count\n")
        for taxon, count in counter.most_common():
            f.write(f"{taxon},{count}\n")


def main():
    args = parse_args()

    taxonomy = load_taxonomy(args.taxonomy)

    excluded_contigs = set()
    if args.exclude_m8:
        excluded_contigs = load_excluded_contigs(args.exclude_m8)

    counts = summarize_m8(
        args.m8,
        taxonomy,
        excluded_contigs
    )

    write_counts(
        counts["phylum"],
        f"{args.out_prefix}.phylum_counts.csv"
    )
    write_counts(
        counts["class"],
        f"{args.out_prefix}.class_counts.csv"
    )
    write_counts(
        counts["family"],
        f"{args.out_prefix}.family_counts.csv"
    )


if __name__ == "__main__":
    main()
