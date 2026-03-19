#!/usr/bin/env python3

import argparse
from collections import Counter

def parse_emapper_file(input_file):
    cog_counts = Counter()
    category_counts = Counter()

    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#'):  # skip header/comments
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue

            # Column 5 = eggNOG_OGs, take last entry after last comma, before last '@'
            cog_field = fields[4]
            if cog_field != '-':
                last_cog = cog_field.split(',')[-1]      # last comma-separated
                cog = last_cog.split('@')[0]             # take before '@'
                cog_counts[cog] += 1

            # Column 7 = COG_category
            category = fields[6]
            if category != '-':
                category_counts[category] += 1

    return cog_counts, category_counts


def write_counts(counts, output_file):
    with open(output_file, 'w') as f:
        f.write("ID\tCount\n")
        for k, v in sorted(counts.items()):
            f.write(f"{k}\t{v}\n")


def main():
    parser = argparse.ArgumentParser(description="Count COGs and COG categories from eggNOG-mapper output")
    parser.add_argument("-i", "--input", required=True, help="eggNOG-mapper annotation file")
    parser.add_argument("-p", "--prefix", required=True, help="Output prefix")
    args = parser.parse_args()

    cog_counts, category_counts = parse_emapper_file(args.input)

    write_counts(cog_counts, f"{args.prefix}_COGs.txt")
    write_counts(category_counts, f"{args.prefix}_COG_categories.txt")


if __name__ == "__main__":
    main()
