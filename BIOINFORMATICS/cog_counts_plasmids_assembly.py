#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Count COGs for genes whose contigs are present in a FASTA file"
    )
    parser.add_argument("-e", "--emapper", required=True,
                        help="eggNOG emapper annotation file")
    parser.add_argument("-f", "--fasta", required=True,
                        help="FASTA file with contig headers")
    parser.add_argument("-o", "--output", required=True,
                        help="Base name for output TSV files")
    return parser.parse_args()


def extract_cog_id(eggnog_ogs):
    try:
        last_block = eggnog_ogs.split(",")[-1]
        return last_block.split("@")[0]
    except Exception:
        return None


def main():
    args = parse_args()

    # --------------------------------------------------
    # 1. Read contig IDs from FASTA
    # --------------------------------------------------
    plasmid_contigs = set()
    with open(args.fasta) as f:
        for line in f:
            if line.startswith(">"):
                plasmid_contigs.add(line[1:].strip())

    # --------------------------------------------------
    # 2. Initialize counters
    # --------------------------------------------------
    cog_counts = defaultdict(lambda: [0, "", ""])
    category_counts = defaultdict(int)
    description_counts = defaultdict(int)

    # --------------------------------------------------
    # 3. Parse emapper and count
    # --------------------------------------------------
    with open(args.emapper) as f:
        reader = csv.reader(f, delimiter="\t")

        for row in reader:
            if not row or row[0].startswith("#"):
                continue

            query = row[0]
            eggnog_ogs = row[4]
            cog_category = row[6]
            description = row[7]

            if eggnog_ogs == "-" or eggnog_ogs == "":
                continue

            cog_id = extract_cog_id(eggnog_ogs)
            if cog_id is None:
                continue

            # final.contigs_1000_prokaryota_k143_4391351_1 → k143_4391351
            info = query.split("_")
            contig = info[3] + "_" + info[4]

            if contig not in plasmid_contigs:
                continue

            # COG-level
            cog_counts[cog_id][0] += 1
            cog_counts[cog_id][1] = cog_category
            cog_counts[cog_id][2] = description

            # Category-level
            if cog_category and cog_category != "-":
                category_counts[cog_category] += 1

            # Description-level
            if description and description != "-":
                description_counts[description] += 1

    # --------------------------------------------------
    # 4. Write outputs
    # --------------------------------------------------

    # 1) COG-level (unchanged)
    with open(f"{args.output}_by_COG.tsv", "w") as out:
        out.write("COG\tCount\tCOG_Category\tDescription\n")
        for cog, (count, cat, desc) in sorted(cog_counts.items()):
            out.write(f"{cog}\t{count}\t{cat}\t{desc}\n")

    # 2) Category-level
    with open(f"{args.output}_by_Category.tsv", "w") as out:
        out.write("COG_Category\tCount\n")
        for cat, count in sorted(category_counts.items()):
            out.write(f"{cat}\t{count}\n")

    # 3) Description-level
    with open(f"{args.output}_by_Description.tsv", "w") as out:
        out.write("Description\tCount\n")
        for desc, count in sorted(description_counts.items()):
            out.write(f"{desc}\t{count}\n")


if __name__ == "__main__":
    main()
