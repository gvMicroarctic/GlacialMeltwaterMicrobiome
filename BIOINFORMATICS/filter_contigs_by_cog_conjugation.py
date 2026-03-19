#!/usr/bin/env python3

import argparse
import re

def read_fasta_headers(fasta_file):
    headers = set()
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                headers.add(line[1:].strip())
    return headers


def read_cog_list(cog_file):
    cogs = set()
    with open(cog_file) as f:
        next(f)  # skip header
        for line in f:
            if line.strip():
                cogs.add(line.split("\t")[0])
    return cogs


def extract_contig_id(query):
    """
    final.contigs_1000_prokaryota_k143_4391351_1
    -> k143_4391351
    """
    m = re.search(r"(k\d+_\d+)", query)
    return m.group(1) if m else None


def extract_cogs(cog_field):
    """
    COG0583@1|root,COG0583@2|Bacteria,2IXKS@203682|Planctomycetes
    -> {COG0583, 2IXKS}
    """
    cogs = set()
    for item in cog_field.split(","):
        cogs.add(item.split("@")[0])
    return cogs


def main(args):
    fasta_headers = read_fasta_headers(args.fasta)
    target_cogs = read_cog_list(args.cog_list)

    selected_contigs = set()

    with open(args.emapper) as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue

            query = cols[0]
            cog_field = cols[4]

            contig_id = extract_contig_id(query)
            if contig_id is None:
                continue

            if contig_id not in fasta_headers:
                continue

            found_cogs = extract_cogs(cog_field)

            if found_cogs & target_cogs:
                selected_contigs.add(contig_id)

    with open(args.output, "w") as out:
        for contig in sorted(selected_contigs):
            out.write(f">{contig}\n")

    print(f"Number of selected contigs: {len(selected_contigs)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter contigs by COG matches from eggNOG annotations"
    )
    parser.add_argument("--fasta", required=True, help="FASTA file with contig sequences")
    parser.add_argument("--cog-list", required=True, help="COG list file")
    parser.add_argument("--emapper", required=True, help="eggNOG emapper annotations file")
    parser.add_argument("--output", required=True, help="Output file with selected contig headers")

    args = parser.parse_args()
    main(args)
