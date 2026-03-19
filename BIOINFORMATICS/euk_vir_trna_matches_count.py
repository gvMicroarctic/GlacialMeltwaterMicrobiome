#!/usr/bin/env python3

import argparse
from collections import defaultdict

# -------------------------------
# Argument parsing
# -------------------------------
parser = argparse.ArgumentParser(
    description="Count contig–vOTU pairs, exclude mapped vOTUs, and add taxonomy"
)

parser.add_argument("-i", "--input", required=True,
                    help="Input votu_euk.m8 file")

parser.add_argument("-e", "--exclude", required=True,
                    help="votus_vs_euk.m8 file (vOTUs in column 2 will be excluded)")

parser.add_argument("--euk-tax", required=True,
                    help="Eukaryotic contig taxonomy CSV")

parser.add_argument("--votu-tax", required=True,
                    help="vOTU taxonomy file (tab-delimited)")

parser.add_argument("-o", "--output", required=True,
                    help="Output TSV file")

args = parser.parse_args()

# -------------------------------
# Read excluded vOTUs
# -------------------------------
excluded_votus = set()

with open(args.exclude) as f:
    for line in f:
        if not line.strip():
            continue
        fields = line.rstrip().split("\t")
        excluded_votus.add(fields[1])

# -------------------------------
# Read eukaryotic contig taxonomy
# contig, kingdom, phylum, class, ...
# -------------------------------
contig_tax = {}

with open(args.euk_tax) as f:
    for line in f:
        if not line.strip():
            continue
        fields = line.rstrip().split(",")
        contig_id = fields[0]
        phylum = fields[2]
        tax_class = fields[3]
        contig_tax[contig_id] = (phylum, tax_class)

# -------------------------------
# Read vOTU taxonomy
# contig, domain, phylum, class, ...
# -------------------------------
votu_tax = {}

with open(args.votu_tax) as f:
    for line in f:
        if not line.strip():
            continue
        fields = line.rstrip().split("\t")
        votu_id = fields[0]
        phylum = fields[2].replace("p__", "")
        tax_class = fields[3].replace("c__", "")
        votu_tax[votu_id] = (phylum, tax_class)

# -------------------------------
# Count contig <-> vOTU pairs
# -------------------------------
pair_counts = defaultdict(int)

with open(args.input) as f:
    for line in f:
        if not line.strip():
            continue

        fields = line.rstrip().split("\t")

        # UniRef field: tRNA-Glu(ctc)_c[239847,239931]_k143_374875
        contig = "_".join(fields[0].split("_")[2:])

        # tRNA-Glu(ctc)_c[239847,239931]_contig_2391
        votu = "_".join(fields[1].split("_")[2:])

        if votu in excluded_votus:
            continue

        pair_counts[(contig, votu)] += 1

# -------------------------------
# Write output
# -------------------------------
with open(args.output, "w") as out:
    out.write(
        "contig\tvotu\tcount\t"
        "contig_phylum\tcontig_class\t"
        "votu_phylum\tvotu_class\n"
    )

    for (contig, votu), count in sorted(pair_counts.items()):
        contig_phylum, contig_class = contig_tax.get(
            contig, ("Unclassified", "Unclassified")
        )
        votu_phylum, votu_class = votu_tax.get(
            votu, ("Unclassified", "Unclassified")
        )

        out.write(
            f"{contig}\t{votu}\t{count}\t"
            f"{contig_phylum}\t{contig_class}\t"
            f"{votu_phylum}\t{votu_class}\n"
        )
