#!/usr/bin/env python3
import csv
import argparse

parser = argparse.ArgumentParser(description="Extract CDS in ±window around megaclusters.")
parser.add_argument("gff", help="Input GFF file")
parser.add_argument("megaclusters", help="Megacluster file")
parser.add_argument("output", help="Output GFF file")
parser.add_argument("-w", "--window", type=int, default=2000,
                    help="Window size to extend megacluster (default 2000 bp)")
args = parser.parse_args()

WINDOW = args.window

# ---------------------------------------------------------
# Load megaclusters
# ---------------------------------------------------------
megaregions = {}  # contig → list of (expanded_start, expanded_end, protein_ids)

with open(args.megaclusters) as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        contig = row["contig"]
        start = int(row["start"])
        end = int(row["end"])
        proteins = row["proteins"].split(",")

        exp_start = max(1, start - WINDOW)
        exp_end = end + WINDOW

        megaregions.setdefault(contig, []).append((exp_start, exp_end, set(proteins)))

# ---------------------------------------------------------
# Read GFF and filter
# ---------------------------------------------------------
output_lines = []
header_lines = []

with open(args.gff) as f:
    for line in f:
        if line.startswith("#"):
            header_lines.append(line)
            continue

        parts = line.strip().split("\t")
        if len(parts) != 9:
            continue

        contig, source, feature, start, end, score, strand, phase, attributes = parts
        if feature != "CDS":
            continue

        start = int(start)
        end = int(end)

        # Extract ID=xxx from attributes
        id_field = None
        for item in attributes.split(";"):
            if item.startswith("ID="):
                id_field = item.split("=", 1)[1]
                break
        if id_field is None:
            continue

        # Check if CDS matches any megacluster region
        if contig in megaregions:
            for exp_start, exp_end, proteins in megaregions[contig]:
                if id_field in proteins or not (end < exp_start or start > exp_end):
                    output_lines.append(line)
                    break

# ---------------------------------------------------------
# Write output GFF
# ---------------------------------------------------------
with open(args.output, "w") as out:
    # Write unique headers only once
    seen_headers = set()
    for h in header_lines:
        if h not in seen_headers:
            out.write(h)
            seen_headers.add(h)

    for l in output_lines:
        out.write(l)

print(f"Written {len(output_lines)} CDS entries to {args.output}")
