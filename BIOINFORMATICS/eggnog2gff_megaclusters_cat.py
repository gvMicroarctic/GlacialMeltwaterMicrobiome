#!/usr/bin/env python3
import csv
import argparse

# -------------------------------
# Arguments
# -------------------------------
parser = argparse.ArgumentParser(
    description="Update GFF using EggNOG annotations with transposase override."
)
parser.add_argument("gff_file", help="Input GFF file")
parser.add_argument("eggnog_file", help="EggNOG annotation file")
parser.add_argument("transposase_file", help="Transposase OG annotation file")
parser.add_argument("output_file", help="Output GFF file")
args = parser.parse_args()

# -------------------------------
# Read transposase COG IDs
# -------------------------------
transposase_cogs = set()

with open(args.transposase_file) as f:
    for line in f:
        if not line.strip():
            continue
        parts = line.rstrip().split("\t")
        if len(parts) >= 2:
            transposase_cogs.add(parts[0])

# -------------------------------
# Read EggNOG annotations
# -------------------------------
eggnog_col5 = {}
eggnog_col7 = {}
eggnog_col8 = {}

with open(args.eggnog_file) as f:
    lines = [line for line in f if not line.startswith("##")]
    reader = csv.DictReader(lines, delimiter="\t")

    # Strip leading # from header if present
    if reader.fieldnames[0].startswith("#"):
        reader.fieldnames = [x.lstrip("#") for x in reader.fieldnames]

    for row in reader:
        protein_id = row["query"]

        # EggNOG column 5 (OGs): take last and remove tax scope
        col5 = row.get("eggNOG_OGs", "")
        if col5:
            last_block = col5.split(",")[-1]
            col5 = last_block.split("@")[0]
        else:
            col5 = "Unclassified"

        eggnog_col5[protein_id] = col5
        eggnog_col7[protein_id] = row.get("COG_category", "Unclassified")
        eggnog_col8[protein_id] = row.get("Description", "Unclassified")

# -------------------------------
# Process GFF
# -------------------------------
output_lines = []
in_fasta = False

with open(args.gff_file) as f:
    for line in f:
        line = line.rstrip()

        if line.startswith("##FASTA"):
            in_fasta = True
            output_lines.append(line)
            continue

        if in_fasta:
            output_lines.append(line)
            continue

        if line.startswith("#"):
            output_lines.append(line)
            continue

        parts = line.split("\t")
        if len(parts) < 9:
            continue

        # Parse attributes
        attributes = parts[8].split(";")
        attr_dict = {a.split("=")[0]: a.split("=")[1] for a in attributes if "=" in a}

        # Store original ID
        original_id = attr_dict.get("ID", "NA")
        attr_dict["ID_old"] = original_id

        protein = original_id

        cog_id = eggnog_col5.get(protein, "Unclassified")
        description = eggnog_col8.get(protein, "Unclassified")

        # -------------------------------
        # Annotation logic
        # -------------------------------
        attr_dict["ID"] = cog_id
        attr_dict["cat1"] = description

        if original_id in transposase_cogs:
            attr_dict["cat"] = "transposase"
        else:
            attr_dict["cat"] = eggnog_col7.get(protein, "Unclassified")

        # Rebuild attribute column
        parts[8] = ";".join(f"{k}={v}" for k, v in attr_dict.items())
        output_lines.append("\t".join(parts))

# -------------------------------
# Write output
# -------------------------------
with open(args.output_file, "w") as out:
    for line in output_lines:
        out.write(line + "\n")

print(f"Updated GFF file written: {args.output_file}")
