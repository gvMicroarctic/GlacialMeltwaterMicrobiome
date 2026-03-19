#!/usr/bin/env python3
import csv
import argparse

# -------------------------------
# Arguments
# -------------------------------
parser = argparse.ArgumentParser(description="Update GFF IDs using EggNOG annotations (two outputs) and preserve sequences.")
parser.add_argument("gff_file", help="Input GFF file")
parser.add_argument("eggnog_file", help="EggNOG annotation file")
parser.add_argument("output_file1", help="Output GFF file using 5th column")
parser.add_argument("output_file2", help="Output GFF file using 8th column")
args = parser.parse_args()

# -------------------------------
# Read EggNOG annotations
# -------------------------------
eggnog_map_col5 = {}  # protein_id -> last value of column 5 before '@'
eggnog_map_col8 = {}  # protein_id -> last value of column 8 before '@'

with open(args.eggnog_file) as f:
    lines = [line for line in f if not line.startswith("##")]
    reader = csv.DictReader(lines, delimiter="\t")

    # Strip leading # from header if present
    if reader.fieldnames[0].startswith("#"):
        reader.fieldnames = [x.lstrip("#") for x in reader.fieldnames]

    for row in reader:
        protein_id = row["query"]
        
        # Column 5 mapping
        col5 = row.get("eggNOG_OGs", "")
        if col5:
            last_block = col5.split(",")[-1]
            last_value = last_block.split("@")[0]
            eggnog_map_col5[protein_id] = last_value

        # Column 8 mapping
        col8 = row.get("Description", "")
        if col8:
            last_block = col8.split(",")[-1]
            last_value = last_block.split("@")[0]
            eggnog_map_col8[protein_id] = last_value

# -------------------------------
# Read GFF and update IDs
# -------------------------------
gff_lines_col5 = []
gff_lines_col8 = []

in_fasta = False  # Flag to detect ##FASTA section

with open(args.gff_file) as f:
    for line in f:
        line = line.rstrip()

        if line.startswith("##FASTA"):
            in_fasta = True
            # Keep the ##FASTA line in both outputs
            gff_lines_col5.append(line)
            gff_lines_col8.append(line)
            continue

        if in_fasta:
            # Everything after ##FASTA is just sequences
            gff_lines_col5.append(line)
            gff_lines_col8.append(line)
            continue

        if line.startswith("#"):
            # Keep all other comments
            gff_lines_col5.append(line)
            gff_lines_col8.append(line)
            continue

        parts = line.split("\t")
        if len(parts) < 9:
            # Not a valid GFF line, skip
            continue

        # extract protein ID from original ID= field
        protein = None
        attributes = parts[8].split(";")
        for attr in attributes:
            if attr.startswith("ID="):
                protein = attr.split("=")[1]
                break

        # Lookup EggNOG mapping
        new_id5 = eggnog_map_col5.get(protein, "Unclassified")
        new_id8 = eggnog_map_col8.get(protein, "Unclassified")

        # Update last column but keep other attributes intact
        parts5 = parts.copy()
        attr_dict5 = {a.split("=")[0]: a.split("=")[1] for a in attributes if "=" in a}
        attr_dict5["ID"] = new_id5
        parts5[8] = ";".join(f"{k}={v}" for k, v in attr_dict5.items())
        gff_lines_col5.append("\t".join(parts5))

        parts8 = parts.copy()
        attr_dict8 = {a.split("=")[0]: a.split("=")[1] for a in attributes if "=" in a}
        attr_dict8["ID"] = new_id8
        parts8[8] = ";".join(f"{k}={v}" for k, v in attr_dict8.items())
        gff_lines_col8.append("\t".join(parts8))

# -------------------------------
# Write output GFFs
# -------------------------------
with open(args.output_file1, "w") as out1:
    for line in gff_lines_col5:
        out1.write(line + "\n")

with open(args.output_file2, "w") as out2:
    for line in gff_lines_col8:
        out2.write(line + "\n")

print(f"Updated GFF files written:\n  {args.output_file1}\n  {args.output_file2}")
