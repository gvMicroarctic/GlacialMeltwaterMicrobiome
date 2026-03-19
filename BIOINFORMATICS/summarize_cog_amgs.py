#!/usr/bin/env python3

import csv
import sys
from collections import Counter

# --- Command Line Argument Setup ---
if len(sys.argv) != 4:
    script_name = sys.argv[0].split('/')[-1]
    print(f"Usage: python {script_name} <gene_list_file> <cog_annotation_file> <output_tsv_file>", file=sys.stderr)
    print("\nExample:", file=sys.stderr)
    print(f"python {script_name} intervening_genes.txt emapper_cog.tsv cog_summary.tsv", file=sys.stderr)
    sys.exit(1)

gene_list_path = sys.argv[1]
cog_annotation_path = sys.argv[2]
output_prefix = sys.argv[3]

# Counters for COG terms and descriptions
cog_counts = Counter()
description_counts = Counter()
cog_descriptions = {}

# New: dictionary for gene-to-COG mapping
gene_to_cog = {}

target_genes = set()

# -------------------------
# Stage 1: Read Gene List
# -------------------------
print(f"Reading gene list from: {gene_list_path}")
with open(gene_list_path, 'r') as f:
    for line in f:
        gene = line.strip()
        if gene:
            target_genes.add(gene)
print(f"Found {len(target_genes)} target genes to annotate.")

# -------------------------
# Stage 2: Process COG Annotations
# -------------------------
print(f"Processing COG annotations from: {cog_annotation_path}")
with open(cog_annotation_path, 'r', newline='') as f:
    reader = csv.reader(f, delimiter='\t')
    
    for row in reader:
        if not row:
            continue
        if row[0].startswith('#'):
            continue

        query_gene = row[0]
        if query_gene in target_genes:
            cog0 = row[4]          # COG ID column
            cog_category = row[6]  # Single-letter COG category
            description = row[7]   # Description
            
            last_block = cog0.split(",")[-1]
            cog = last_block.split("@")[0].strip()
            
            cog_counts[cog] += 1
            description_counts[description] += 1
            cog_descriptions[cog] = (description, cog_category)

            # New: record gene-to-COG mapping
            gene_to_cog[query_gene] = cog

# -------------------------
# Stage 3: Write COG summary
# -------------------------
cog_out = output_prefix + ".cog.txt"
with open(cog_out, 'w', newline='') as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    writer.writerow(["COG_Term", "Count", "Description", "COG_category"])
    sorted_counts = sorted(cog_counts.items(), key=lambda item: item[1], reverse=True)
    for cog, count in sorted_counts:
        description = cog_descriptions.get(cog, ('-', '-'))
        writer.writerow([cog, count, description[0], description[1]])
print(f"Wrote summary for {len(cog_counts)} unique COG terms to {cog_out}.")

# -------------------------
# Stage 4: Write Description summary
# -------------------------
cog_desc_out = output_prefix + ".cog_description.txt"
with open(cog_desc_out, 'w', newline='') as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    writer.writerow(["Description", "Count"])
    sorted_desc = sorted(description_counts.items(), key=lambda item: item[1], reverse=True)
    for description, count in sorted_desc:
        writer.writerow([description, count])
print(f"Wrote summary for {len(description_counts)} unique descriptions to {cog_desc_out}.")

# -------------------------
# Stage 5: Write gene-to-COG mapping
# -------------------------
gene_cog_out = output_prefix + ".gene_to_cog.txt"
with open(gene_cog_out, 'w', newline='') as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    writer.writerow(["Gene", "COG_Term"])
    for gene, cog in gene_to_cog.items():
        writer.writerow([gene, cog])
print(f"Wrote gene-to-COG mapping for {len(gene_to_cog)} genes to {gene_cog_out}.")
