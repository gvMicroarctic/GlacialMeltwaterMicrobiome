#!/usr/bin/env python3

import sys
import csv
from collections import defaultdict

# --- STEP 0: Parse command-line arguments ---
if len(sys.argv) != 4:
    print("Usage: python cog_flank_matrix.py <cog_file> <flank_file> <output_prefix>")
    sys.exit(1)

cog_file = sys.argv[1]
flank_file = sys.argv[2]
output_prefix = sys.argv[3]

# ===============================================================
# STEP 1 — GET COG DATA
# ===============================================================
gene_to_cog_cat = {}
gene_to_cog = {}
cog_to_description = {}
cog_to_cog_cat = {}

with open(cog_file) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue

        parts = line.strip().split("\t")
        if len(parts) < 8:
            continue
        gene = parts[0]

        # ---- A) COG CATEGORY FROM COLUMN 7----
        cog_cat = parts[6].strip()
        if not cog_cat or cog_cat == "-":
            cog_cat = "No_ANNOTATION"
        gene_to_cog_cat[gene] = cog_cat

        # ---- B) LAST COG FROM COLUMN 5 ----
        col5 = parts[4]
        last_block = col5.split(",")[-1]
        cog = last_block.split("@")[0].strip()
        if cog == "":
            cog = "No_ANNOTATION"
        gene_to_cog[gene] = cog
        
        # ---- c) DESCRIPTION FROM COLUMN 8 ----
        cog_description = parts[7].strip()
        if not cog_description or cog_description == "-":
            cog_description = "No_DESCRIPTION"
        cog_to_description[cog] = cog_description
        
        cog_to_cog_cat[cog] = cog_cat

print(f"Loaded {len(gene_to_cog_cat)} COG categories.")
print(f"Loaded {len(gene_to_cog)} COGs.")
print(f"Loaded {len(cog_to_description)} COG descriptions.")

# ===============================================================
# STEP 2 — PARSE FLANK FILE AND COUNT
# ===============================================================
cog_cat_counts = defaultdict(lambda: defaultdict(int))
cog_counts = defaultdict(lambda: defaultdict(int))

all_cog_categories = set()
all_cogs = set()
all_flanks = []

with open(flank_file) as f:
    reader = csv.DictReader(f, delimiter="\t")
    flank_columns = [c for c in reader.fieldnames if c.startswith("-") or c.startswith("+")]
    all_flanks = flank_columns

    for row in reader:
        for flank in flank_columns:
            genes_str = row[flank].strip()

            if genes_str == "-":
                cog_cat_counts["No_PROTEIN"][flank] += 1
                cog_counts["No_PROTEIN"][flank] += 1
                all_cog_categories.add("No_PROTEIN")
                all_cogs.add("No_PROTEIN")
                continue

            genes = [g.strip() for g in genes_str.split(",") if g.strip()]

            for gene in genes:
                # ---- Count COG CATEGORY ----
                cog_cat = gene_to_cog_cat.get(gene, "No_EGGNOG")
                cog_cat_counts[cog_cat][flank] += 1
                all_cog_categories.add(cog_cat)

                # ---- Count COG ----
                cog = gene_to_cog.get(gene, "No_EGGNOG")
                cog_counts[cog][flank] += 1
                all_cogs.add(cog)

# ===============================================================
# STEP 3A — WRITE COG CATEGORY MATRIX
# ===============================================================
cog_cat_out = output_prefix + ".cog_category.txt"

with open(cog_cat_out, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    header = ["COG_category"] + all_flanks
    writer.writerow(header)

    for cog_cat in sorted(all_cog_categories):
        row = [cog_cat] + [cog_cat_counts[cog_cat].get(flank, 0) for flank in all_flanks]
        writer.writerow(row)

# ===============================================================
# STEP 3B — WRITE COG MATRIX
# ===============================================================
cog_out = output_prefix + ".cog.txt"

with open(cog_out, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    header = ["COG", "Description", "COG_CAT"] + all_flanks
    writer.writerow(header)

    for cog in sorted(all_cogs):
        if cog == "No_ANNOTATION" or cog == "No_EGGNOG" or cog == "No_PROTEIN":
            descr = "No_DESCRIPTION"
            cog_cat = "No_ANNOTATION"
        else:
            descr = cog_to_description.get(cog, "No_DESCRIPTION")
            cog_cat = cog_to_cog_cat.get(cog, "No_ANNOTATION")
            
        row = [cog, descr, cog_cat] + [cog_counts[cog].get(flank, 0) for flank in all_flanks]
        writer.writerow(row)

# ===============================================================
print(f"\nWrote COG matrix     → {cog_cat_out}")
print(f"Wrote TAXON matrix   → {cog_out}")
