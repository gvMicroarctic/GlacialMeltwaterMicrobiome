#!/usr/bin/env python3

import glob
import os
import re

######################################
# 1. Read GFF files → build lookup
######################################

gff_pattern = (
    "/storage/varliero/icevirome_prj/MAG/pathofact/"
    "PathoFact_MAGs_*/PathoFact_intermediate/Prodigal/bin_*.gff"
)

gff_files = glob.glob(gff_pattern)

# (bin_name, orf_id) -> composite_id
gff_lookup = {}

for gff in gff_files:
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            cols = line.rstrip().split("\t")
            if len(cols) < 9:
                continue

            # get contig name
            contig_name = cols[0]

            # get protein (ORF) name
            m = re.search(r"ID=([^;]+)", cols[8])
            if not m:
                continue
            orf_id = m.group(1)          # e.g. 1_1
            
            # get bin name
            bin_num = contig_name.split("_")[3]
            bin_name = "bin_" + bin_num
            
            # create composite name
            composite = f"{contig_name}_{orf_id}"

            gff_lookup[(bin_name, orf_id)] = composite

print(f"Loaded {len(gff_lookup)} ORFs from GFF files")

######################################
# 2. Merge PathoFact TSVs + add column
######################################

tsv_pattern = (
    "/storage/varliero/icevirome_prj/MAG/pathofact/"
    "PathoFact_MAGs_*/PathoFact_report/"
    "Toxin_gene_library_bin_*_report.tsv"
)

outfile = "/storage/varliero/icevirome_prj/MAG/pathofact/PathoFact_toxins.tsv"
tsv_files = sorted(glob.glob(tsv_pattern))

if not tsv_files:
    print("No prediction files found.")
    exit(1)

print(f"Found {len(tsv_files)} TSV files")

with open(outfile, "w") as out:
    header_written = False

    for f in tsv_files:
        # extract bin name from filename
        # e.g. Toxin_gene_library_bin_k143_2149083_bin_19_report.tsv
        base = os.path.basename(f)
        bin_num = base.replace("Toxin_gene_library_bin_", "").replace("_report.tsv", "")
        bin_name = "bin_" + bin_num

        with open(f) as infile:
            for i, line in enumerate(infile):
                line = line.rstrip()

                if i == 0:
                    if not header_written:
                        out.write(line + "\tComposite_ORF_ID\n")
                        header_written = True
                    continue

                fields = line.split("\t")

                 # get protein (ORF) name
                orf = fields[1]

                composite = gff_lookup.get(
                    (bin_name, orf),
                    "NA"
                )h

                out.write(line + "\t" + composite + "\n")

print("Done.")
