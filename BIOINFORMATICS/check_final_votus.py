#!/usr/bin/env python3

# Clean vOTUs from transposons, plasmids, and add hallmark info

votu = {}

# 1. Transposase file (dram)
with open("./viral_contig/dramv/votu_transposase_dram.txt") as f1:
    for line in f1:
        contig_base0 = line.split("_")
        contig_base = contig_base0[0] + "_" + contig_base0[1]
        gene = line.split("\t")[9] + ";" + line.split("\t")[11] + ";" + line.split("\t")[23] + ";" + line.split("\t")[29]
        if contig_base not in votu:
            votu[contig_base] = {}
        votu[contig_base]["transposase_dram"] = gene
        
# 2. Transposase file (eggnog mapper)
with open("./viral_contig/COG/votu_transposase_cog.txt") as f1:
    for line in f1:
        contig_base0 = line.split("_")
        contig_base = contig_base0[0] + "_" + contig_base0[1]
        gene = line.split("\t")[7]
        if contig_base not in votu:
            votu[contig_base] = {}
        votu[contig_base]["transposase_cog"] = gene

# 3. Genomad file
with open("./viral_contig/genomad_plasmid/viral_contigs_unique_final_5000_summary/viral_contigs_unique_final_5000_plasmid_summary.tsv") as f2:
    header = f2.readline()  # skip header
    for line in f2:
        if not line.startswith("contig_"):
            continue
        contig_base = line.split("\t")[0]
        if contig_base not in votu:
            votu[contig_base] = {}
        votu[contig_base]["genomad"] = "1"

# 4. Plasme file
with open("./viral_contig/viral_contigs_unique_plasmids.txt_report.csv") as f3:
    header = f3.readline()  # skip header
    for line in f3:
        if not line.startswith("contig_"):
            continue
        contig_base = line.split("\t")[0]
        if contig_base not in votu:
            votu[contig_base] = {}
        votu[contig_base]["plasme"] = "1"

# 5. Hallmark file
with open("./viral_contig/virsorter2_2/final-viral-score.tsv") as f4:
    header = f4.readline()  # skip header
    for line in f4:
        if not line.startswith("contig_"):
            continue
        fields = line.strip().split("\t")
        contig_base = fields[0].split("||")[0]  # remove ||full suffix
        hallmark = fields[9]  # "hallmark" column (0-based index 10)
        if contig_base in votu:
            votu[contig_base]["hallmark"] = hallmark

# 6. Write output
with open("./viral_contig/votu_manual_curation.txt", "w") as out:
    out.write("contig\ttransposase_dram\ttransposase_cog\tgenomad\tplasme\thallmark\n")
    for contig in sorted(votu.keys()):
        trans_dram = votu[contig].get("transposase_dram", "0")
        trans_cog = votu[contig].get("transposase_cog", "0")
        geno = votu[contig].get("genomad", "0")
        plasm = votu[contig].get("plasme", "0")
        hallmark = votu[contig].get("hallmark", "NA")
        out.write(f"{contig}\t{trans_dram}\t{trans_cog}\t{geno}\t{plasm}\t{hallmark}\n")
