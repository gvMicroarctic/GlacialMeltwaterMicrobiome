#!/usr/bin/env python3

import sys
from collections import defaultdict, Counter

if len(sys.argv) != 4:
    print(f"Usage: {sys.argv[0]} <genes.gff> <hits.m8> <prefix>")
    sys.exit(1)

gff_file = sys.argv[1]
m8_file = sys.argv[2]
prefix = sys.argv[3]

all_genes = defaultdict(list)
transposase_hits = set()

# ---- STEP 1: read GFF ----
with open(gff_file) as g:
    for line in g:
        if line.startswith("#"):
            continue
        cols = line.strip().split("\t")
        if len(cols) < 9 or cols[2] != "CDS":
            continue

        contig = cols[0]
        start = int(cols[3])
        end = int(cols[4])

        attributes = cols[8]
        protein = None
        for field in attributes.split(";"):
            if field.startswith("ID=") or field.startswith("protein_id="):
                protein = field.split("=")[1]
                break
        if not protein:
            continue

        all_genes[contig].append((start, end, protein))

for contig in all_genes:
    all_genes[contig].sort(key=lambda x: x[0])

# ---- STEP 2: read m8 ----
with open(m8_file) as m:
    for line in m:
        if line.strip():
            protein = line.split("\t")[0]
            transposase_hits.add(protein)

# ---- STEP 3: detect contiguous clusters ----
cluster_info = []
for contig, genes in all_genes.items():
    current_cluster = []
    
    for (start, end, protein) in genes:
        if protein in transposase_hits:
            current_cluster.append((start, end, protein))
        else:
            if current_cluster:
                s = current_cluster[0][0]
                e = current_cluster[-1][1]
                prots = [p for _,_,p in current_cluster]
                cluster_info.append((contig, s, e, prots))
                current_cluster = []
    if current_cluster:
        s = current_cluster[0][0]
        e = current_cluster[-1][1]
        prots = [p for _,_,p in current_cluster]
        cluster_info.append((contig, s, e, prots))

# ---- STEP 4: write clusters file ----
clusters_file = f"{prefix}_contigs.txt"
with open(clusters_file, "w") as out:
    out.write("cluster_id\tcontig\tstart\tend\tn_transposases\tproteins\n")
    for i, (contig, s, e, prots) in enumerate(cluster_info, 1):
        out.write(f"{i}\t{contig}\t{s}\t{e}\t{len(prots)}\t{','.join(prots)}\n")

# ---- STEP 5: pivot per bin ----
pivot = defaultdict(Counter)

for (_, contig, _, _, prots) in [(i+1, *info) for i, info in enumerate(cluster_info)]:
    # bin = first two parts separated by "_" (can adjust if needed)
    bin_id = "_".join(contig.split("_")[:2])
    pivot[bin_id][len(prots)] += 1

# find all cluster sizes
sizes = sorted({size for p in pivot.values() for size in p})

pivot_file = f"{prefix}_mags.txt"
with open(pivot_file, "w") as out:
    out.write("bin\t" + "\t".join(map(str, sizes)) + "\n")
    for bin_id in sorted(pivot.keys()):
        row = [str(pivot[bin_id].get(size, 0)) for size in sizes]
        out.write(f"{bin_id}\t" + "\t".join(row) + "\n")

print(f"Wrote: {clusters_file}")
print(f"Wrote: {pivot_file}")
