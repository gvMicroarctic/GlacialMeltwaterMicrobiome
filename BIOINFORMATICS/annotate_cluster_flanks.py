#!/usr/bin/env python3

import sys
from collections import defaultdict

if len(sys.argv) != 4:
    print(f"Usage: {sys.argv[0]} <clusters.tsv> <genes.gff> <prefix>")
    sys.exit(1)

clusters_file = sys.argv[1]
gff_file = sys.argv[2]
prefix = sys.argv[3]

WINDOWS = [0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000]

# ---- read cluster file ----
clusters = []
with open(clusters_file) as f:
    header = next(f)
    for line in f:
        parts = line.strip().split("\t")
        cluster_id, contig, start, end = parts[:4]
        clusters.append((int(cluster_id), contig, int(start), int(end), line.strip()))

# ---- read GFF ----
genes = defaultdict(list)

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

        # extract protein
        protein = None
        for field in cols[8].split(";"):
            if field.startswith("ID=") or field.startswith("protein_id="):
                protein = field.split("=")[1]
                break
        if not protein:
            continue

        genes[contig].append((start, end, protein))

# ---- annotate flanking regions ----

outfile = f"{prefix}.txt"
with open(outfile, "w") as out:
    flank_cols = []

    # build header columns in order
    for i in range(len(WINDOWS)-1):
        flank_cols.append(f"-{WINDOWS[i]}_{WINDOWS[i+1]}")
    for i in range(len(WINDOWS)-1):
        flank_cols.append(f"+{WINDOWS[i]}_{WINDOWS[i+1]}")

    out.write("cluster_id\tcontig\tcluster_start\tcluster_end\ttransposase_number\ttransposases\t" + "\t".join(flank_cols) + "\n")

    for cluster_id, contig, cstart, cend, cline in clusters:
        flank_hits = {col: [] for col in flank_cols}

        # upstream negative windows → use protein END
        for i in range(len(WINDOWS)-1):
            w1, w2 = WINDOWS[i], WINDOWS[i+1]
            window_min = cstart - w2
            window_max = cstart - w1
            label = f"-{w1}_{w2}"
            
            for gstart, gend, prot in genes[contig]:
                #print(f"out {cluster_id} and {contig} and {gstart} and {gend} and {window_min} and {window_max} and {prot}")
                
                if not (gend < window_min or gstart > window_max):
                    #print(f"in {window_min} and {window_max}")
                    flank_hits[label].append(prot)

        # downstream positive windows → use protein START
        for i in range(len(WINDOWS)-1):
            w1, w2 = WINDOWS[i], WINDOWS[i+1]
            window_min = cend + w1
            window_max = cend + w2
            label = f"+{w1}_{w2}"

            for gstart, gend, prot in genes[contig]:
                if not (gend < window_min or gstart > window_max):
                    flank_hits[label].append(prot)

        # build output line
        row = [cline] + [",".join(flank_hits[col]) if flank_hits[col] else "-" for col in flank_cols]
        out.write("\t".join(row) + "\n")

print(f"flank annotation written to: {outfile}")
