#!/usr/bin/env python3

import pandas as pd
import sys
import re

# ------------------------------------------------------------
# Command line arguments
# ------------------------------------------------------------
if len(sys.argv) != 6:
    print("\nUsage:")
    print("  python analyze_pathogenic_flanks.py flanks.txt pathofact.tsv genes.gff summary_output.tsv intervals_output.tsv\n")
    sys.exit(1)

FLANKS = sys.argv[1]
PATHOFACT = sys.argv[2]
GFF = sys.argv[3]
OUT_SUMMARY = sys.argv[4]
OUT_INTERVALS = sys.argv[5]

# ------------------------------------------------------------
# Helper: Normalize contig names
# ------------------------------------------------------------
def normalize_contig(name):
    return re.sub(r'^bin_\d+_(.*)', r'\1', str(name))

# ------------------------------------------------------------
# 1. Load flanks
# ------------------------------------------------------------
flanks = pd.read_csv(FLANKS, sep="\t")
flanks["contig"] = flanks["contig"].apply(normalize_contig)

# ------------------------------------------------------------
# 2. Load PathoFact predictions
# ------------------------------------------------------------
patho = pd.read_csv(PATHOFACT, sep="\t")
patho["Toxin_prediction"] = patho["Toxin_prediction"].astype(str).str.lower()
patho["NormContig"] = patho["Contig"].apply(normalize_contig)

# Keep only pathogenic proteins
pathogenic = patho[patho["Toxin_prediction"] == "pathogenic"].copy()
# Unique key combining contig + ORF to avoid collisions
pathogenic["UniqueORF"] = pathogenic.apply(
    lambda row: f"{row['NormContig']}_{row['ORF']}", axis=1
)

# ------------------------------------------------------------
# 3. Parse GFF for protein coordinates
# ------------------------------------------------------------
gff_coords = {}
with open(GFF) as g:
    for line in g:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        contig, _, feature, start, end, *_ , info = parts
        if feature != "CDS":
            continue
        m = re.search(r"ID=([^;]+)", info)
        if not m:
            continue
        protein_id = m.group(1)
        norm_contig = normalize_contig(contig)
        unique_id = f"{norm_contig}_{protein_id}"  # unique key per contig
        gff_coords[unique_id] = {
            "contig": norm_contig,
            "start": int(start),
            "end": int(end)
        }

# ------------------------------------------------------------
# 4. Main analysis
# ------------------------------------------------------------
summary_rows = []
interval_rows = []

interval_size = 1000
flank_distance = 20000  # +/- 20 kb

for _, row in flanks.iterrows():
    cluster_id = row["cluster_id"]
    contig = normalize_contig(row["contig"])
    cstart = int(row["cluster_start"])
    cend = int(row["cluster_end"])

    # define flank region for summary (+/- 20kb)
    flank_start = max(0, cstart - flank_distance)
    flank_end = cend + flank_distance

    # create interval labels relative to cluster
    intervals = {}

    # upstream intervals (-20000 ... -1000), smaller number first
    for i in range(flank_distance // interval_size, 0, -1):
        start_label = (i - 1) * interval_size
        end_label = i * interval_size
        label = f"-{start_label}_{end_label}"
        intervals[label] = []

    # downstream intervals (+0 ... +20000)
    for i in range(flank_distance // interval_size + 1):
        start_label = i * interval_size
        end_label = (i + 1) * interval_size
        label = f"+{start_label}_{end_label}"
        intervals[label] = []

    # pathogenic proteins in this contig
    contig_pathogenic = pathogenic[pathogenic["NormContig"] == contig]

    pathogenic_hit_count = 0

    for _, prow in contig_pathogenic.iterrows():
        pid = prow["UniqueORF"]  # unique key
        if pid not in gff_coords:
            continue
        pstart = gff_coords[pid]["start"]
        pend = gff_coords[pid]["end"]

        # Only count proteins in summary if within +/- 20 kb of cluster
        if pend >= flank_start and pstart <= flank_end:
            pathogenic_hit_count += 1

        # Assign proteins to intervals
        for label in intervals:
            interval_start_label, interval_end_label = map(int, label.replace("+","").replace("-","").split("_"))

            if label.startswith("-"):
                # upstream intervals
                interval_start = cstart - interval_end_label
                interval_end = interval_start + interval_size - 1
            else:
                # downstream intervals
                interval_start = cend + interval_start_label + 1
                interval_end = cend + interval_end_label

            if pend >= interval_start and pstart <= interval_end:
                intervals[label].append(pid)

    # build summary row
    summary_rows.append({
        "cluster_id": cluster_id,
        "contig": contig,
        "pathogenic_hits": pathogenic_hit_count
    })

    # build interval row
    interval_row = {k: ",".join(v) if v else "-" for k,v in intervals.items()}
    interval_rows.append({
        "cluster_id": cluster_id,
        "contig": contig,
        "cluster_start": cstart,
        "cluster_end": cend,
        **interval_row
    })

# ------------------------------------------------------------
# 5. Write output files
# ------------------------------------------------------------
pd.DataFrame(summary_rows).to_csv(OUT_SUMMARY, sep="\t", index=False)
pd.DataFrame(interval_rows).to_csv(OUT_INTERVALS, sep="\t", index=False)

print("\n✔ Output written:")
print(f"  {OUT_SUMMARY}")
print(f"  {OUT_INTERVALS}\n")
