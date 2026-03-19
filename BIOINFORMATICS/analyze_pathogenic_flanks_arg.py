#!/usr/bin/env python3

import pandas as pd
import sys

# ------------------------------------------------------------
# Command line arguments
# ------------------------------------------------------------
if len(sys.argv) != 5:
    print("\nUsage:")
    print("  python analyze_arg_flanks.py flanks.txt merged_arg_table.tsv summary_output.tsv intervals_output.tsv\n")
    sys.exit(1)

FLANKS = sys.argv[1]
MERGED_ARG = sys.argv[2]
OUT_SUMMARY = sys.argv[3]
OUT_INTERVALS = sys.argv[4]

# ------------------------------------------------------------
# 1. Load flanks
# ------------------------------------------------------------
flanks = pd.read_csv(FLANKS, sep="\t")

print(f"{flanks}")

# Explicit contig formatting inline
# 'bin_100_NODE_184_bin_100_2' -> 'NODE_184_bin_100_2'
flanks["NormContig"] = ["_".join(str(c).split("_")[2:]) for c in flanks["contig"]]

print(f"{flanks['NormContig']}")

# ------------------------------------------------------------
# 2. Load merged ARG table
# ------------------------------------------------------------
merged = pd.read_csv(MERGED_ARG, sep="\t")

print(f"{merged}")

# Extract start and end from new_protein_name
# Assumes format: ..._<start>_<end>
merged["start"] = merged["new_protein_name"].apply(lambda x: int(str(x).split("_")[-2]))
merged["end"]   = merged["new_protein_name"].apply(lambda x: int(str(x).split("_")[-1]))

# Explicit contig formatting inline
merged["NormContig"] = ["_".join(str(c).split("_")[:-2]) for c in merged["new_protein_name"]]

print(f"yo {merged['NormContig']}")

# Unique protein ID for tracking
merged["UniqueID"] = merged["new_protein_name"]

# ------------------------------------------------------------
# 3. Main analysis
# ------------------------------------------------------------
summary_rows = []
interval_rows = []

interval_size = 1000
flank_distance = 20000  # +/- 20 kb

for _, row in flanks.iterrows():
    cluster_id = row["cluster_id"]
    # Explicit contig formatting inline
    contig_parts = str(row["contig"]).split("_")
    if len(contig_parts) > 2:
        contig = "_".join(contig_parts[2:])
    else:
        contig = str(row["contig"])
    cstart = int(row["cluster_start"])
    cend = int(row["cluster_end"])

    # Define flank region
    flank_start = max(0, cstart - flank_distance)
    flank_end = cend + flank_distance

    # Create interval labels
    intervals = {}
    for i in range(flank_distance // interval_size, 0, -1):
        label = f"-{(i-1)*interval_size}_{i*interval_size}"
        intervals[label] = []
    for i in range(flank_distance // interval_size + 1):
        label = f"+{i*interval_size}_{(i+1)*interval_size}"
        intervals[label] = []

    # ARG proteins in this contig
    contig_args = merged[merged["NormContig"] == contig]

    hit_count = 0
    for _, prow in contig_args.iterrows():
        pstart = int(prow["start"])
        pend = int(prow["end"])
        pid = prow["UniqueID"]

        # Summary: count proteins within +/- 20 kb
        if pend >= flank_start and pstart <= flank_end:
            hit_count += 1

        # Assign to intervals
        for label in intervals:
            interval_start_label, interval_end_label = map(int, label.replace("+","").replace("-","").split("_"))
            if label.startswith("-"):
                interval_start = cstart - interval_end_label
                interval_end = interval_start + interval_size - 1
            else:
                interval_start = cend + interval_start_label + 1
                interval_end = cend + interval_end_label

            if pend >= interval_start and pstart <= interval_end:
                intervals[label].append(pid)

    # Build summary row
    summary_rows.append({
        "cluster_id": cluster_id,
        "contig": contig,
        "arg_hits": hit_count
    })

    # Build interval row
    interval_row = {k: ",".join(v) if v else "-" for k, v in intervals.items()}
    interval_rows.append({
        "cluster_id": cluster_id,
        "contig": contig,
        "cluster_start": cstart,
        "cluster_end": cend,
        **interval_row
    })

# ------------------------------------------------------------
# 4. Write output
# ------------------------------------------------------------
pd.DataFrame(summary_rows).to_csv(OUT_SUMMARY, sep="\t", index=False)
pd.DataFrame(interval_rows).to_csv(OUT_INTERVALS, sep="\t", index=False)

print("\n✔ Output written:")
print(f"  Summary: {OUT_SUMMARY}")
print(f"  Intervals: {OUT_INTERVALS}\n")

# ------------------------------------------------------------
# 5. Count proteins in +/-10 kb and +/-20 kb for reporting
# ------------------------------------------------------------
proteins_in_10kb = set()
proteins_in_20kb = set()

for _, row in flanks.iterrows():
    cluster_id = row["cluster_id"]
    contig_parts = str(row["contig"]).split("_")
    if len(contig_parts) > 2:
        contig = "_".join(contig_parts[2:])
    else:
        contig = str(row["contig"])
    cstart = int(row["cluster_start"])
    cend = int(row["cluster_end"])

    flank_10_start = max(0, cstart - 10000)
    flank_10_end   = cend + 10000
    flank_20_start = max(0, cstart - 20000)
    flank_20_end   = cend + 20000

    contig_args = merged[merged["NormContig"] == contig]

    for _, prow in contig_args.iterrows():
        pstart = int(prow["start"])
        pend   = int(prow["end"])
        pid    = prow["UniqueID"]
        
        if pend >= flank_10_start and pstart <= flank_10_end:
            proteins_in_10kb.add(pid)
        if pend >= flank_20_start and pstart <= flank_20_end:
            proteins_in_20kb.add(pid)

print(f"\n✔ Total proteins within ±10 kb of any cluster: {len(proteins_in_10kb)}")
print(f"  Proteins: {sorted(proteins_in_10kb)}\n")
print(f"✔ Total proteins within ±20 kb of any cluster: {len(proteins_in_20kb)}")

# ------------------------------------------------------------
# 6. Extract rows for ±10 kb proteins from merged ARG table
# ------------------------------------------------------------
proteins_10kb_list = list(proteins_in_10kb)
proteins_10kb_rows = merged[merged["UniqueID"].isin(proteins_10kb_list)]

# Print results to terminal
print(f"\n✔ Rows for ±10 kb proteins ({len(proteins_10kb_rows)} rows):")
print(proteins_10kb_rows)

# Optionally, save to TSV
output_10kb_file = "args_within_10kb.tsv"
proteins_10kb_rows.to_csv(output_10kb_file, sep="\t", index=False)
print(f"\n✔ Saved ±10 kb protein rows to {output_10kb_file}")
