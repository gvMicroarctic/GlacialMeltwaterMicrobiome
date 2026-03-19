#!/usr/bin/env python3
import csv
import argparse
from collections import defaultdict

# --------------------------------------------
# Parse arguments
# --------------------------------------------
parser = argparse.ArgumentParser(description="Detect megaclusters with min/max span thresholds.")
parser.add_argument("input", help="Input cluster file")
parser.add_argument("output", help="Output megacluster file")
parser.add_argument("minspan", type=int, help="Minimum required span (e.g. 2000)")
parser.add_argument("maxspan", type=int, help="Maximum allowed span (e.g. 10000)")
args = parser.parse_args()

MIN_SPAN = args.minspan
MAX_SPAN = args.maxspan

# --------------------------------------------
# Load clusters grouped by contig
# --------------------------------------------
clusters_by_contig = defaultdict(list)

with open(args.input, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        row["start"] = int(row["start"])
        row["end"] = int(row["end"])
        row["n_transposases"] = int(row["n_transposases"])
        row["proteins"] = row["proteins"].split(",")
        clusters_by_contig[row["contig"]].append(row)

# --------------------------------------------
# Build megaclusters with min/max span rules
# --------------------------------------------
megaclusters = []
megacluster_id = 1

for contig, rows in clusters_by_contig.items():

    rows.sort(key=lambda r: r["start"])

    # Start first megacluster
    current = {
        "contig": contig,
        "start": rows[0]["start"],
        "end": rows[0]["end"],
        "n_transposases": rows[0]["n_transposases"],
        "proteins": list(rows[0]["proteins"]),
        "count": 1
    }

    for r in rows[1:]:

        new_end = max(current["end"], r["end"])
        span = new_end - current["start"]

        if span <= MAX_SPAN:
            # merge
            current["end"] = new_end
            current["n_transposases"] += r["n_transposases"]
            current["proteins"].extend(r["proteins"])
            current["count"] += 1
        else:
            # finalize only if meets min requirements
            final_span = current["end"] - current["start"]
            if current["count"] >= 2 and final_span >= MIN_SPAN:
                megaclusters.append((megacluster_id, current))
                megacluster_id += 1

            # start a new megacluster
            current = {
                "contig": contig,
                "start": r["start"],
                "end": r["end"],
                "n_transposases": r["n_transposases"],
                "proteins": list(r["proteins"]),
                "count": 1
            }

    # finalize last megacluster
    final_span = current["end"] - current["start"]
    if current["count"] >= 2 and final_span >= MIN_SPAN:
        megaclusters.append((megacluster_id, current))
        megacluster_id += 1

# --------------------------------------------
# Write output
# --------------------------------------------
with open(args.output, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "megacluster_id", "contig", "start", "end",
        "span", "n_transposases", "n_clusters", "proteins"
    ])

    for mcid, mc in megaclusters:
        span = mc["end"] - mc["start"]
        writer.writerow([
            mcid,
            mc["contig"],
            mc["start"],
            mc["end"],
            span,
            mc["n_transposases"],
            mc["count"],
            ",".join(mc["proteins"])
        ])

print(f"Saved {len(megaclusters)} megaclusters to {args.output}")
