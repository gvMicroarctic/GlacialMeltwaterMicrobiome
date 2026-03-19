#!/usr/bin/env python3

import sys
import csv

PRIORITY = {
    "PFAM": 1,
    "KEGG": 2,
    "TIGR": 3,
    "SwissProt": 4
}

def main(infile, outfile):
    best_hit = {}

    with open(infile, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            comp_id = row["ORF_ID"]
            db = row["Database"]

            # skip unknown databases
            if db not in PRIORITY:
                continue

            if comp_id not in best_hit:
                best_hit[comp_id] = row
            else:
                # lower number = higher priority
                if PRIORITY[db] < PRIORITY[best_hit[comp_id]["Database"]]:
                    best_hit[comp_id] = row

    with open(outfile, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(best_hit.values())


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: deduplicate_pathofact.py <input.tsv> <output.tsv>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
