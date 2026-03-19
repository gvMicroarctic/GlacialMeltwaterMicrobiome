#!/usr/bin/env python3

import sys
import csv
from collections import defaultdict

def main(intervals_file, toxins_file, out_name, out_desc):

    # -----------------------------
    # 1. Read PathoFact annotations
    # -----------------------------
    orf_to_name = {}
    orf_to_desc = {}

    with open(toxins_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cid = row["Composite_ORF_ID"]
            orf_to_name[cid] = row["NAME"]
            orf_to_desc[cid] = row["Description"]

    # ------------------------------------
    # 2. Read interval file and count hits
    # ------------------------------------
    name_counts = defaultdict(lambda: defaultdict(int))
    desc_counts = defaultdict(lambda: defaultdict(int))

    with open(intervals_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        interval_cols = reader.fieldnames[4:]  # after cluster_start/end

        for row in reader:
            for interval in interval_cols:
                cid = row[interval]
                if cid == "-" or cid == "":
                    continue

                if cid not in orf_to_name:
                    continue

                name = orf_to_name[cid]
                desc = orf_to_desc[cid]

                name_counts[name][interval] += 1
                desc_counts[desc][interval] += 1

    # --------------------
    # 3. Write output files
    # --------------------
    def write_output(outfile, counts_dict, first_col_name):
        with open(outfile, "w", newline="") as out:
            writer = csv.writer(out, delimiter="\t")
            writer.writerow([first_col_name] + interval_cols)

            for key in sorted(counts_dict):
                row = [key] + [counts_dict[key].get(i, 0) for i in interval_cols]
                writer.writerow(row)

    write_output(out_name, name_counts, "NAME")
    write_output(out_desc, desc_counts, "Description")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(
            "Usage:\n"
            "  count_toxins_by_interval.py "
            "<flank_toxin_intervals.tsv> "
            "<PathoFact_toxins.tsv> "
            "<output_NAME.tsv> "
            "<output_Description.tsv>"
        )
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
