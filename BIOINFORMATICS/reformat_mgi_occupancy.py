#!/usr/bin/env python3

import argparse
import pandas as pd


def load_mag_map(file):
    mag_map = {}

    with open(file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            bin_id, mag = line.split()
            mag_map[bin_id] = mag

    return mag_map


def main():

    parser = argparse.ArgumentParser(description="Reformat MGI table with MAG and contig names.")
    parser.add_argument("-i", required=True, help="Input mgi_occupancy_report.csv")
    parser.add_argument("-m", required=True, help="MAG correspondence file")
    parser.add_argument("-o", required=True, help="Output file")

    args = parser.parse_args()

    df = pd.read_csv(args.i)

    mag_map = load_mag_map(args.m)

    new_rows = []

    for _, row in df.iterrows():

        contig_full = row["mgi_contig"]

        parts = contig_full.split("_")
        contig = f"{parts[0]}_{parts[1]}"      # NODE_1
        bin_id = f"bin_{parts[-1]}"            # bin_112

        mag = mag_map.get(bin_id, "NA")

        new_rows.append({
            "MAG": mag,
            "contig": contig,
            "mgi_range": row["mgi_range"],
            "has_cluster": row["has_cluster"],
            "cluster_count": row["cluster_count"],
            "cluster_ids": row["cluster_ids"]
        })

    out_df = pd.DataFrame(new_rows)

    out_df.to_csv(args.o, index=False)


if __name__ == "__main__":
    main()
