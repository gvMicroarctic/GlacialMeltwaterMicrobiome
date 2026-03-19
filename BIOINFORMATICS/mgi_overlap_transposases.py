#!/usr/bin/env python3

import sys
import pandas as pd


def parse_mgi(mgi_file):
    mgi_data = []

    with open(mgi_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')

            if len(parts) < 3:
                continue

            contig = parts[0]

            # Extract start and end from "1-164" format
            start, end = map(int, parts[2].split('-'))

            mgi_data.append({
                'contig': contig,
                'start': start,
                'end': end
            })

    return pd.DataFrame(mgi_data)


def check_clusters(mgi_df, cluster_file):

    clusters = pd.read_csv(cluster_file, sep='\t')
    results = []

    for _, row in clusters.iterrows():

        # Match contig names (handling prefixes if necessary)
        # Note: MGI file uses 'NODE_1_bin_112'
        # Cluster uses 'bin_112_NODE_1_bin_112'

        c_contig = row['contig']
        c_start = row['cluster_start']
        c_end = row['cluster_end']
        c_id = row['cluster_id']

        # Look for MGIs on the same contig where intervals overlap
        # Overlap logic: (StartA <= EndB) and (EndA >= StartB)

        match = mgi_df[
            (mgi_df['contig'].apply(lambda x: x in c_contig)) &
            (mgi_df['start'] <= c_end) &
            (mgi_df['end'] >= c_start)
        ]

        if not match.empty:
            results.append({
                'cluster_id': c_id,
                'contig': c_contig,
                'overlap': True,
                'mgi_count': len(match)
            })
        else:
            results.append({
                'cluster_id': c_id,
                'contig': c_contig,
                'overlap': False,
                'mgi_count': 0
            })

    return pd.DataFrame(results)


if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: python check_overlap.py <mgi_positions.txt> <transposase_clusters.txt>")
        sys.exit(1)

    mgi_file = sys.argv[1]
    cluster_file = sys.argv[2]

    mgi_df = parse_mgi(mgi_file)
    results_df = check_clusters(mgi_df, cluster_file)

    print("\n--- Overlap Analysis Results ---")
    print(results_df.to_string(index=False))

    total_overlapping = results_df['overlap'].sum()

    print(f"\nTotal clusters found in MGI intervals: {total_overlapping} out of {len(results_df)}")
