#!/usr/bin/env python3

import sys
import pandas as pd

def parse_mgi(mgi_file):
    mgi_data = []
    with open(mgi_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3: continue
            contig = parts[0]
            start, end = map(int, parts[2].split('-'))
            mgi_data.append({'contig': contig, 'mgi_start': start, 'mgi_end': end})
    return pd.DataFrame(mgi_data)

def check_mgi_occupancy(mgi_df, cluster_file):
    clusters = pd.read_csv(cluster_file, sep='\t')
    results = []

    for _, mgi in mgi_df.iterrows():
        m_contig = mgi['contig']
        m_start = mgi['mgi_start']
        m_end = mgi['mgi_end']

        # Find clusters on the same contig with coordinate overlap
        # Logic: (StartA <= EndB) and (EndA >= StartB)
        match = clusters[
            (clusters['contig'].str.contains(m_contig)) & 
            (clusters['cluster_start'] <= m_end) & 
            (clusters['cluster_end'] >= m_start)
        ]

        results.append({
            'mgi_contig': m_contig,
            'mgi_range': f"{m_start}-{m_end}",
            'has_cluster': len(match) > 0,
            'cluster_count': len(match),
            'cluster_ids': ",".join(map(str, match['cluster_id'].tolist())) if not match.empty else "None"
        })

    return pd.DataFrame(results)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python mgi_occupancy_check.py <mgi_positions.txt> <transposase_clusters.txt>")
        sys.exit(1)

    mgi_df = parse_mgi(sys.argv[1])
    occupancy_results = check_mgi_occupancy(mgi_df, sys.argv[2])

    total_mgis = len(occupancy_results)
    occupied_mgis = occupancy_results['has_cluster'].sum()
    empty_mgis = total_mgis - occupied_mgis
    
    print("\n--- MGI Occupancy Summary ---")
    print(f"Total MGI intervals: {total_mgis}")
    print(f"MGIs with Transposase Clusters: {occupied_mgis} ({ (occupied_mgis/total_mgis)*100:.2f}%)")
    print(f"MGIs without Transposase Clusters: {empty_mgis} ({ (empty_mgis/total_mgis)*100:.2f}%)")
    
    # Save detailed report
    occupancy_results.to_csv("mgi_occupancy_report.csv", index=False)
    print("\nDetailed report saved to 'mgi_occupancy_report.csv'")
