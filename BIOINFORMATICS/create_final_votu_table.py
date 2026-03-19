#!/usr/bin/env python3

import pandas as pd

# Define file paths
fasta_path = "./viral_contig/viral_contigs_unique_final_5000_clean.fasta"
checkv_path = "./viral_contig/checkv/quality_summary.tsv"
virsorter_path = "./viral_contig/AMG/virsorter2/final-viral-score.tsv"
taxonomy_path = "./viral_contig/votu_taxonomy_clean.txt"
bacphlip_path = "./viral_contig/bacphlip_results_checkv_100.txt"

# 1. Get sequence headers from FASTA
# This creates the master list of contigs
contig_ids = []
with open(fasta_path, 'r') as f:
    for line in f:
        if line.startswith(">"):
            contig_ids.append(line.strip().replace(">", "").split()[0])

df_final = pd.DataFrame({'contig_id': contig_ids})

# 2. Add CheckV info
cv_cols = ["contig_id", "contig_length", "gene_count", "viral_genes", "host_genes", 
           "checkv_quality", "completeness", "completeness_method", "contamination"]
df_cv = pd.read_csv(checkv_path, sep='\t')[cv_cols]
df_final = pd.merge(df_final, df_cv, on='contig_id', how='left')

# 3. Add VirSorter2 info
vs_cols = ["seqname", "max_score", "max_score_group", "hallmark"]
df_vs = pd.read_csv(virsorter_path, sep='\t')[vs_cols]
df_vs.columns = ['contig_id', 'max_score', 'max_score_group', 'hallmark']
df_final = pd.merge(df_final, df_vs, on='contig_id', how='left')

# 4. Add Taxonomy info (combining columns with ;)
df_tax = pd.read_csv(taxonomy_path, sep='\t', header=None)
# Join all columns from index 1 to the end with a semicolon
df_tax['taxonomy'] = df_tax.iloc[:, 1:].fillna('').apply(lambda x: ';'.join(x.astype(str)), axis=1)
df_tax = df_tax[[0, 'taxonomy']].rename(columns={0: 'contig_id'})
df_final = pd.merge(df_final, df_tax, on='contig_id', how='left')

# 5. Add Bacphlip info
df_bp = pd.read_csv(bacphlip_path, sep='\t', header=None, names=['contig_id', 'bacphlip'])
df_final = pd.merge(df_final, df_bp, on='contig_id', how='left')

# Save the final table
output_file = "./viral_contig/viral_contig_metadata_final.tsv"
df_final.to_csv(output_file, sep='\t', index=False)

print(f"Process complete. Table saved to: {output_file}")
