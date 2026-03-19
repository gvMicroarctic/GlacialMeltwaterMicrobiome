#!/usr/bin/env python
import csv
import sys
import os # Import os for output file path handling

# Initialize a single dictionary to hold results from both files
# { contig : set([gene1, gene2, ...]) }
combined_hallmarks = {}
output_file_path = "contig_hallmark_summary.tsv" # Define the output file name

# --- Step 1: genomad file ---
print("Processing Genomad Hallmarks (File 1).")
file1_path = "./viral_contig/AMG/genomad/viral_contigs_unique_final_5000_clean_summary/viral_contigs_unique_final_5000_clean_virus_genes.tsv"
#hallmark genes if 1 in virus_hallmark column

with open(file1_path) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if row.get("virus_hallmark") == "1":
            gene_full = row["gene"]  # e.g., contig_6_1
            contig = "_".join(gene_full.split("_")[:2])  # contig_6
            
            # Add gene to the combined dictionary
            if contig not in combined_hallmarks:
                combined_hallmarks[contig] = set()
            combined_hallmarks[contig].add(gene_full)

# --- Step 2: virsorter2 file ---
print("Processing VirSorter2 Hallmarks (File 2).")
file2_path = "./viral_contig/AMG/virsorter2/for-dramv/viral-affi-contigs-for-dramv.tab"
#hallmark genes if 0 in fourth-to-last value

with open(file2_path) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith(">"):
            continue

        parts = line.split("|")
        
        if len(parts) >= 3: 
            third_last = parts[-4]
            if third_last == "0":  # hallmark found
                gene_full = parts[0]  # contig_2__1
                contig = gene_full.split("__")[0]  # contig_2
                gene_mapped = gene_full.replace("__", "_")  # contig_2_1
                
                # Add gene to the combined dictionary
                if contig not in combined_hallmarks:
                    combined_hallmarks[contig] = set()
                combined_hallmarks[contig].add(gene_mapped)

# --- Final Step: Write Combined Results to a TSV File ---
# Use 'w' mode to create or overwrite the file
with open(output_file_path, 'w', newline='') as outfile:
    # Use csv.writer to easily handle tab separation
    writer = csv.writer(outfile, delimiter='\t')
    
    # Write the header row
    writer.writerow(["contig_name", "hallmark_gene_count", "hallmark_gene_list"])
    
    contigs_written = 0
    
    # Iterate over the combined results
    for contig, genes_set in combined_hallmarks.items():
        gene_count = len(genes_set)
        
        # Convert the set of genes into a sorted, semicolon-separated string
        gene_list_str = ";".join(sorted(list(genes_set)))
        
        # Write the row to the file
        writer.writerow([contig, gene_count, gene_list_str])
        contigs_written += 1

print(f"Wrote summary for {contigs_written} contigs to {output_file_path}.")