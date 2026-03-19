#!/usr/bin/env python
import csv
import re
import sys
import os

if len(sys.argv) != 3:
    # Print usage information if the wrong number of arguments is provided
    print(f"Usage: python {sys.argv[0]} <input_tsv_file> <output_txt_file>", file=sys.stderr)
    print("\nExample:", file=sys.stderr)
    print(f"python {sys.argv[0]} contig_hallmark_summary.tsv intervening_genes.txt", file=sys.stderr)
    sys.exit(1)

input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

# List to store ALL intervening genes found across ALL contigs
all_intervening_genes = []

## 1. Process Input File (Read from Command Line Argument)
print(f"## Processing input file: {input_file_path}")

try:
    with open(input_file_path, 'r', newline='') as infile:
        # Use csv.DictReader for easier access by column name
        reader = csv.DictReader(infile, delimiter='\t')
        
        for row in reader:
            # Basic header check (expanded check for robustness)
            required_headers = ["contig_name", "hallmark_gene_count", "hallmark_gene_list"]
            if not all(k in row for k in required_headers):
                print(f"Error: Input file must contain headers: {', '.join(required_headers)}", file=sys.stderr)
                sys.exit(1)

            contig = row["contig_name"]
            
            try:
                gene_count = int(row["hallmark_gene_count"])
            except ValueError:
                print(f"Warning: Skipping contig {contig} due to non-integer gene count.", file=sys.stderr)
                continue
            
            # Check the condition: only proceed if more than one hallmark gene is present
            if gene_count > 1:
                hallmark_genes = set(row["hallmark_gene_list"].split(';'))
                
                # --- Extract Gene Indices and Range ---
                indices = []
                contig_prefix = None
                
                # Extract the numeric index (Y) from the end of the gene name (contig_X_Y)
                for gene in hallmark_genes:
                    match = re.search(r'_(\d+)$', gene)
                    if match:
                        indices.append(int(match.group(1)))
                        # Determine the contig prefix (e.g., 'contig_30') once
                        if contig_prefix is None:
                            contig_prefix = "_".join(gene.split('_')[:-1])

                if not indices:
                    print(f"Warning: Could not extract indices for contig {contig}. Skipping.", file=sys.stderr)
                    continue

                min_index = min(indices)
                max_index = max(indices)
                
                # --- Identify Missing Genes ---
                current_missing_genes = []
                
                # Iterate through all possible gene indices in the range [min_index, max_index]
                for i in range(min_index, max_index + 1):
                    potential_gene = f"{contig_prefix}_{i}"
                    
                    # Check if this potential gene is NOT in the hallmark_genes set
                    if potential_gene not in hallmark_genes:
                        current_missing_genes.append(potential_gene)
                
                # Add the missing genes for this contig to the master list
                if current_missing_genes:
                     all_intervening_genes.extend(current_missing_genes)

except FileNotFoundError:
    print(f"Error: Input file not found at {input_file_path}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"An unexpected error occurred while reading the input file: {e}", file=sys.stderr)
    sys.exit(1)

## 2. Write Output File (Plain Text, Single Column)

print(f"Writing gene list to output file: {output_file_path}")

try:
    with open(output_file_path, 'w') as outfile:
        # Write each gene name followed by a newline character
        outfile.write('\n'.join(all_intervening_genes))

    print(f"Wrote {len(all_intervening_genes)} intervening gene names to {output_file_path}.")

except Exception as e:
    print(f"An unexpected error occurred while writing the file: {e}", file=sys.stderr)
    sys.exit(1)
