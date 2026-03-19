#!/usr/bin/env python3
import pandas as pd
import argparse
import sys

def main(annotation_file, gff_file, outprefix):
    # -----------------------------
    # Read annotation table
    # -----------------------------
    try:
        df = pd.read_csv(annotation_file, sep="\t", header=0)
    except Exception as e:
        print(f"Error reading annotation file: {e}")
        sys.exit(1)
    
    # -----------------------------
    # Count total genes in GFF
    # -----------------------------
    # Note: We filter for 'gene' or 'CDS' to avoid overcounting 
    # every feature (exons, mRNAs, etc.)
    total_genes = 0
    try:
        with open(gff_file) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split("\t")
                if len(parts) > 2 and parts[2].lower() == "gene":
                    total_genes += 1
        
        # Fallback: if no 'gene' features found, count 'CDS'
        if total_genes == 0:
            with open(gff_file) as f:
                total_genes = sum(1 for line in f if "\tCDS\t" in line)
    except Exception as e:
        print(f"Error reading GFF file: {e}")
        sys.exit(1)

    if total_genes == 0:
        print("Warning: No genes/CDS found in GFF. Relative abundance will be NaN.")
        total_genes = float('nan')
    
    print(f"Total genes identified in GFF: {total_genes}")
    
    # -----------------------------
    # Helper function for counting
    # -----------------------------
    def get_counts(dataframe, column_name):
        # count occurrences
        counts = dataframe[column_name].value_counts().reset_index()
        # Rename columns (handles different pandas version behaviors)
        counts.columns = [column_name, "abs_abundance"]
        
        # Calculate relative abundance
        counts["rel_abundance"] = counts["abs_abundance"] / total_genes
        return counts

    # Process NAME
    if "NAME" in df.columns:
        name_counts = get_counts(df, "NAME")
        name_counts.to_csv(f"{outprefix}_by_NAME.tsv", sep="\t", index=False)
    
    # Process Description
    if "Description" in df.columns:
        desc_counts = get_counts(df, "Description")
        desc_counts.to_csv(f"{outprefix}_by_Description.tsv", sep="\t", index=False)
    
    print("Done! Output files created successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count toxins and calculate relative abundance.")
    parser.add_argument("-a", "--annotation", required=True, help="PathoFact toxins TSV file")
    parser.add_argument("-g", "--gff", required=True, help="GFF file with all genes")
    parser.add_argument("-o", "--outprefix", required=True, help="Prefix for output files")
    
    args = parser.parse_args()
    main(args.annotation, args.gff, args.outprefix)