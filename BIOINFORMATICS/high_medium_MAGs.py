#!/usr/bin/env python3

import shutil
import argparse
import os
import re

#retrieve only medium and high quality MAGs

#set up command-line argument parsing
parser = argparse.ArgumentParser(description="Filter MAGs by completeness and contamination, and rename them.")
parser.add_argument("input_file", help="Path to CheckM quality report TSV file")
parser.add_argument("source_folder", help="Folder containing MAG .fa files")
parser.add_argument("output_folder", help="Folder to copy and rename high- and medium-quality MAGs into")
args = parser.parse_args()

#ensure output folder exists
os.makedirs(args.output_folder, exist_ok=True)

count_in = 0
with open(args.input_file, "r") as file_in:
    #skip header
    file_in.readline()
    for line in file_in:
        info = line.strip().split("\t")
        
        #CheckM report columns: [0] Bin ID, [1] Completeness, [2] Contamination
        bin_id = info[0]
        try:
            completeness = float(info[1])
            contamination = float(info[2])
        except (ValueError, IndexError):
            print(f"Skipping line with malformed data: {line.strip()}")
            continue

        if completeness >= 50 and contamination < 10:
            source_file = os.path.join(args.source_folder, bin_id + ".fa")
            print(f"Processing: {source_file}")
            
            if os.path.exists(source_file):
                #determine the new filename structure ---
                
                #split the original Bin ID (e.g., 'sample.001') to form the new filename
                base_name_parts = bin_id.split(".")
                
                #create the new file name: <base_name>.<number>.fa
                if len(base_name_parts) >= 2:
                    new_file_name = f"{base_name_parts[0]}_{base_name_parts[1]}.fa"
                    mag_name_for_header = f"{base_name_parts[0]}_{base_name_parts[1]}"
                else:
                    #fallback for unexpected bin IDs (e.g., if it's just 'sample')
                    new_file_name = f"{bin_id}_filtered.fa"
                    mag_name_for_header = bin_id

                dest_file = os.path.join(args.output_folder, new_file_name)
                
                #open and process the file line by line
                
                #open the destination file in the output folder for writing
                #we do not use shutil.copy() here because we are modifying the content.
                try:
                    with open(source_file, 'r') as in_file, open(dest_file, "w") as file_out:
                        for line in in_file:
                            line = line.strip()
                            if line.startswith('>'):
                                #header Modification Logic from your second script ---
                                
                                #remove anything after "_length" or "flag"
                                if "length" in line:
                                    # Example: >NODE_1_length_1000... -> >NODE_1
                                    header = line.split("_length")[0]
                                elif "flag" in line:
                                    #example: >NODE_2_flag_0_ ... -> >NODE_2
                                    header0 = line.split("flag")[0]
                                    header = header0.strip()
                                else:
                                    #if neither is present, use the whole header line
                                    header = line 
                                    
                                #new header format: >[OriginalHeader]_[MAGName]
                                contig_name = f'{header}_{mag_name_for_header}'
                                file_out.write(f'{contig_name}\n')
                            else:
                                #write sequence line as-is
                                file_out.write(f'{line}\n')
                                
                    count_in += 1
                
                except IOError as e:
                    print(f"Error reading/writing file {bin_id}: {e}")
            else:
                print(f"Warning: {source_file} does not exist!")

print(f"Number of MAGs with high- and medium-quality: {count_in}")
