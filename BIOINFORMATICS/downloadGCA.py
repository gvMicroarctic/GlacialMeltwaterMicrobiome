#!/usr/bin/env python3

#Download specific genomes from NCBI

import argparse
import subprocess
import os

#get variables from command line
parser = argparse.ArgumentParser(description='Download specific genomes from NCBI.')
parser.add_argument('-i', type=str, help='Input file.')
parser.add_argument('-o', type=str, help='Output folder.')
parser.add_argument('-l', type=str, help='Log output file.')
args = parser.parse_args()

#name variables from command line
input_file = args.i
output_folder = args.o
log_output_file = args.l

#open and loop through files in input directory
retrieve_gca = {}
with open(input_file, 'r') as file_in:
	for line in file_in:
		line = line.strip()
		if line:  # skip empty lines
			retrieve_gca[line] = {}

#retrieve taxonomy for genomes of interest
taxonomy = {}
with open("../pp_db/GCA2taxonomy.txt", 'r') as file_in:
	next(file_in) #skip first line 
	for line in file_in:
		line = line.strip()
		info = line.split("\t")
		#gcas
		gcas = info[2].split(";")
		#taxonomy
		taxa = info[1].split(";")
		domain = taxa[0][3:].lower()
		phylum = taxa[1][3:].replace(" ", "_")
		family = taxa[4][3:].replace(" ", "_")
		genus = taxa[5][3:].replace(" ", "_")
		#save information for each gca
		for gca in gcas:
			taxonomy[gca] = {}
			taxonomy[gca]["domain"] = domain
			taxonomy[gca]["file_name"] = phylum + "_" + family + "_" + genus

#open log file
file_log = open(log_output_file, 'w')

#download genomes
for gca in retrieve_gca:
	if gca not in taxonomy:
		gca = gca.replace("GCF", "GCA")
	if gca in taxonomy:
		download_command = f"ncbi-genome-download {taxonomy[gca]['domain']} --section genbank --assembly-accession {gca} --format fasta -o {output_folder}"
		process = subprocess.Popen(download_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		output, error = process.communicate()
		return_code = process.returncode
		if return_code == 0:
			file_log.write(f"Genome downloaded successfully: {gca}\n")
			file_log.flush()
		else:
			file_log.write(f"Error in downloading genome: {gca}\n")
			file_log.flush()
	
		#rename files
		rename_command = f"cp {output_folder}/genbank/{taxonomy[gca]['domain']}/{gca}/*.fna.gz {output_folder}/{taxonomy[gca]['file_name']}_{gca}.fna.gz"
		process = subprocess.Popen(rename_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		output, error = process.communicate()
		return_code = process.returncode
		if return_code == 0:
			file_log.write(f"Genome renamed successfully: {gca}\n")
			file_log.flush()
		else:
			file_log.write(f"Error in renaming genome: {gca}\n")
			file_log.flush()
	else:
		file_log.write(f"Error because gca not found: {gca}\n")

print(f'{count_all} genome accessions.')
print(f'{count_unique} unique genome accessions.')