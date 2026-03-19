#!/usr/bin/env python3

import argparse

#MGIs to genes

#get variables from command line
parser = argparse.ArgumentParser(description='Screen bed file to look for MGIs.')
parser.add_argument('-mgi', type=str, help='Input MGI positions.')
parser.add_argument('-gff', type=str, help='Input gff annotation.')
parser.add_argument('-o', type=str, help='Output MGI positions and genes.')
args = parser.parse_args()

#name variables from command line
input_file_mgi = args.mgi
input_file_gff = args.gff
output_file = args.o

#Get CDS annotation information from gff file
genes = {}
with open(input_file_gff, 'r') as file_in:
	for line in file_in:
	
		if line.startswith(">"):
			break  #exit the loop immediately if >

		if not line.startswith("#"): 
			line = line.strip()
			info = line.split("\t")
			
			#format contig name
			contig0 = info[0].split("_")
			contig = contig0[2] + "_" + contig0[3] + "_" + contig0[4] + "_" + contig0[5]
			
			#format gene name
			gene0 = info[8].split(';')[0]
			gene = gene0.replace('ID=', '')

			#in our study, considered match if overlap with mgi min and max (gene does not have to be necessarly completely inside)
			gene_min = float(info[3])
			gene_max = float(info[4]) 

			if contig not in genes:
				genes[contig] = {}
			
			#save gene start and end position
			genes[contig][gene] = {}
			genes[contig][gene] = (gene_min, gene_max) #gene info = start and end

#Open output file
file_out = open(output_file, 'w')

#Look through all MGIs
with open(input_file_mgi, 'r') as file_in:
	for line in file_in:
		line = line.strip()
		info = line.split("\t")
		contig = info[0]
		
		#start and end positions
		pos = info[2]
		mgi_min = float(pos.split("-")[0])
		mgi_max = float(pos.split("-")[1])
		
		#check if genes within the mgi area
		if contig in genes:
			overlaps = { gene:mgi_min < gene_max and gene_min < mgi_max
				for gene,(gene_min,gene_max) in genes[contig].items() }
			
			for gene, does_overlap in overlaps.items():
				if does_overlap:
					file_out.write(f'{line}\t{gene}\t{int(genes[contig][gene][0])}-{int(genes[contig][gene][1])}\n')