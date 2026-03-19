#!/usr/bin/env python3

#import os
#from Bio import SeqIO
#import pandas as pd

#format output from CAT taxonomy assigner

#open CAT taxonomy file (only_official names)
input_file = "./assembly/cat_taxonomy_only_official.txt"
output_file = "./assembly/cat_taxonomy_only_official_formatted.csv"

with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:

	#add header to output file
	out_file.write("contig,domain,phylum,class,order,family,genus,species\n")

	for line in in_file:
		line = line.strip()
		info = line.split("\t")
		if not line.startswith("#"):
			line_clean0 = line.replace("*", "")
			line_clean1 = line_clean0.replace("no support", "Unclassified")
			line_clean = line_clean1.replace("NA", "Unclassified")
			info = line_clean.split("\t")
			if len(info) < 12: #contig not taxonomically assigned: "ORF has no hit to database" or "no taxid found"
				contig = info[0]
				domain_t = "Unclassified"
				phylum_t = "Unclassified"
				class_t = "Unclassified"
				order_t = "Unclassified"
				family_t = "Unclassified"
				genus_t = "Unclassified"
				species_t = "Unclassified"
			else:
				contig = info[0]
				domain_t = info[5]
				phylum_t = info[6]
				class_t = info[7]
				order_t = info[8]
				family_t = info[9]
				genus_t = info[10]
				species_t = info[11]
			
			#use a list to hold the fields and then join them with tabs
			output_data = [contig, domain_t, phylum_t, class_t, order_t, family_t, genus_t, species_t]
			
			#join the list elements with a tab and add a newline character
			output_line = ",".join(output_data) + "\n"
			
			#write the processed line to the output file
			out_file.write(output_line)




