#!/usr/bin/env python3

#Taxonomy from blastn output (m8)

from taxonomy_assignment_functions import common_taxonomy

#save information from m8 file
db_seq = {}
contig2seq = {}
with open("./viral_contig/blastn/blastn_IMGVR_viral_taxonomy.m8", 'r') as in_file:
	for line in in_file:
		line = line.strip()
		info = line.split("\t")
		
		#IMG_VR sequence
		seq = info[1].split('|')[0]
		db_seq[seq] = {} #db_seq = {}
		
		#contig
		contig_all = info[0].split('_')
		contig = contig_all[0] + "_" + contig_all[1]

		if contig not in contig2seq:
			contig2seq[contig] = {}
		contig2seq[contig][seq] = {} #contig - db_seq = {}
		
#get taxonomic information
taxonomy_all = {}
with open("../IMG_VR_2022-12-19_7/IMGVR_all_Sequence_information.tsv", 'r') as in_file:
	lines = in_file.readlines()
	for line in lines[1:]:
		line = line.strip()
		info = line.split("\t")
		if info[0] in db_seq:
			
			tax_info = info[14].split(";")
			
			taxonomy_all[info[0]] = "d__Viruses" #protein = taxonomic path
			
			if "p__" in info[14]:
				taxonomy_all[info[0]] = taxonomy_all[info[0]] + ";" + tax_info[2]
			else:
				taxonomy_all[info[0]] = taxonomy_all[info[0]] + ";" + "p__Unclassified"
				
			if "c__" in info[14]:
				taxonomy_all[info[0]] = taxonomy_all[info[0]] + ";" + tax_info[3]
			else:
				taxonomy_all[info[0]] = taxonomy_all[info[0]] + ";" + "c__Unclassified"
			
			if "o__" in info[14]:
				taxonomy_all[info[0]] = taxonomy_all[info[0]] + ";" + tax_info[4]
			else:
				taxonomy_all[info[0]] = taxonomy_all[info[0]] + ";" + "o__Unclassified"	
				
			if "f__" in info[14]:
				taxonomy_all[info[0]] = taxonomy_all[info[0]] + ";" + tax_info[5]
			else:
				taxonomy_all[info[0]] = taxonomy_all[info[0]] + ";" + "f__Unclassified"	

#output file
out_file = open("./viral_contig/blastn/taxonomy_viruses_blastn.txt", "w")

#merge information and print classification
for contig in contig2seq:
	
	#create list containing all taxonomic paths
	paths = []
	for seq in contig2seq[contig]:
		paths.append(taxonomy_all[seq])
	
	#get taxonomic classification
	classification = common_taxonomy(paths)
	out_file.write(f'{contig}\t{classification}\n')
