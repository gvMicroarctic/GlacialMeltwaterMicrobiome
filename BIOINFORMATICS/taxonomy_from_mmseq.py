#!/usr/bin/env python3

#Taxonomy from blast output (m8) obtained with mmseq

from taxonomy_assignment_functions import common_taxonomy

#save information from m8 file
db_group = {}
contig2group = {}
with open("./viral_contig/mmseqs2/mmseq_vog_taxonomy_filtered.m8", 'r') as in_file:
	for line in in_file:
		line = line.strip()
		info = line.split("\t")
		
		#vog viral group
		group = info[1].split('.')[0]
		db_group[group] = {} #db_group = {}
		
		#contig
		contig = info[0]
		if contig not in contig2group:
			contig2group[contig] = {}
		contig2group[contig][group] = {} #contig - db_group = {}
		
#get taxonomic information
taxonomy = {}
with open("../taxonomy/taxid2taxonomy.txt", 'r') as in_file:
	lines = in_file.readlines()
	for line in lines[1:]:
		line = line.strip()
		info = line.split("\t")
		taxonomy[info[0]] = info[1] #taxid = taxonomy

#open file with correspondances between accession numbers and taxids
taxonomy_all = {}
with open("../mmseqs2_db/vog.lca.tsv", 'r') as in_file:
	for line in in_file:
		line = line.strip()
		info = line.split("\t")
		
		if info[0] in db_group:
		
			tax_info = taxonomy[info[4]].split(";")
			
			#get taxonomic information for all db_proteins
			taxonomy_all[info[0]] = tax_info[0] + ";" + tax_info[1] + ";" + tax_info[2] + ";" + tax_info[3] + ";" + tax_info[4]

#output file
out_file = open("./viral_contig/mmseqs2/taxonomy_viruses_mmseq.txt", "w")

#merge information and print classification
for contig in contig2group:
	
	#create list containing all taxonomic paths
	paths = []
	for group in contig2group[contig]:
		paths.append(taxonomy_all[group])
	
	#get taxonomic classification
	classification = common_taxonomy(paths)
	out_file.write(f'{contig}\t{classification}\n')
	
