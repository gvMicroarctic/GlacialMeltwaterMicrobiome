#!/usr/bin/env python3

#Get viral taxonomy from blastp output (m8)

from taxonomy_assignment_functions import common_taxonomy

#save information from m8 file
db_protein = {}
contig2protein = {}
with open("./viral_contig/blastp/blastp_refseq_viral_taxonomy_filtered.m8", 'r') as in_file:
	for line in in_file:
		line = line.strip()
		info = line.split("\t")

		db_protein[info[1]] = {} #db_protein = {}

		#contig
		contig_all = info[0].split('_')
		contig = contig_all[0] + "_" + contig_all[1]
		#print(f"{contig}")
		if contig not in contig2protein:
			contig2protein[contig] = {}
		contig2protein[contig][info[1]] = {} #contig - db_protein = {}

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
with open("../viruses_refseq/prot.accession2taxid", 'r') as in_file:
	for line in in_file:
		line = line.strip()
		info = line.split("\t")

		if info[1] in db_protein:
		
			#print(f"one {info[2]}")
			
			#print(f"three {taxonomy[info[2]]}")

			tax_info = taxonomy[info[2]].split(";")
			
			#print(f"two {info[1]} and {info[2]} and {tax_info}")

			#get taxonomic information for all db_proteins
			taxonomy_all[info[1]] = tax_info[0] + ";" + tax_info[1] + ";" + tax_info[2] + ";" + tax_info[3] + ";" + tax_info[4]

#output file
out_file = open("./viral_contig/blastp/taxonomy_viruses_blastp.txt", "w")

#merge information and print classification
for contig in contig2protein:
	
	#create list containing all taxonomic paths
	paths = []
	for protein in contig2protein[contig]:
		paths.append(taxonomy_all[protein])
	
	#get taxonomic classification
	classification = common_taxonomy(paths)
	out_file.write(f'{contig}\t{classification}\n')
