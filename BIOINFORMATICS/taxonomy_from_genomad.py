#!/usr/bin/env python3

#Taxonomy from genomad output

#open output file
file_out = open("./viral_contig/genomad/taxonomy_viruses_genomad.txt", 'w')

#get and format taxonomic information
with open("./viral_contig/genomad/viral_contigs_unique_final_5000_clean_summary/viral_contigs_unique_final_5000_clean_virus_summary.tsv", 'r') as file_in:
	lines = file_in.readlines()
	for line in lines[1:]:
		line = line.strip()
		info = line.split("\t")
		
		#format contig name
		contig0 = info[0]
		if "|" in contig0:
			info_c = contig0.split("|")
			contig = info_c[0]
		else:
			contig = contig0
			
		taxonomy = info[10]
		
		if taxonomy == "Unclassified":
			file_out.write(f'{contig}\td__Viruses;p__Unclassified;c__Unclassified;o__Unclassified;f__Unclassified\n')
		else:
			taxa = taxonomy.split(";")
			
			domain_r = "d__" + taxa[0]

			if taxa[3] == "":
				phylum_r = "p__Unclassified"
			else:
				phylum_r = "p__" + taxa[3]
				
			if taxa[4] == "":
				class_r = "c__Unclassified"
			else:
				class_r = "c__" + taxa[4]
			
			if taxa[5] == "":
				order_r = "o__Unclassified"
			else:
				order_r = "o__" + taxa[5]
				
			if taxa[6] == "":
				family_r = "f__Unclassified"
			else:
				family_r = "f__" + taxa[6]
			
			file_out.write(f'{contig}\t{domain_r};{phylum_r};{class_r};{order_r};{family_r}\n')
