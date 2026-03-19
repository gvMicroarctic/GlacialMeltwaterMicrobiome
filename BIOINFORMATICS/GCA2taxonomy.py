#!/usr/bin/env python3

#associate taxonomy to GCA genomes

#load taxid names
gca = {}
with open("../pp_db/assembly_summary_genbank.txt", 'r') as file:
	lines = file.readlines()
	for line in lines[2:]: #do not read the first two lines
		line = line.strip()
		info = line.split("\t")
		gca[info[6]] = {}
with open("../pp_db/assembly_summary_genbank.txt", 'r') as file:
	lines = file.readlines()
	for line in lines[2:]: #do not read the first two lines
		line = line.strip()
		info = line.split("\t")
		gca[info[6]][info[0]] = {}
			
#load taxonomy
taxonomy = {}
with open("../pp_db/taxid2taxonomy.txt", 'r') as file:
	for line in file:
		line = line.strip()
		info = line.split("\t")
		taxonomy[info[0]] = info[1]
		
#combine information
with open("../pp_db/GCA2taxonomy.txt", 'w') as file:
	file.write(f'taxid\ttaxonomy\tGCAs')
	for taxid in gca:
		if taxid in taxonomy: #a few GCA taxids are in taxonomy dataset, exclude them
			file.write(f'\n{taxid}\t{taxonomy[taxid]}\t')
			for accession in gca[taxid]:
				file.write(f'{accession};')
