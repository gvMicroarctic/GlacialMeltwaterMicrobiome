#!/usr/bin/env python3

#Associate taxid to taxonomy from NCBI datasets

#load taxid names
taxid_names = {}
with open("../pp_db/names.dmp", 'r') as in_file:
	for line in in_file:
		line = line.strip()
		info = line.split("\t")
		if info[6] == "scientific name":
			taxid_names[info[0]] = info[2]
			

#load taxon ranks and taxid relations
taxid_ranks = {}
taxid_rel = {}
with open("../pp_db/nodes.dmp", 'r') as in_file:
	for line in in_file:
		line = line.strip()
		info = line.split("\t")
		taxid_ranks[info[0]] = info[4] #child = rank
		taxid_rel[info[0]] = info[2] #child = parent
		
#combine information
data = {}
for child in taxid_rel:
	data[child] = {}
    
for child in taxid_rel:
	data[child][taxid_ranks[child]] = taxid_names[child] #taxid - rank = name
	parent = taxid_rel[child]
	while parent != "1" and parent != "131567":
		if (taxid_ranks[parent] != "clade"):
			data[child][taxid_ranks[parent]] = taxid_names[parent] #taxid - rank = name
			parent = taxid_rel[parent]
			#print("four")
		elif ("clade" not in data[child]) or (taxid_ranks[parent] == "clade" and "group" in taxid_names[parent]):
			data[child][taxid_ranks[parent]] = taxid_names[parent] #taxid - rank = name
			parent = taxid_rel[parent]
			#print("three")
		else:
			parent = taxid_rel[parent]

#print
ranks = ["domain", "phylum", "class", "order", "family","genus", "species", "clade"]
with open("../pp_db/taxid2taxonomy.txt", 'w') as out_file:
	out_file.write(f'taxid\ttaxonomy')
	for taxid in data:
		out_file.write(f'\n{taxid}\t')
		for rank in ranks:
			rank_first = rank[0]
			if rank == "superkingdom":
				rank_first = "d"	
			if rank == "clade":
				rank_first = "j"	
			if rank in data[taxid]:
				out_file.write(f'{rank_first}__{data[taxid][rank]};')
			else:
				out_file.write(f'{rank_first}__Unclassified;')
				
#add merged entries to taxid2taxonomy.txt
taxonomy = {}
with open("../pp_db/taxid2taxonomy.txt", 'r') as in_file:
	for line in in_file:
		line = line.strip()
		info = line.split("\t")
		taxonomy[info[0]] = info[1]

with open("../pp_db/taxid2taxonomy.txt", 'a') as out_file:
	with open("../pp_db/merged.dmp", 'r') as in_file:
		for line in in_file:
			line = line.strip()
			info = line.split("\t")
			
			out_file.write(f'\n{info[0]}\t{taxonomy[info[2]]}')
		
			
	
		
		
		
		
		
				
				
