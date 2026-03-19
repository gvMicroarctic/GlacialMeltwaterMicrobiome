#!/usr/bin/env python3

#Merge viral taxonomy

#define functions

#create dictionary with taxonomic paths (no Unclassified)
def create_screen_path(c, dict_c, rank):

	dict_p = {}
	
	for tool in dict_c[c]:
		
		path = dict_c[c][tool]
		
		#trim from rank
		if rank == "all":
			path_new = path
		else:
			path_new = path.split(rank)[0]
		
		#remove unclassifies	
		if "d__Unclassified" in path_new:
			path_new = path_new.replace("d__Unclassified;","")
				
		if "p__Unclassified" in path_new:
			path_new = path_new.replace("p__Unclassified;","")
				
		if "c__Unclassified" in path_new:
			path_new = path_new.replace("c__Unclassified;","")
				
		if "o__Unclassified" in path_new:
			path_new = path_new.replace("o__Unclassified;","")
				
		if "f__Unclassified" in path_new:
			path_new = path_new.replace("f__Unclassified","")
				
		dict_p[path_new] = {}
		
	#return dictionary with taxonomic paths without unclassified
	return dict_p
	
#create dictionary with taxonomic paths (no Unclassified)
def path_in(path):

	if "d__" not in path:
		path = path + ";d__Unclassified"
					
	if "p__" not in path:
		path = path + ";p__Unclassified"
					
	if "c__" not in path:
		path = path + ";c__Unclassified"
					
	if "o__" not in path:
		path = path + ";o__Unclassified"
					
	if "f__" not in path:
		path = path + ";f__Unclassified"
				
	path_info = path.split(";")
				
	for taxon in path_info:
		if "d__" in taxon:
			domain_t = taxon
	for taxon in path_info:
		if "p__" in taxon:
			phylum_t = taxon
	for taxon in path_info:
		if "c__" in taxon:
			class_t = taxon
	for taxon in path_info:
		if "o__" in taxon:
			order_t = taxon
	for taxon in path_info:
		if "f__" in taxon:
			family_t = taxon

	#return string to print
	return domain_t, phylum_t, class_t, order_t, family_t
	
#load final vOTUs
final_contigs = {}
with open("./viral_contig/viral_contigs_unique_final_5000_clean.fasta", 'r') as in_file:
	for line in in_file:
		line = line.strip()
		if line.startswith(">"):
			contig = line[1:]
			final_contigs[contig] = {}

contigs = {}
with open("./viral_contig/blastp/taxonomy_viruses_blastp.txt", 'r') as in_file:
	for line in in_file:
		line = line.strip()
		info = line.split("\t")
		if info[0] in final_contigs:
			if info[0] not in contigs:
				contigs[info[0]] = {}
				contigs[info[0]]["blastp"] = info[1]
			else:
				contigs[info[0]]["blastp"] = info[1]
			
with open("./viral_contig/blastn/taxonomy_viruses_blastn.txt", 'r') as in_file:
	for line in in_file:
		line = line.strip()
		info = line.split("\t")
		if info[0] in final_contigs:
			if info[0] not in contigs:
				contigs[info[0]] = {}
				contigs[info[0]]["blastn"] = info[1]
			else:
				contigs[info[0]]["blastn"] = info[1]
		
with open("./viral_contig/mmseqs2/taxonomy_viruses_mmseq.txt", 'r') as in_file:
	for line in in_file:
		line = line.strip()
		info = line.split("\t")
		if info[0] in final_contigs:
			if info[0] not in contigs:
				contigs[info[0]] = {}
				contigs[info[0]]["mmseq"] = info[1]
			else:
				contigs[info[0]]["mmseq"] = info[1]
			
with open("./viral_contig/genomad/taxonomy_viruses_genomad.txt", 'r') as in_file:
	for line in in_file:
		line = line.strip()
		info = line.split("\t")
		if info[0] in final_contigs:
			if info[0] not in contigs:
				contigs[info[0]] = {}
				contigs[info[0]]["genomad"] = info[1]
			else:
				contigs[info[0]]["genomad"] = info[1]
		
#open output file
out_file = open("./viral_contig/votu_taxonomy.txt", 'w')

#set up dictionary to save contigs for further steps
contigs_new = {}

#initiate dictionary to save contig names
contigs_name = {}

#repeat on contigs that have been not assigned yet, removing from lower taxa
for rank in ("all", "f__", "o__", "c__", "p__", "d__"):
	
	contigs_tax = contigs_name.copy()
	contigs_name = {}
	
	#for each contig
	for contig in contigs:
		
		#only if in contig left to retrieve
		if (contig in contigs_tax) or (rank == "all"):
		
			#create dictionary containing taxonomic paths without Unclassified
			screen_path =  create_screen_path(contig, contigs, rank)
			#get taxonomic path that is the longest
			path_max = max(screen_path, key = len)
			#set the status
			status = "out"
		
			#dictionary with keys to remove
			path_rem = {}
		
			while status == "out": 
				
				#check if all shorter taxonomic paths are within the longest one
				status = "in"
				for path_t in screen_path:
					if path_t not in path_max:
						#set the status
						status = "out"
	
				#get final annotation
				if status == "in":
			
					#get taxonomy associated to final path to print
					domain_t, phylum_t, class_t, order_t, family_t = path_in(path_max)
					#print final path
					
					#domain can be Unclassified even if the rest is classified because of how blastn and blastp work
					if domain_t != "d__Unclassified":
						out_file.write(f'{contig}\t{domain_t}\t{phylum_t}\t{class_t}\t{order_t}\t{family_t}\n')
					else:
						out_file.write(f'{contig}\td__Viruses\t{phylum_t}\t{class_t}\t{order_t}\t{family_t}\n')

				#stay in while loop and repeat		
				elif status == "out":
				
					#feed dictionary with keys to remove
					path_rem[path_max] = {}
					#remove keys
					screen_path_new = screen_path.copy()
					for rem in path_rem:
						if rem in screen_path_new:
							del screen_path_new[rem]
			
					#get new path_max
					if len(screen_path_new) != 0:
						path_max = max(screen_path_new, key = len)
				
					else:
			
						#set the status
						status = "in"
						#save contigs in dictionary so that carried to further steps (remove families)
						contigs_name[contig] = {}

#print "Unclassified" also for contigs that are not present in any of the taxonomic classifications
for contig in final_contigs:
	if contig not in contigs:
		out_file.write(f'{contig}\td__Viruses\tp__Unclassified\tc__Unclassified\to__Unclassified\tf__Unclassified\n')
