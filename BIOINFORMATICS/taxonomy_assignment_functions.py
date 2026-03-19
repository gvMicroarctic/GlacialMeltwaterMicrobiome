#!/usr/bin/env python3

#Functions to create taxonomic assignments from blastp, blastn, mmseq and genomad

def classify_taxon(classification, taxon, combined, count_total, taxon_path, trimmed_path, count_f, rank):
	proportion = combined[rank][taxon] / count_total

	if proportion > 0.80:  # Save if a taxon is present with more than 80% occurrence
		classification = taxon_path[taxon]

	else:  #check agreement at higher taxonomic levels
	
		combined_count = 0
		for new_path in trimmed_path:
			if new_path in taxon_path[taxon]: ###mistake
				combined_count += trimmed_path[new_path]

		proportion = combined_count / count_total

		if proportion > 0.8:  # Save if > 80%
			classification = taxon_path[taxon]

		if proportion > 0.5:  # If more than one family has high abundance
			count_f += 1

	return count_f, classification

def common_taxonomy(paths):

	combined = {"domain": {}, "phylum": {}, "class": {}, "order": {}, "family": {}}
	taxon_path = {}
	count_total = 0
	trimmed_path = {}

	for path in paths:

		#split taxonomic path into different taxa
		tax_info = path.split(";")
		taxon_d = tax_info[0]
		taxon_p = tax_info[1]
		taxon_c = tax_info[2]
		taxon_o = tax_info[3]
		taxon_f = tax_info[4]

		#create dictionary with counts associated to each taxon

		#initialise
		if taxon_d not in combined["domain"]:
			combined["domain"][taxon_d] = 0
		if taxon_p not in combined["phylum"]:
			combined["phylum"][taxon_p] = 0
		if taxon_c not in combined["class"]:
			combined["class"][taxon_c] = 0
		if taxon_o not in combined["order"]:
			combined["order"][taxon_o] = 0
		if taxon_f not in combined["family"]:
			combined["family"][taxon_f] = 0
		
		#count
		combined["domain"][taxon_d] += 1
		combined["phylum"][taxon_p] += 1
		combined["class"][taxon_c] += 1
		combined["order"][taxon_o] += 1
		combined["family"][taxon_f] += 1

		#total
		count_total += 1

		#create taxon associated with different taxonomic paths with unclassied lower levels
		taxon_path[taxon_p] = taxon_d + ";" + taxon_p + ";c__Unclassified;o__Unclassified;f__Unclassified"
		taxon_path[taxon_c] = taxon_d + ";" + taxon_p + ";" + taxon_c + ";o__Unclassified;f__Unclassified"
		taxon_path[taxon_o] = taxon_d + ";" + taxon_p + ";" + taxon_c + ";" + taxon_o + ";f__Unclassified"
		taxon_path[taxon_f] = taxon_d + ";" + taxon_p + ";" + taxon_c + ";" + taxon_o + ";" + taxon_f
		
		#remove tailing unclassified taxa and save in a new list (new_paths)
		new_path = path.replace("d__Unclassified;p__Unclassified;c__Unclassified;o__Unclassified;f__Unclassified", "")
		new_path = new_path.replace(";p__Unclassified;c__Unclassified;o__Unclassified;f__Unclassified", "")
		new_path = new_path.replace(";c__Unclassified;o__Unclassified;f__Unclassified", "")
		new_path = new_path.replace(";o__Unclassified;f__Unclassified", "")
		new_path = new_path.replace(";f__Unclassified", "")

		#save the paths without unclassified paths in dictionary: new path = counts
		if new_path not in trimmed_path:
			trimmed_path[new_path] = 0 #initialise
		trimmed_path[new_path] += 1 

	classification = "NA"

	#save taxonomy if already agreement at the family level
	count_f = 0
	rank = "family"
	for taxon in combined[rank]:
		if taxon != "f__Unclassified":
			count_f, classification = classify_taxon(classification, taxon, combined, count_total, taxon_path, trimmed_path, count_f, rank)
	if count_f > 1:
		classification = "NA"

	#save taxonomy if agreement at order, class or phylum level
	for rank in ("order", "class", "phylum"):
		count_f = 0
		if (classification == "NA"):
			for taxon in combined[rank]:
				if rank == "order":
					if taxon != "o__Unclassified":
						count_f, classification = classify_taxon(classification, taxon, combined, count_total, taxon_path, trimmed_path, count_f, rank)
				if rank == "class":
					if taxon != "c__Unclassified":
						count_f, classification = classify_taxon(classification, taxon, combined, count_total, taxon_path, trimmed_path, count_f, rank)
				if rank == "phylum":
					count_f, classification = classify_taxon(classification, taxon, combined, count_total, taxon_path, trimmed_path, count_f, rank)
			if count_f > 1:
				classification = "NA"

	#save taxonomy if not agreement at previous taxonomic levels
	if (classification == "NA"):
		classification = "d__Viruses;p__Unclassified;c__Unclassified;o__Unclassified;f__Unclassified"

	return classification