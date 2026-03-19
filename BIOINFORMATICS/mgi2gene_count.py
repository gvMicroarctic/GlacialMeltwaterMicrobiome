#!/usr/bin/env python3

#Get number of MGIs per MAG

#open file obtained from mgi2gene.py
mag_count = {}
#N.B: I am retrieving both internal and external MGIs
with open('./MAG/MGI/mgi_positions.txt', 'r') as file_in: 

	for line in file_in:
		line = line.strip()
		info = line.split("\t")
		
		#mag
		mag0 = info[0].split("_")
		mag = mag0[2] + "_" + mag0[3]
		
		#create dictionary
		if mag not in mag_count:
			mag_count[mag] = {}
			mag_count[mag]["in"] = 0
			mag_count[mag]["ex"] = 0
			mag_count[mag]["tot"] = 0
			
		if info[1] == "internal_mgi":
			mag_count[mag]["in"] += 1
			mag_count[mag]["tot"] += 1
			
		elif info[1] == "external_mgi":
			mag_count[mag]["ex"] += 1
			mag_count[mag]["tot"] += 1
	
#open output file to print information
file_out = open("./MAG/MGI/mgi2gene_count.csv", 'w') # just this name to change

#print information
file_out.write(f'MAG,internal_mgi,external_mgi,total_mgi\n')
for mag in mag_count:
	file_out.write(f'{mag},{mag_count[mag]["in"]},{mag_count[mag]["ex"]},{mag_count[mag]["tot"]}\n')
