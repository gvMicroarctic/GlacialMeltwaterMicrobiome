#!/usr/bin/env python3

#Merge all MAG files into one

import os
import re

#specify the directory path
folder_path = './MAG/medium_high_drep_95/dereplicated_genomes/'

#list all files in the directory
files = os.listdir(folder_path)

#open output file
file_out = open("./MAG/all_mags.fasta", "w") 

#iterate through each file
for file_name in files:
	
	#create the full file path by joining the folder path and the file name
	file_path = os.path.join(folder_path, file_name)
	
	#retrieve only MAG name
	if file_name[0].isalpha(): #if it starts with letters
		name = re.search(r'(\w+).(\d+)', file_name)[0]
	else:
		name = re.search(r'(\d+)', file_name)[0]
	
	#open file
	with open(file_path, 'r') as in_file:
		for line in in_file:
			line = line.strip()
			if line.startswith('>'):
				if "length" in line:
					header = line.split("_length")[0]
				else:
					header0 = line.split("flag")[0]
					header = header0.strip()
				contig_name = header
				file_out.write(f'{contig_name}\n')
			else:
				file_out.write(f'{line}\n')
