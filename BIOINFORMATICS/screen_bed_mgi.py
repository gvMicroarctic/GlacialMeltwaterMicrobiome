#!/usr/bin/env python3

import argparse

#Screen bed file to look for MGIs
#Bed file: contig_mag - start position (not included) - end position (included) - coverage: 

#get variables from command line
parser = argparse.ArgumentParser(description='Screen bed file to look for MGIs.')
parser.add_argument('-i', type=str, help='Input bed file.')
parser.add_argument('-con', type=str, help='Output contig coverage file.')
parser.add_argument('-mag', type=str, help='Output MAG coverage file.')
parser.add_argument('-mgi', type=str, help='Output MGI positions.')
args = parser.parse_args()

#name variables from command line
input_file = args.i
output_file_con = args.con
output_file_mag = args.mag
output_file_mgi = args.mgi

#Get information to calculate mean contig average
contig_all = {}
mag_all = {}
with open(input_file, 'r') as file_in:
	for line in file_in:
		line = line.strip()
		info = line.split("\t")
		contig = info[0]
		start = int(info[1]) + 1
		end = int(info[2]) + 1 #+1 so then last number is also included
		coverage = float(info[3])
		mag0 = contig.split("_")
		mag = "bin_" + mag0[3]
		
		if contig not in contig_all:
			contig_all[contig] = {}
			contig_all[contig]["total_cov"] = 0
			contig_all[contig]["total_pos"] = 0
			
		if mag not in mag_all:
			mag_all[mag] = {}
			mag_all[mag]["total_cov"] = 0
			mag_all[mag]["total_pos"] = 0
		
		for pos in range(start, end):
			contig_all[contig]["total_cov"] += coverage 
			contig_all[contig]["total_pos"] += 1 
			
			mag_all[mag]["total_cov"] += coverage 
			mag_all[mag]["total_pos"] += 1 

#Output files
file_out_con = open(output_file_con, 'w')
file_out_mag = open(output_file_mag, 'w')
file_out_mgi = open(output_file_mgi, 'w')

#Calculate and print mean contig average
for contig in contig_all:
	mean = contig_all[contig]["total_cov"]/contig_all[contig]["total_pos"]
	file_out_con.write(f'{contig}\t{mean}\t{contig_all[contig]["total_pos"]}\n')
	
#Calculate and print coverage for MAGs; save coverage as needed to mine for MGIs
coverage_mean = {}
for mag in mag_all:
	mean = mag_all[mag]["total_cov"]/mag_all[mag]["total_pos"]
	file_out_mag.write(f'{mag}\t{mean}\t{mag_all[mag]["total_pos"]}\n')
	coverage_mean[mag] = mean
	
#Bellas et al. (2020) outlined two methods to identify MGIs
#1. if mean coverage >= 5x: MGI if coverage dropped to > 25% of the mean over >= 100 bp region
#2. if mean coverage 2-5x: MGI if coverage dropped over >= 200 bp region
#I am calculating means at MAG levels as 
	
#Save all positions that have a coverage < the 25% of the contig average
coverage_low = {} #save positions with coverage lower than threshold
with open(input_file, 'r') as file_in:
	for line in file_in:
		line = line.strip()
		info = line.split("\t")
		contig = info[0]
		mag0 = contig.split("_")
		mag = "bin_" + mag0[3]
		
		start = int(info[1]) + 1
		end = int(info[2]) + 1 #+1 so then last number is also included
		coverage = float(info[3])
		
		#if mag mean coverage is >= 5
		if coverage_mean[mag] >= 5:
			
			#create threshold coverage value
			thr_cov = coverage_mean[mag]/4
			
			#save all stretches of low coverage bases (regardless from length)
			if contig not in coverage_low:
				coverage_low[contig] = {}
				stretch = 1
				coverage_low[contig][stretch] = {}
		
			for pos in range(start, end):
				if coverage < thr_cov:
				
					if stretch not in coverage_low[contig]:
						coverage_low[contig][stretch] = {}
				
					coverage_low[contig][stretch][pos] = {}
				
				else:
					stretch += 1
					
		#if mag mean coverage between 2 and 5
		elif coverage_mean[mag] > 2 and coverage_mean[mag] < 5:
			
			#save all stretches of low coverage bases (regardless from length)
			if contig not in coverage_low:
				coverage_low[contig] = {}
				stretch = 1
				coverage_low[contig][stretch] = {}
		
			for pos in range(start, end):
				if coverage == 0:
				
					if stretch not in coverage_low[contig]:
						coverage_low[contig][stretch] = {}
				
					coverage_low[contig][stretch][pos] = {}
				
				else:
					stretch += 1
				
#Save and print all MGIs
for contig in coverage_low:

	mag0 = contig.split("_")
	mag = "bin_" + mag0[3]
	
	#get threshold on lengths
	if coverage_mean[mag] >= 5:
		thr_len = 100
	elif coverage_mean[mag] > 2 and coverage_mean[mag] < 5:
		thr_len = 200
	
	#select and print MGIs
	for stretch in coverage_low[contig]:
		if len(coverage_low[contig][stretch]) >= thr_len:
			file_out_mgi.write(f'{contig}')
			count_in = 0
			for pos in coverage_low[contig][stretch]:
				count_in += 1
				if count_in == 1:
					pos_start = pos
				#file_out_mgi.write(f'{pos},')
			
			#print(pos)
			if (pos == contig_all[contig]["total_pos"]) or (pos_start == 1):
				#N.B. pos_start for bed_format should be "pos_start -1"; to keep in mind
				file_out_mgi.write(f'\texternal_mgi\t{pos_start}-{pos}\n')
			else:
				file_out_mgi.write(f'\tinternal_mgi\t{pos_start}-{pos}\n')
				
