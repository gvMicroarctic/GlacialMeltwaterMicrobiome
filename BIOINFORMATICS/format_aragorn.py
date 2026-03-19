#!/usr/bin/env python3

#Format aragorn output and create fasta file

import argparse

#get variables from command line
parser = argparse.ArgumentParser(description='Format aragorn output and create fasta file.')
parser.add_argument('-i', type=str, help='Input aragorn file.')
parser.add_argument('-o', type=str, help='Output aragorn fasta file.')
parser.add_argument('-c1', type=str, help='Start of contig name 1.')
parser.add_argument('-c2', type=str, help='Start of contig name 2.')
args = parser.parse_args()

#name variables from command line
input_file = args.i
output_file = args.o
contig_name1 = args.c1
contig_name2 = args.c2

#output file
file_out = open(output_file, "w")

#open aragorn file
contigs = {}
seq = 0
with open(input_file, 'r') as in_file:
	for line in in_file:
		line = line.strip()
		info = line.split("\t")
		if line.startswith(contig_name1) or line.startswith(contig_name2):
			contig = line
		elif line.startswith(">"):
			seq = 1
			header =  line.replace(" ", "_")
			file_out.write(f'{header}_{contig}\n')
		elif line == "":
			seq = 0
		else:
			if seq == 1:
				sequence = line.upper()
				file_out.write(f'{sequence}\n')
