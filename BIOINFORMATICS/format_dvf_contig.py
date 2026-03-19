#!/usr/bin/env python3

#open dvf file and retrieve only contigs with scores >= 0.9 and pvalue <= 0.05

#define functions to work with fasta files
def read_fasta(file_path):
    """Reads a FASTA file and returns a dictionary of sequences."""
    sequences = {}
    current_sequence = ""
    
    with open(file_path, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                current_sequence0 = line[1:]
                current_sequence = current_sequence0.split(" ")[0]
                sequences[current_sequence] = ""
            else:
                sequences[current_sequence] += line
    
    return sequences

def select_sequence(fasta_file_path, selected_sequence_names):
    """Selects sequences from a FASTA file based on a list of sequence names."""
    all_sequences = read_fasta(fasta_file_path)
    
    selected_sequences = {name: sequence for name, sequence in all_sequences.items() if name in selected_sequence_names}
    
    return selected_sequences
    
#open dvf score files
dvf_file = open("./viral_contig/dvf/final.contigs_5000.fa_gt1bp_dvfpred.txt", "r")
dvf_file.readline() #discard header

#define empty list to save contig names
contig_names = []

for line in dvf_file: #iterate through the file line by line
	columns = line.split("\t") #split in different columns
	contig0 = columns[0]
	contig = contig0.split(" ")[0]
	score = float(columns[2])
	pvalue = float(columns[3])
	if ((score >= 0.9) and (pvalue <= 0.05)): #save only contig names were score >=0.9 and pvalue <= 0.05 
		contig_names.append(contig)
		
#open contig file and retrieve only quality checked sequences
selected_sequences = select_sequence("./assembly/final.contigs_5000.fa", contig_names)

#open output file and save information in it
output_file = open("./viral_contig/dvf/dvf.fasta", 'w')

for name, sequence in selected_sequences.items():
    output_file.write(f">{name}\n{sequence}\n")
