#!/usr/bin/env python3

#open seeker file and retrieve only contigs assigned to phages

#define functions
def read_fasta(file_path):
    """Reads a FASTA file and returns a dictionary of sequences."""
    sequences = {}
    current_sequence = ""
    
    with open(file_path, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
            	current_sequence = line[1:].split()
            	sequences[current_sequence[0]] = ""
            else:
                sequences[current_sequence[0]] += line
    
    return sequences

def select_sequence(fasta_file_path, selected_sequence_names):
    """Selects sequences from a FASTA file based on a list of sequence names."""
    all_sequences = read_fasta(fasta_file_path)
    
    selected_sequences = {name: sequence for name, sequence in all_sequences.items() if name in selected_sequence_names}
    
    return selected_sequences
    
    
#open input seeker file
dvf_file = open("viral_contig/seeker/seeker_trimmed.txt", "r")

#define empty list to save contig names
contig_names = []

for line in dvf_file: #iterate through the file line by line
	columns = line.split("\t") #split in different columns
	contig = columns[0]
	definition = columns[1]
	score = float(columns[2])
	if ((definition == "Phage") and (score >= 0.9)): #save only contig names that were assigned to "Phage" 
		contig_names.append(contig)
		
#open contig file and retrieve only quality checked sequences
selected_sequences = select_sequence("./assembly/final.contigs_5000.fa", contig_names)

#open
output_file = open("./viral_contig/seeker/seeker.fasta", 'w')

for name, sequence in selected_sequences.items():
    output_file.write(f">{name}\n{sequence}\n")
