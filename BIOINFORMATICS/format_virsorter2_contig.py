#!/usr/bin/env python3

#open virsorter2 file and retrieve only contigs with max_score >= 0.7

#define functions
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
    
    
#open input virsorter file
dvf_file = open("./viral_contig/virsorter2/final-viral-score.tsv", "r")

#define empty list to save contig names
contig_names = []
next(dvf_file) #skip first line
for line in dvf_file: #iterate through the file line by line
	columns = line.split("\t") #split in different columns
	contig = columns[0]
	max_score = float(columns[3])
	if (max_score >= 0.9): #save only contig names that were assigned to "Phage" 
		contig_names.append(contig)
		
#open contig file and retrieve only quality checked sequences
selected_sequences = select_sequence("./viral_contig/virsorter2/final-viral-combined.fa", contig_names)

#open
output_file = open("./viral_contig/virsorter2/virsorter2.fasta", 'w')

for name, sequence in selected_sequences.items():
    output_file.write(f">{name}\n{sequence}\n")
