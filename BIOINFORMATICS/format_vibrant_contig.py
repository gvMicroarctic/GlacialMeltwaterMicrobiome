#!/usr/bin/env python3

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

#read as header - sequence the vibrant output
all_sequences = read_fasta("./viral_contig/vibrant/VIBRANT_final.contigs_5000/VIBRANT_phages_final.contigs_5000/final.contigs_5000.phages_combined.fna")
    
#open
output_file = open("./viral_contig/vibrant/vibrant.fasta", 'w')

#print sequences into file
for contig in all_sequences:
    output_file.write(f">{contig}\n{all_sequences[contig]}\n")
