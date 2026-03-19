#!/usr/bin/env python3

#merge all viral contigs

#define functions
def read_fasta_seq(file_path):
    """Reads a FASTA file and returns sequences."""
    sequences = []
    
    with open(file_path, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if not line.startswith('>'):
                sequences.append(line)
    
    return sequences
    
#open contig file and retrieve only quality checked sequences - DNA
sequences_virsorter_dna = read_fasta_seq("./viral_contig/virsorter2/virsorter2.fasta")  
sequences_vibrant_dna = read_fasta_seq("./viral_contig/vibrant/vibrant.fasta")  
sequences_dvf_dna = read_fasta_seq("./viral_contig/dvf/dvf.fasta")  
sequences_seeker_dna = read_fasta_seq("./viral_contig/seeker/seeker.fasta") 

#merge the contig sequences into one list
sequences = sequences_virsorter_dna + sequences_vibrant_dna + sequences_dvf_dna + sequences_seeker_dna
unique_sequences = list(dict.fromkeys(sequences))

#save new unique contigs with new contig names
output_file = open("./viral_contig/viral_contigs.fasta", 'w')
name = 0
for sequence in unique_sequences:
	name +=1
	output_file.write(f">contig_{name}\n{sequence}\n")
