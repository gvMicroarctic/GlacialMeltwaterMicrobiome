#!/bin/bash
#
#SBATCH --job-name=bowtie2    
#SBATCH --output=bowtie2.txt
#SBATCH --ntasks=20
#SBATCH --time=50:00:00
#SBATCH --mem=60000
#SBATCH --mail-user=gilda.varliero@wsl.ch
#SBATCH --mail-type=SUBMIT,END,FAIL

#run mapper
bowtie2 -x ${INDEX} -1 ${FOR} -2 ${REV} -S ${OUT}.sam --end-to-end -X 1000 --threads 40

#convert sam to bam, sort bam and index bam
samtools view -bS ${OUT}.sam | samtools sort -o ${OUT}_si.bam && samtools index ${OUT}_si.bam

#remove sam
rm ${OUT}.sam
