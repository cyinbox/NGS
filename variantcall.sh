#!/bin/bash
# get SRA:./fastq-dump --split-files SRR5617496
# Example usage: ./variantcall2.sh NC_007795.fasta SRR5617496_1.fastq SRR5617496_2.fastq
# dos2unix -iso -n variantcall.sh variantcall2.sh

#Step 1: index the reference genome
bwa index $1

#Step 2: use BWA-MEM to align paired-end sequences. Briefly, the algorithm works by seeding alignments with
bwa mem -M -t 16 $1 $2 $3 > genome.sam

#picard uses faiindex and dictionary (steps 3-7) to access and safety check access to the reference files
#Step 3: index the reference genome again for faindex, it allow efficient random access to the reference bases
samtools faidx $1

#Step 4: converte the sam file to a bam file for fast processing speed
samtools view -S -b genome.sam > genome.bam

#Step 5: sort the bam file
samtools sort genome.bam > genome_sorted.bam

#Step 6: index the sorted bam file
samtools index genome_sorted.bam

#Step 7: create fasta sequence dictionary: the dictionary of the contig names and sizes 
java -jar picard.jar CreateSequenceDictionary REFERENCE=$1 OUTPUT=NC_007795.dict

#Step 8: add or read groups into the sorted bam file
java -jar picard.jar AddOrReplaceReadGroups  I=genome_sorted.bam O=genome_sortedG.bam  RGID=4 RGLB=PtA02-T1 RGPL=IlluminaMiSeq  RGPU=unit1  RGSM=19

#Step 9: build index on the sorted bam file
java -jar picard.jar BuildBamIndex I=genome_sortedG.bam

#Step 10: variant calling by GATK
java -jar GenomeAnalysisTk.jar -T UnifiedGenotyper -R $1 -I genome_sortedG.bam  -o genome.vcf -filterMBQ

#Step 11: coverage analysis
bedtools genomecov -ibam genome.bam -d>genome_coverage.txt

#Step 12: delete temp files
rm genome.sam
rm genome.bam
rm genome_sorted.bam
rm genome_sortedG.bam
#rm $2
#rm $3
