# -*- coding: utf-8 -*-
"""
Created on Sat Aug 25 17:46:00 2018

Changchuan Yin, Ph.D.
University of Illinois at Chicago
Email: cyin1@uic.edu

Citation:
Yin., C., & Yau., S.S.-T (2018). Whole genome SNP genotyping of methicillin-resistant Staphylococcus aureus. 
Communications in Information and Systems, 2018

"""
#------------------------------------------------------------------------------------------------
# Program to compute SNP variants from a .vcf file, the vcf file was created by mapping NGS reads 
# of a genome onto a reference genome. This test is to show deletion variant in the testing genome.
#-----------------------------------------------------------------------------------------------


import matplotlib.pyplot as plt
import gene_snp as gs

#-------------------------------------------------------------------------------------------------
# 1. Input vcf file
vcfFileName='SRR5714648.vcf' # S.aureus H41 NGS (MSSA) NGS reads mapped onto S.aureus NCTC8325 (NC_007795)
names=vcfFileName.split('.')
h5FileName=names[0]+'.h5'
varants_png=names[0]+'.png'
varants_pdf=names[0]+'.pdf'

#-------------------------------------------------------------------------------------------
# 2. Get variants from the vcf file
callset,variants,variantSNPs = gs.getVariants(vcfFileName)

print(callset)
x=sorted(callset.keys())
print('sorted keys:',x)

chrom = 'variants/CHROM'
print('\n')
print('Variant Chrom',callset[chrom])

sample=callset['samples']
print('Samples:',sample)

#Reference
ref = callset['variants/REF']
print('variant REF:',ref)
refs=str(ref)
print('Length reference',len(refs))

qual=callset['variants/QUAL']
print('Variant QUAL',qual)

#-------------------------------------------------------------------------------------------
# 3. Compute the histogram of the variants
winSize=250
[x,y]=gs.getSNPHistogram(callset, winSize=250)

#-------------------------------------------------------------------------------------------
# 4. Plot the variant histogram of variants 
#title='SRR5617496 variant histogram, window size 250 bp'
fig, ax = plt.subplots(figsize=(24, 6))
ax.stem(x, y,'-.')
ax.set_xlabel('Chromosome position (bp)')
ax.set_ylabel('Variant count')
#ax.set_title(title)
plt.savefig(varants_png, dpi = 300) 

#------------------------------------------------------------------------------------------
# 5. Get the zeros regions (deletions) of the variant histograms
idxZeros = gs.getZeros(y) #The indices of zero element
print('Zero positions in reference genome:\n',idxZeros)

seData=gs.getEndItems(idxZeros)
print('Zero positions mapped on the reference geneome [start, end, length]\n',seData)  

[idx,zeros]=gs.getLongestZeros(seData)
print('Zeros:',zeros)
startP=zeros[0]*winSize
endP=zeros[1]*winSize
length=endP-startP
print(str(startP),str(endP))
print('Longest zeros (deletion) in genome position start at',str(startP),' and end at:',str(endP), 'length:',length)
