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
#------------------------------------------------------------------------------
# Program to compute SNP variants from a .vcf file, the vcf file was created by mapping NGS reads 
# of a genome onto a reference genome.
#-----------------------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import gene_snp as gs

#------------------------------------------------------------------------------
# 1.Test vcf file
vcfFileName='SRR4019421.vcf' # S.aureus TX-AML1601921 (MRSA) NGS reads mapped onto S.aureus NCTC8325 (NC_007795)
names=vcfFileName.split('.')
h5FileName=names[0]+'.h5'
varants_png=names[0]+'.png'
varants_pdf=names[0]+'.pdf'

#-----------------------------------------------------------------------------
# 2.Get variants from the vcf file
callset,variants,variantSNPs = gs.getVariants(vcfFileName)

print(callset)
x=sorted(callset.keys())
print('sorted keys:',x)

chrom = 'variants/CHROM'
print('\n')
print('Chromomes',callset[chrom])

sample=callset['samples']
print('test sample',sample)

#Reference
ref = callset['variants/REF']
print('TEST REF',ref)
refs=str(ref)
print('Len1',len(refs))

qual=callset['variants/QUAL']
print('QUAL',qual)

#------------------------------------------------------------------------------
# 3.Get the hisgrogram of the variants
winSize=250
[x,y]=gs.getSNPHistogram(callset, winSize=250)

#-----------------------------------------------------------------------------
# 4. Plot the variant histogram
#title=names + 'variant histogram, window size 250 bp'
fig, ax = plt.subplots(figsize=(24, 6))
ax.stem(x, y,'-.')
ax.set_xlabel('Chromosome position (bp)')
ax.set_ylabel('Variant count')
#ax.set_title(title)

plt.savefig(varants_png, dpi = 300) 