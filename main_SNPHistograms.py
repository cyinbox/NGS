# -*- coding: utf-8 -*-
"""
Created on Sat Aug 25 17:46:00 2018

Changchuan Yin
Dept. of Mathematics, Statistics and Computer Science.
University of Illinois at Chicago
Email: cyin1@uic.edu

Citation:
Yin., C., & Yau., S.S.-T (2018). Whole genome single nucleotide polymorphism genotyping of Staphylococcus aureus. 
Communications in Information and Systems, 2018

"""
#------------------------------------------------------------------------------
# Program to compute SNP variants from a .vcf file
# - Histogram of the SNPs 
# - High and low SNP regions
#-----------------------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import gene_snp as gs

#------------------------------------------------------------------------------
# 0. Input vcf File
#
vcfFileName='SRR5617496.vcf'#:'S.aureus USA300-22862_R1'
names=vcfFileName.split('.')
h5FileName=names[0]+'.h5'
varants_png=names[0]+'.png'
varants_eps=names[0]+'.eps'
varants_pdf=names[0]+'.pdf'

#-----------------------------------------------------------------------------
# 1.Get variants from the vcf file
callset,variants,variantSNPs = gs.getVariants(vcfFileName)

print(callset)
x=sorted(callset.keys())
print('Sorted keys:',x)

chrom = 'variants/CHROM'
print('\n')
print('Chromomes',callset[chrom])

sample=callset['samples']
print('Sample',sample)

#Reference
ref = callset['variants/REF']
print('REF',ref)
refs=str(ref)
print('Length:',len(refs))

qual=callset['variants/QUAL']
print('QUAL',qual)

#------------------------------------------------------------------------------
# 2. Get and plot the hisgrogram of the variants
winSize=250
[pos,SNPs]=gs.getSNPHistogram(callset, winSize=250)

title=str(names) + 'variant histogram, window size 250 bp'
fig, ax = plt.subplots(figsize=(24, 6))
ax.stem(pos,SNPs,'-.')
ax.set_xlabel('chromosome position (bp)',fontsize=20)
ax.set_ylabel('variant count',fontsize=20)
ax.set_title(title)

plt.savefig(varants_png, dpi = 300) 
plt.savefig(varants_eps, dpi = 300) 
plt.show()

#----------------------------------------------------------------------------------
# 3. Get the high SNPs regions
posHighSNPs = []
posLowSNPs = []
threshold  = 10  # When SNP number is larger than threshod 10, the SNP position is recorded
for pos_, SNPs_ in zip(pos, SNPs):
    if SNPs_>= threshold:
      print(int(pos_))
      posHighSNPs.append(pos_)
    else:
      posLowSNPs.append(pos_)  

print('Positions of high SNPs are',posHighSNPs)

#---------------------------------------------------------------------------------
# 4. Results of high and low SNPs 
'''
# Using [pos,SNPs]=gs.getSNPHistogram(callset, winSize=250) for typical MRSA genomes,
# the following regions that have high SNPs after checking the positions of high SNPs
#289625-294125 size: 4500
#405625-407375 size: 1750
#550375-554875 size: 4500

#The following regions are selected for phylogenetic tree (high SNP) since each region is long
1462625-1506875 size: 44250 
1832875-1859875 size: 27000 
1924875-1962375 size: 37500 
2034875-2085375 size: 50500 
'''
'''
# In gene_snp.py file, add these two functions.
def getVariantSNPsH(variantSNPs):
 variantSNPsH = {k: variantSNPs[k] for k in filter(lambda x: 1462625 <= x <= 1506875 or 1832875 <= x <= 1859875 or 1924875 <= x <= 1962375 or 2034875 <= x <= 2085375  , variantSNPs.keys())}
 return variantSNPsH

# Note: The SNPs Lowh regions are defined from file: main_SNPHistograms.py
def getVariantSNPsL(variantSNPs):
 variantSNPsL = {k: variantSNPs[k] for k in filter(lambda x: x<1462625 or 1506875<x<1832875 or 1859875<x<1924875 or 1962375< x <2034875 or x> 2085375, variantSNPs.keys())}
 return variantSNPsL
'''