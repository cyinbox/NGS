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
#-----------------------------------------------------------------------------------------------
# Program to compute the pair-wise distances between SNP variants from .vcf files, the vcf files
# were created by mapping NGS reads of a genome onto a reference genome.  
# The distance is for constructing a phylogenetic tree.
#-----------------------------------------------------------------------------------------------

import gene_snp as gs

#----------------------------------------------------------------------------------------------
# 1. The VCF files for computing the distances. VCF file name: genome name
genes = {'SRR7295360.vcf':'S.aureus HPV107','SRR7295358.vcf':'S.aureus NRS2','ERR2541869.vcf':'S.aureus ST88','ERR2275539.vcf':'S.aureus MSSA-HWC2014'}

n=len(genes)
keys = genes.keys()
keys = list(keys)
names = genes.values()
names = list(names)
print(names)

# 2. Compute the pair-wise Jaccard distances of variants, the output is a text file (distance matrix), distances.txt
# Note: this distance file is the input to the MATLAB program, phylogeneticTreeByDistances.m, for phylogenetic tree build.
# https://github.com/cyinbox/NGS/blob/master/phylogeneticTreeByDistances.m

#------------------------------------------------------------------------------
# 2.1 High SNP distances
file = open('distancesHighSNPs.txt', 'w')

if len(names) > len(set(names)):
   print('Not unique strain names!')  # Only unique names are allowed in phylogenetic tree

else:  
   for i in range(0,n):
     keyA=keys[i]
     nameA=names[i]
     print(nameA)
     callsetA,variantsA,variantSNPsA = gs.getVariants(keyA)
     variantSNPsH_a = gs.getVariantSNPsH(variantSNPsA) 
   
     for j in range(i+1,n):
       keyB = keys[j]
       nameB = names[j]
       print(nameB)
       callsetB,variantsB,variantSNPsB = gs.getVariants(keyB)
       variantSNPsH_b = gs.getVariantSNPsH(variantSNPsB)
      
       dist=gs.getJaccabDistSNP(variantSNPsH_a,variantSNPsH_b)
       print(keyA,keyB,dist)
       strDist = nameA+','+nameB+','+str(dist)+'\n'
       file.write(strDist)
       
file.close()
print('High SNP distances computed successfully.')

#------------------------------------------------------------------------------
# 2.2. Low SNP distances
file = open('distancesLowSNPs.txt', 'w')

if len(names) > len(set(names)):
   print('Not unique strain names!')  # Only unique names are allowed in phylogenetic tree

else:  
   for i in range(0,n):
     keyA = keys[i]
     nameA = names[i]
     print(nameA)
     callsetA,variantsA,variantSNPsA = gs.getVariants(keyA)
     variantSNPsL_a = gs.getVariantSNPsL(variantSNPsA) 
   
     for j in range(i+1,n):
       keyB = keys[j]
       nameB = names[j]
       print(nameB)
       callsetB,variantsB,variantSNPsB = gs.getVariants(keyB)
       variantSNPsL_b = gs.getVariantSNPsL(variantSNPsB)
      
       dist = gs.getJaccabDistSNP(variantSNPsL_a,variantSNPsL_b)
       print(keyA,keyB,dist)
       strDist = nameA+','+nameB+','+str(dist)+'\n'
       file.write(strDist)
       
file.close()
print('Low SNP distances computed successfully.')


