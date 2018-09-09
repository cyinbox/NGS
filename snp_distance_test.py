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
#-----------------------------------------------------------------------------------------------
# Program to compute SNP variants from a .vcf file, the vcf file was created by mapping NGS reads 
# of a genome onto a reference genome. This test is to compute the pair-wise Jaccard distances of 
# variants.
#-----------------------------------------------------------------------------------------------


import gene_snp as gs

#----------------------------------------------------------------------------------------------
# 1.Test genome vcf files
'''
genes={'SRR7295360.vcf':'S.aureus HPV107','SRR7295358.vcf':'S.aureus NRS2','ERR2541869.vcf':'S.aureus ST88','ERR2275539.vcf':'S.aureus MSSA-HWC2014',\
      'SRR5714648.vcf':'S.aureus MSSA-H41','SRR6675978.vcf':'S.aureus MSSA-ST97','ERR2541870.vcf':'S.aureus Nigeria','SRR1014703.vcf':'S.aureus USA400-BAA1752 CC1-ST1',\
      'SRR1014705.vcf':'S.aureus CC80-24329-ST153','SRR1014706.vcf':'S.aureus USA700-NRS386 CC72-ST72','SRR1014707.vcf':'S.aureus USA700-GA-442 CC72-ST72','SRR1014709.vcf':'S.aureus USA800-NY-313 CC5-ST83',\
      'SRR1014711.vcf':'S.aureus USA100-OR-293 CC5-ST5','SRR1014712.vcf':'S.aureus USA100-NY-76 CC5-ST5','SRR1014713.vcf':'S.aureus USA100-NRS382 CC5-ST5','SRR1014714.vcf':'S.aureus USA100-NY-54 CC5-ST105',\
      'SRR1014715.vcf':'S.aureus USA100-CA-126 CC5-ST5','SRR1014716.vcf':'S.aureus USA100-CA-248 CC5-ST5','SRR1014717.vcf':'S.aureus USA1000-CA-629 CC59-ST87','ERR377327.vcf':'S.aureus USA300-RU-CAMP-29c',\
      'ERR2276455.vcf':'S.aureus USA300-RU-CAMP-P29a','ERR2276459.vcf':'S.aureus USA300-RU-CAMP-P29b','ERR2541868.vcf':'S.aureus ST88-a','ERR2541871.vcf':'S.aureus ST88-b',\
      'SRR6304957.vcf':'S.aureus PtA04-T3','SRR6304955.vcf':'S.aureus PtA02-T1','SRR5617496.vcf':'S.aureus USA300-22862_R1','SRR1014694.vcf':'S.aureus USA300-R-54',\
      'SRR1014708.vcf':'S.aureus USA800-NRS387 CC5-ST5','SRR4019421.vcf':'S.aureus TX-AML1601921','SRR1014724.vcf':'S.aureus USA600-CA-347 CC45-ST45','SRR1014722.vcf':'S.aureus USA600-BAA1754 CC45-ST45',\
      'SRR1014700.vcf':'S.aureus USA500-NRS385E','SRR1014720.vcf':'S.aureus USA200-NRS383 CC30-ST346','ERR2562460.vcf':'S.aureus CC398-ST899','SRR6475664.vcf':'S.aureus ST398','ERR1795563.vcf':'S.aureus MSSA-SHAIPI'}
'''
genes={'SRR7295360.vcf':'S.aureus HPV107','SRR7295358.vcf':'S.aureus NRS2','ERR2541869.vcf':'S.aureus ST88','ERR2275539.vcf':'S.aureus MSSA-HWC2014'}

n=len(genes)
keys=genes.keys()
keys=list(keys)
names=genes.values()
names=list(names)

# 2.Compute the pair-wise Jaccard distances of variants, the output is a text file (distance matrix), distances.txt
# Note: this distance file is the input to the MATLAB program, phylogeneticTreeByDistances.m, for phylogenetic tree build.
# https://github.com/cyinbox/NGS/blob/master/phylogeneticTreeByDistances.m

file = open('distances.txt.txt', 'w')

if len(names) > len(set(names)):
   print('Not unique strain names!')  # Only unique names are allowed in phylogenetic tree

else:  
   for i in range(0,n):
     keyA=keys[i]
     nameA=names[i]
     callsetA,variantsA,variantSNPsA = gs.getVariants(keyA)
  
     for j in range(i+1,n):
       keyB=keys[j]
       nameB=names[j]
       callsetB,variantsB,variantSNPsB = gs.getVariants(keyB)
   
       dist=gs.getJaccabDistSNP(variantSNPsA,variantSNPsB)
       print(keyA,keyB,dist)
       strDist=nameA+','+nameB+','+str(dist)+'\n'
       file.write(strDist)

file.close()
