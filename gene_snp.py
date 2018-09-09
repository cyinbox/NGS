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
# Library to compute SNP variants from a .vcf file, the vcf file was created by mapping NGS reads 
# of a genome onto a reference genome. 
# Required module library: scikit-allel/https://scikit-allel.readthedocs.io/en/latest/
# pip install scikit-allel
#-----------------------------------------------------------------------------------------------

import numpy as np
import h5py
from pathlib import Path
import allel

#----------------------------------------------------------------------------
# Helper functions 
def getEndItems(data):

 x=getSubList(data)
 se=[]
 seData=[]
 
 for y in x:
   if len(y)>1:
    se.append(y[0])
    se.append(y[-1])
    m=y[-1]-y[0]+1
    se.append(m)
    seData.append(se)
    se=[]

 
 return seData

def getZeros(aList):
 idxs = []
 
 for i in range(len(aList)):
    if aList[i] == 0: 
        idxs.append(i)

 return idxs

#get zero positions
def getLongestZeros(x):
 listM=[]
 
 for i in x:
    listM.append(i[2])
 maxV=max(listM)
 idx=listM.index(maxV)

 zeroV=x[idx]
 
 return [idx,zeroV]

#Test: arrays, each array start, end and length of zeros
'''
x=[[12, 14, 3], [19, 22, 4], [24, 45, 22], [32, 49, 18]]

[idx,zeros]=getLongestZeros(x)
print(idx,zeros)
'''
def splitList(n):
    #return the list index"""
    return [(x+1) for x,y in zip(n, n[1:]) if y-x != 1]

def getSubList(myList):
    
   #split the list base on the index"""
   myIndex = splitList(myList)
   output = list()
   prev = 0
   
   for index in myIndex:
        newList = [ x for x in myList[prev:] if x < index]
        output.append(newList)
        prev += len(newList)
   output.append([ x for x in myList[prev:]])
   return output

#data =[ 1, 4,5,6,7, 8,10, 15,16,17,18, 22, 25,26,27,28]
#print(getSubList(data))

#--------------------------------------------------------------------------
# create variants 
def createVariants(poss,refs,alts):
    
 variants={}
 i=0
 for snp in zip(refs,alts):
    snpx=snp[0]+'->'+snp[1]
    pos=poss[i]
    i=i+1
    variants[pos]=snpx
  
 return variants

#--------------------------------------------------------------------------
# Get the Jaccab distance of two variants
# Inputs: variant, variants2 (positions and SNPs at the positions)
# Output: distance
def getJaccabDistSNP(variantsA,variantsB):
    
 key1=variantsA.keys()
 key2=variantsB.keys()
 key12=list(set().union(key1,key2))

 #count how many SNPs in the same corresponding positions
 commonSNP=0
 for key, value in variantsA.items():
    if key in key2:
      if value==variantsB[key]:
        commonSNP=commonSNP+1
       
 distance=1-commonSNP/len(key12)
 #print('SNP distance=',distance)
 return distance

#Test
'''
poss=[100,102,250,300,302,400,450,1000]
refs=['A','C','G','A','A','C','G','A']
alts=['T','G','C','T','T','G','C','T']

variants=createVariants(poss,refs,alts)
print('Test',variants)

poss2=[100,102,250,300,302,400,450,1000,1001]
refs2=['A','C','G','A','A','C','G','A','G']
alts2=['G','G','C','T','T','G','C','C','A']
variants2=createVariants(poss2,refs2,alts2)

dist=getJaccabDistSNP(variants,variants2)
print('Distance:',dist)
'''
  
#--------------------------------------------------------------------------  
# Convert vcf file to H5
# Input: vcf file name
# output: H5 file with the same name generated in the same folder            
def convertVCFToH5(vcfFileName):
    
 names=vcfFileName.split('.')
 h5FileName=names[0]+'.h5'
 
 h5File = Path(h5FileName)
 
 if h5File.is_file():
   print(' ')
 else:
   print('Convertion')
   allel.vcf_to_hdf5(vcfFileName,  h5FileName, fields='*', overwrite=True)  #The saved data can be accessed via the h5py library, e.g.:


#--------------------------------------------------------------------------  
# get variants from vcf file
# Input: vcf file
# Outputs:callset,variants,variantSNPs
def getVariants(vcfFileName):
    
 convertVCFToH5(vcfFileName) # Need to conver to H5 format to use VariantChunkedTable
 names=vcfFileName.split('.')
 h5FileName=names[0]+'.h5'
 callset = h5py.File(h5FileName, mode='r')
 
 chrom = 'variants'
 variants = allel.VariantChunkedTable(callset[chrom],index='POS')#['variants'], names=['POS', 'REF', 'ALT'],index='POS')
 poss=variants['POS']
 refs=variants['REF']
 alts=variants['ALT'][:, 0]

 variantSNPs={} #make a new format of variants: pos:A->T etc.
 i=0
 
 for snp in zip(refs,alts):
    snpx=snp[0]+'->'+snp[1]
    pos=poss[i]
    i=i+1
    variantSNPs[pos]=snpx 

 return callset,variants,variantSNPs

#--------------------------------------------------------------------------
# Get the distribution density (histogram) of SNP variants 
# Input: callset, window size
# Output: distribution density of SNP variants (y) at each position (x)
def getSNPHistogram(callset, winSize):
    
   pos = allel.SortedIndex(callset['variants/POS'])
   bins = np.arange(0, pos.max(), winSize)
    
   # use window midpoints as x coordinate
   x = (bins[1:] + bins[:-1])/2
    
   # compute variant density in each window
   y, _ = np.histogram(pos, bins=bins)
   #y = y / windowSize
    
   return [x,y]
  
#----------------------------------------------------------------------------
# Represent SNP in two dimensions
