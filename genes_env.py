#!/usr/bin/python
# -*- coding: utf-8 -*-
# 25/02/15
# genes_env.py 1.0, part of the Regulus suite 
# python 2.7
# Copyright © 2013 IBENS/Dyogen : Magali Naville and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or magali.naville@ens-lyon.fr
# Related Publication: Naville et al. (2015) Nature Communications 
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France


# Lists all the genes that are neighbour of each CNE in Human, along with their 'status' in the other species possessing the considered CNE.

### Usage ###
# genes_env.py chrX Examples/fam_chrX.list 

import math
import gzip
import bz2
import sys

#################
### Functions ###
#################

## Reads a bed file containing gene coordinates (in gz or bz2 format) and put them in a dictionary. ##

def GetCoord(fichier, dico, esp, chro):
	if fichier[len(fichier)-1] == '2' :
		file = bz2.BZ2File(fichier,'r')
	else :
		file = gzip.open(fichier,'r')
	coo = file.readlines()
	file.close()
	
	if esp not in dico :
		dico[esp] = {}
	
	for c in coo :
		lc = c.strip().split()
		if chro == '' :
			dico[esp][lc[4]] = [lc[0],int(lc[1]),int(lc[2])] # chr, deb, fin
		elif lc[0] == chro.split('chr')[1] or lc[0] == chro :
			dico[esp][lc[4]] = [lc[0],int(lc[1]),int(lc[2])]

## Reads a bed files containing CNE coordinates (in gz or bz2 format) and put them in a dictionary. ## 

def GetCoordCNE(fichier, dico, esp):
	if fichier[len(fichier)-1] == '2' :
		file = bz2.BZ2File(fichier,'r')
	else :
		file = gzip.open(fichier,'r')
	coo = file.readlines()
	file.close()
	
	if esp not in dico :
		dico[esp] = {}
	
	for c in coo :
		lc = c.strip().split()
		dico[esp][lc[4].split('_')[1]] = [lc[0].replace('chr',''),int(lc[1]),int(lc[2])] # chr, deb, fin


###########
### Data ###
###########

CHR = sys.argv[1]

## Lists elements for which we search target genes ##
file = open(sys.argv[2],'r')
scne = file.readlines()
file.close()

SCNE = {}
for s in scne :
  ls = s.strip().split()
  if ls[1][0:2] == 'hg' :
    SCNE['U'+ls[0]] = ls[1:len(ls)]


## CNEs coordinates ##
CNE = {}

## Boreoeutheria ancestral genes + current genes in the different species ##
file = bz2.BZ2File('ancGenes.Boreoeutheria.list.bz2','r')
anc = file.readlines()
file.close()

ANC = {}
for a in anc :
  la = a.strip().split()
  for i in la :
    if i[0:5] == 'ENSG0' :
    	ANC[i] = la 
      
## Lists genes homologous to the considered human gene (up to Coelomata).##
file = bz2.BZ2File('ancGenes.Coelomata.list.bz2','r')
homo = file.readlines()
file.close()

for h in homo :
	lh = h.strip().split()
	for i in lh :
		if i[0:5] == 'ENSG0' :
			if i in ANC :
				for j in range(1,len(lh)) :
					if lh[j] not in ANC[i] :
						ANC[i].append(lh[j])

## Human gene coordinates ##
COO = {}
GetCoord('genesST.Homo.sapiens.list.bz2', COO, 'hg19', CHR)

## Species list ##
file = open('46_species2.list','r')
liste_esp = file.readlines()
file.close()

LIST_SPECIES = {}
for l in liste_esp :
	ll = l.strip().split()
	LIST_SPECIES[ll[0]] = ll[1]

## Weighting factors and maximum distance allowed for each species (functions of genome size and assembly quality) ##
file = open('46_species_sizes.list','r')
size = file.readlines()
file.close()

SIZE = {}
for s in size :
	ls = s.strip().split('\t')
	SIZE[ls[6]] = int(1000000*float(ls[4]))

###########
### Main ###
###########

OUT = ''

## Bed of CNEs in Human ##
file = gzip.open('hg19_'+CHR+'.bed.gz')
hum = file.readlines()
file.close()

file = gzip.open('genes_env_'+CHR+'.list.gz','w')

for h in hum :
	lh = h.strip().split()
	el = lh[4].split('_')[1] # CNE name (without species name)
	UNIQ = {} # List of ancestral genes already analyzed, to avoid redundancies linked to paralogous modern genes.
	
	# Species presenting an orthologous CNE :
	SPECIES = []
	for i in SCNE[el] :
		if len(i.split('_')) == 2 :
			if i.split('_')[0] in LIST_SPECIES :
				SPECIES.append(i.split('_')[0])
	
	# CNE and gene coordinates of corresponding species.
	for s in SPECIES :
		if s not in CNE and s in LIST_SPECIES :
			GetCoordCNE(s+'_'+CHR+'.bed.gz',CNE,s)
			GetCoord('genesST.'+LIST_SPECIES[s]+'.list.bz2',COO,s,'')
	
	# Neighbouring genes in Human
	
	TEST = {} # 1 if the gene is close to the CNE in at least one species, 0 otherwise. 
	
	for j in COO['hg19'] :
		if COO['hg19'][j][0]==CNE['hg19'][el][0] and (math.fabs(CNE['hg19'][el][1]-COO['hg19'][j][2])<=1200000 or (CNE['hg19'][el][1]>=COO['hg19'][j][1] and CNE['hg19'][el][2]<= COO['hg19'][j][2]) or math.fabs(COO['hg19'][j][1]-CNE['hg19'][el][2])<=1200000) : # test de distance
			if j in ANC : # if the gene possess an ancestor in Boreorutheria
				if ANC[j][0] not in UNIQ : # Test to check if ancestral gene has not been analyzed yet
					UNIQ[ANC[j][0]] = {}
					TEST[ANC[j][0]] = 0
					for k in range(1,len(ANC[j])): # Look at each modern gene
						for s in SPECIES : # Test to analyze only genes from species presenting the CNE.
							if ANC[j][k] in COO[s] and el in CNE[s] :
								if s not in UNIQ[ANC[j][0]] :
									UNIQ[ANC[j][0]][s] = []
								
								if COO[s][ANC[j][k]][0] != CNE[s][el][0] : # if different chromosomes 
									UNIQ[ANC[j][0]][s].append(ANC[j][k]+'_3')
								elif COO[s][ANC[j][k]][0] == CNE[s][el][0] :
									if (CNE[s][el][1]>=COO[s][ANC[j][k]][1] and CNE[s][el][2]<=COO[s][ANC[j][k]][2]) or (CNE[s][el][1]<COO[s][ANC[j][k]][1] and CNE[s][el][2]>=COO[s][ANC[j][k]][1]) or (CNE[s][el][1]<=COO[s][ANC[j][k]][2] and CNE[s][el][2]>COO[s][ANC[j][k]][2]) :
										UNIQ[ANC[j][0]][s].append(ANC[j][k]+'_1_0')
										TEST[ANC[j][0]] = 1
									elif CNE[s][el][2]<COO[s][ANC[j][k]][1] : # if CNE_end < gene_start
										dist = COO[s][ANC[j][k]][1]-CNE[s][el][2]
										if dist <= SIZE[s] : # distance lower than the threshold
											UNIQ[ANC[j][0]][s].append(ANC[j][k]+'_1_'+str(dist))
											TEST[ANC[j][0]] = 1
										else : # distance greater than the threshold
											UNIQ[ANC[j][0]][s].append(ANC[j][k]+'_2')
									else : # if CNE_start > gene_end 
										dist = CNE[s][el][1]-COO[s][ANC[j][k]][2]
										if dist <= SIZE[s] :
											UNIQ[ANC[j][0]][s].append(ANC[j][k]+'_1_'+str(dist))
											TEST[ANC[j][0]] = 1
										else :
											UNIQ[ANC[j][0]][s].append(ANC[j][k]+'_2')
	
	for s in SPECIES :
		if el in CNE[s] :
			file.write(el+'\t'+s+'\t'+str(len(SPECIES))+'\t')
			for u in UNIQ :
				if TEST[u] == 1 :
					file.write(u+'-')
					if s not in UNIQ[u] :
						file.write('0\t')
					else :
						out = ''
						for p in UNIQ[u][s] :
							out = out+p+'-'
						file.write(out[0:len(out)-1]+'\t')
		file.write('\n')
	
	file.write('\n')

file.close()
