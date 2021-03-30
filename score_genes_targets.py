#!/usr/bin/python
# -*- coding: utf-8 -*-
# 25/02/15
# score_genes_targets.py version 1.0, part of the Regulus suite
# python 2.7
# Copyright Â© 2015 IBENS/Dyogen : Magali Naville and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or magali.naville@ens-lyon.fr
# Related Publication: Naville et al. (2015) Nature Communications
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France


# Calculate the association score of each gene found in the environment of a family of orthologous CNEs.

### Usage ###
# score_genes_targets.py genes_env_chr*.list.gz

import gzip
import sys
import bz2

###############
### Functions ###
###############

## Function that sorts vectors according to the values of column 1 and then column 3 (distance CNE-gene) ##
def cmpval(x,y):
    if float(x[1])>float(y[1]):
        return -1
    elif float(x[1])==float(y[1]):
        if float(x[3])<float(y[3]):
	    return -1
	elif float(x[3])==float(y[3]):
	    return 0
	else:
	    return 1
    else:
        return 1

###########
### Data ###
###########

## Weighting factors depending of the degree of rearrangement of each species genome compared to the human one. ##
file = open('species_synteny.list','r')
synt = file.readlines()
file.close()

SY = {} # taux de syntenie
COUV = {} # couverture 
for s in range(1,len(synt)) :
	ls = synt[s].split('\t')
	SY[ls[11].strip()] = ls[5].replace(',','.')
	COUV[ls[11].strip()] = ls[7].replace(',','.')

## Names of genes in Human corresponding to ancestral genes ##
NEW = {}
file = bz2.BZ2File('ancGenes.Boreoeutheria.list.bz2','r')
new = file.readlines()
file.close()

for i in new :
	li = str(i).split()
	NEW[li[0]] = []
	for k in range(1,len(li)):
		if li[k][0:5] == 'ENSG0' :
			NEW[li[0]].append(li[k])

## Current gene names ##
COUR = {} 
file = bz2.BZ2File('names.Homo.sapiens.list.bz2','r')
I = file.readlines()
file.close()

for i in I :
	li = i.strip().split()
	COUR[li[0]] = li[1]

## Genes found on the considered chromosome ##
X = []
file = bz2.BZ2File('genesST.Homo.sapiens.list.bz2','r')
genes = file.readlines()
file.close()

for i in genes :
	li = i.split()
	if li[0] == sys.argv[1].split('chr')[1].split('.')[0] :
		X.append(li[4])
print sys.argv[1].split('chr')[1].split('.')[0]

###########
### Main ###
###########

file = gzip.open(sys.argv[1],'r') 
cne = file.readlines()
file.close()

CNE = {} # CNE dictionary 
KEYS = [] # CNE identifiers
SIZE = {} # dictionary of the number of species in which each CNE is present
SCORE = {} # Dictionary CNE/genes/[scores0 scores1 scores2 dist_moy_scores1 Nbre_1_outgroups]
test_tar = 0 # Test to check that at least one CNE has a putative target
LIMIT = 1000000

for c in cne :
  lc = c.split('\t')
  #print lc
  if len(lc)>1 : # pour eviter les '\n'
    if lc[1] != 'choHof1' and lc[1] != 'eriEur1' and lc[1] != 'macEug1' and lc[1] != 'sorAra1' and lc[1] != 'tarSyr1' and lc[1] != 'hg19' and lc[1] != 'petMar1' : # Species of too low genome coverage, not taken into account + human (common denominator)
	if lc[0] not in CNE : # CNE not referenced in the dictionary yet
		CNE[lc[0]] = {}
		SIZE[lc[0]] = lc[2]
		KEYS.append(int(lc[0].split('U')[1]))
	CNE[lc[0]][lc[1]] = lc[3:len(lc)-1] # CNE[SCNE#][esp]

KEYS.sort()

for n in  KEYS : # for each CNE
	c = 'U'+str(n)
	SCORE[c] = {} # dictonary of each gene score for a given CNE
	for i in CNE[c] : # pour chaque espece i
		for j in CNE[c][i] :
			lj = j.split('-')
			if lj[1] == "0" :
				score = 0
				dist = 0
			else :
				test_tar = 1
				score = 0
				dist = LIMIT
				for k in range(1,len(lj)) :
					lk = lj[k].split('_')
					if lk[1] == '1' :
						score = 1
						if int(lk[2]) < dist :
							dist = lk[2]
					if lk[1] == '2' and (score == 0 or score == 3) :
						score = 2
					if lk[1] == '3' and score == 0 :
						score = 3
			if SCORE[c].__contains__(lj[0]) == 0 :
				SCORE[c][lj[0]] = [0,0,0,0,0,0] 
			if score == 1 :
				SCORE[c][lj[0]][1] += float(SY[i])
				SCORE[c][lj[0]][3] += int(dist)
				SCORE[c][lj[0]][4] += 1
			elif score == 2 :
				SCORE[c][lj[0]][2] += 1/float(SY[i])
			else : # score de 0 ou de 3
				SCORE[c][lj[0]][0] += float(COUV[i])/float(SY[i])

LIGNES = []

if test_tar == 1 :

	file = gzip.open('score_genes_'+sys.argv[1].split('env_')[1].split('.')[0]+'.list.gz','w')

	for n in  KEYS :
		c = 'U'+str(n)
		file.write(c+'\t'+str(SIZE[c])+'\t')
		ligne = c+'\t'+str(SIZE[c])+'\t'
		UN = [] # genes presenting a number of "1" scores different from 0
		for i in SCORE[c] :
			if SCORE[c][i][1] != 0 :
				UN.append([i,float(SCORE[c][i][1])-float(SCORE[c][i][2])-float(SCORE[c][i][0]),SCORE[c][i][1],SCORE[c][i][2],SCORE[c][i][0],float(SCORE[c][i][3])/float(SCORE[c][i][4]),SIZE[c]])

		UN.sort(cmpval)
		for u in UN :
			file.write(u[0]+' %.2f %.0f ' %(float(u[1]), float(u[5]))) # absolute score (1-2-3-0) mean distance
			ligne = ligne + u[0]+' %.2f %.0f ' %(float(u[1]), float(u[5]))
			names = ''
			for j in NEW[u[0]] :
				if X.__contains__(j) :
					if j in COUR:
						names = names + COUR[j] + '/'
					else:
						names = names + j + '/'
			
			file.write(names[0:len(names)-1]+' ')
			ligne = ligne + names[0:len(names)-1]+' '
		file.write('\n')
		ligne = ligne + '\n'
		
		LIGNES.append(ligne)
	file.close()

	file = gzip.open('targets_'+sys.argv[1].split('env_')[1].split('.')[0]+'.list.gz','w')
	
	for ligne in LIGNES :
		ls = ligne.split()
		#print ls
		gen = []
		GEN = {}
		if len(ls) > 2 :
			i = 2
			while i < len(ls) :
				gen.append(ls[i])
				GEN[ls[i]] = [ls[i+1],ls[i+2],ls[i+3]] # absolute score, mean distance, human gene name
				i += 4
			tar = []
			n = 0

			while GEN[gen[n]][0] == GEN[gen[0]][0] :
				tar.append(gen[n])
				if n < len(gen)-1 :
					n += 1
				else :
					break
			
			file.write(ls[0]+'\t'+ls[1]+'\t')
			if len(tar) == len(gen) :
				file.write('-\t')
			else :
				file.write(str(float(GEN[gen[n-1]][0])-float(GEN[gen[n]][0]))+'\t')
			file.write(GEN[gen[0]][0]+'\t')
			for t in tar :
				file.write(GEN[t][2]+' ')
			file.write('\n')
		# to also indicate elements with no target
		else :
			file.write(ls[0]+'\t'+ls[1]+'\t-\t-\t-\n')
		
	file.close()
