#!/usr/bin/python
# -*- coding: utf-8 -*-
# Scans a maf multiple alignement to search regions of particular conservation.

import sys
import re
import os
from optparse import OptionParser
from bisect import bisect_left
import gzip

Champs=[]
score=[]
col=[]
colproc=[]  # processed column
currwin=[]  # current window of length W
lstat=[]
extsart=[]
exend=[]
fused=[]
scorid=[]
lis1=[]
lis2=[]
y=0
labnr=0

CD=os.getcwd()

annot={}
block={}

n=0
ID=0
z=0
t=0 

##############
# Parameters #
##############

usage = "usage: %prog [MAF file] [options]"

##########
# Module #
##########

# Module to read options
parser = OptionParser(usage=usage)

parser.add_option("-w", dest="W", default="50", help="Minimum size of anchor window to initiate a search [default = 50 bp]")
parser.add_option("-x", dest="X", default="3", help="Max cutoff for consecutive column mismatches [default = 3 columns]")
parser.add_option("-i", dest="ID", default="90", help="Minimum %ID of alignment in anchor window [default = 90 %]")
parser.add_option("-s", dest="MAS", default="300", help="Minimum score of a block in the MAF file before it is considered for a search [default = 300]")
parser.add_option("-f", dest="F", default="0", help="Filters out output profiles that contain low complexity sequences and/or that fall below the length provided as argument and/or that overlap small cap nucleotides (repeat masked) in one of the species given in the -c option. Low complexity means that at least one base has a frequency below 5%. [default = 0 bp]")
parser.add_option("-l", dest="L", default="80", help="Threshold for low complexity sequences filtering; sequences containing more thant L% of the same nucleotide are discarded [default = 80%]")
parser.add_option("-t", dest='T',default="60", help="Minimum similarity of a column to be considered as conserved [default = 60 %]")
parser.add_option("-n", dest='N',default="6", help="Minimum number of species required in the alignment block to be considered [default = 6]")
parser.add_option("-c", dest="CAT", default="", help="Catalogue of species (at least 2) to consider when searching the MAF file for the regions of interest. Argument is a coma separated list of species, with species names as indicated in the MAF file. [compulsory ; default = None]")
parser.add_option("-a", dest="annot", default="", help="If this option is set, a directory called SM_annot must be in the current directory with annotation files in UCSC format (>10 tab separated columns) called species_annot.ucsc, where species is replaced by the species names as indicated in the MAF file. Annotation coordinates are O-based, in UCSC fashion. The argument is a coma separated list of species for which annotations files are provided in the SM_annot directory. When this option is set, no output will be produced in regions of the multiple alignements of the MAF file that overlap the annotations in these genomes. Generally these annotations are non redundant (e.g. ensembl genes) and will be used as provided. This option can be combined with -r for other species.  [default = none].")
parser.add_option("-r", dest="reduced", default="", help="If this option is set, a directory called SM_annot must be in the current directory with annotation files in rr format. This format is produced by the remred.py script, which removes redundancy in large redundant annotation files such as mammalian EST alignment files. Annotation coordinates are O-based, in UCSC fashion. The argument is a coma separated list of species for which annotations files are provided in the SM_annot directory. When this option is set, no output will be produced in regions of the multiple alignements of the MAF file that overlap the annotations in these genomes. This option can be combined with -a for other species. In general, using this option instead of -a reduces computation time [default = none]")

(options, args) = parser.parse_args()

W=int(options.W)
ID=int(options.ID)
X=int(options.X)
MAS=int(options.MAS)
#M=str(options.M)
F=int(options.F)
L=int(options.L)
T=int(options.T)
N=int(options.N)
esp0 = ''

## Module to parse the catalogue of species

if len(options.CAT) == 0:
    print >> sys.stderr, "ERROR: Please provide at least two valid species names"
else:
    scat= str(options.CAT).split(",")
    lsc=len(scat)

# Module to read annotation files in UCSC format
if len(options.annot) > 0:
    lis1=str(options.annot).split(",")
    for sp in lis1:
        if sp not in scat:
            print >> sys.stderr, "ERROR: species", sp, "not in species catalogue"
            sys.exit(0)
        else:
            fle= CD+"/SM_annot/"+str(sp)+"_annotation.ucsc"
            print >> sys.stderr, "Reading annotation file", fle
            for ligne in open(fle,'r').xreadlines():
                champs=ligne.split("\t")
                cle=sp + "." + champs[0]
                annot[cle]=eval(champs[1])
	    print >> sys.stderr," Done... "    

# Module to read annotation files in 'rr' (non redundant) format
if len(options.reduced) > 1:
    lis2=str(options.reduced).split(",")
    for sp in lis2:
        if sp not in scat:
            print >> sys.stderr, "ERROR: species", sp, "not in species catalogue"
            sys.exit(0)
        else:
	    fle= CD+"/SM_annot/"+str(sp)+"_annotation.rr"
            print >> sys.stderr, "Reading non redundant annotation file", fle
            for ligne in open(fle,'r').xreadlines():
                champs=ligne.split("\t")
                cle=sp + "." + champs[0]
                annot[cle]=eval(champs[1])
        print >> sys.stderr," Done... "     


if len(options.annot) > 0 or len(options.reduced) > 0:
    for i in annot:
        z=z+len(annot[i])

    print >> sys.stderr,"Read ALL annotation files with ", z, "exons"


lis=lis1+lis2
print >> sys.stderr, "reading", sys.argv[1], "with W=", W, "ID=", ID, "MAS=", MAS, "X=", X, "F=", F, "T=", T, "N=", N
print >> sys.stderr, "considering the following", lsc, "species", scat

#############
# Functions #
#############

# Function to check the strand of sequence and provide coordinate on + strand in 1-based counting
def minusflip(ch):
    if ch[4]=="-":
        ch[2]=int(ch[5])-(int(ch[2]) + 1  + int(ch[3])) + 2
    return ch


# Function to extract a column from an alignment block
def getcol(b,i):
    c=[]
    for j in ESP:
        c.append(b[j][6][i:i+1])
    return c


# Function to compute some stats on a given column of the multiple alignment, modified to take into account possible mismatches (return "1" if identity criterion is ok => "pseudo-identity")
def colstat(l):
    ide=0
    cons=0
    D = {'A':0,'T':0,'G':0,'C':0,'-':0,'N':0,'ATGCN':0}

    for i in range(len(l)) :
        D[str(l[i]).upper()]+=1
	if str(l[i]) != '-' :
		D['ATGCN'] += 1

    if D['-'] < float(len(ESP))*50.0/100.0 :
	for b in ['A','T','G','C'] :
		if D[b] >= float(D['ATGCN'])*float(T)/100.0 :
			ide = 1
		if D[b] == D['ATGCN'] :
			cons = 1

    l.append(ide) # 1 if "pseudo" conserved column
    l.append(cons) # 1 if strictly conserved column
    return l

# Function to compute global stats on window of processed columns (size W)
# returns a list of stats = wstats
# stats currently include:
# nr of identical columns, id status of 1st column, id status of last column

def winstat(w):
    wid=0       # nr of identical (or 'pseudo-identical', for the relaxed version) columns in window
    wstats=[]   # list to contain all window stats
    for i in range(len(w)):
        wid = wid + w[i][len(ESP)]
    widpc=100.0*(float(wid)/float(len(w)))
    wstats.append(widpc)

    wstats.append(w[0][len(ESP)])#1st score of the window
    wstats.append(w[len(w)-1][len(ESP)])#last score of the window
    
    return wstats

# Function to filter out profiles that show low complexity (at least one base with frequency < 5%), overlap repeats or are too short

def myfilter(p):

    l = 0

    b1=p.count("A")+p.count("a")
    b2=p.count("G")+p.count("g")
    b3=p.count("C")+p.count("c")
    b4=p.count("T")+p.count("t")
    b5=p.count(".")
    #b6=p.count("~")

    if len(p) < F:
        l=1

    if b1 > float(L)*float(len(p)-b5)/100.0 or b2 > float(L)*float(len(p)-b5)/100.0 or b3 > float(L)*float(len(p)-b5)/100.0 or b4 > float(L)*float(len(p)-b5)/100.0:
    	l=2

    return l


# Function to test if a window passes the required selection criterions
# takes a list of window stats, a window size, and criterions as input
# return 1 if window is OK, 0 otherwise

def testwin(ls,ws,c1):
    if ls[0] >= c1:
        r=1
    else:
        r=0
    return r

# Function to test if a selected window overlaps annotations
# input is annotation dict, block list,start of window,end of window
# returns a 3 item list: 0|1 if window [overlaps|does not overlap], trimmed pos left, trimmed pos right, window start and end
def overlap(a,b,s,e):
    br=0
    #for j in lis:
    for j in ESP :
        if br==1:
            break 
        else:
	    if a.__contains__(b[j][1]) :
	      gaps=b[j][6].count("-",0,s)       # no. gaps up to window start
	      gape=b[j][6].count("-",0,e)       # no. gaps up to window end

	      if b[j][4]=="+":
		  ws=int(b[j][2])+int(s)-gaps       # true win start without gaps
		  we=int(b[j][2])+int(e)-gape       # true win end without gaps

	      else:
		  ws=int(b[j][2]) + (int(b[j][3]) - (int(e) - gape))
		  we=int(b[j][2]) + (int(b[j][3]) - (int(s) - gaps))

	      fl=[1,int(s),int(e),ws,we]
  #            print fl

	      # Looking for the position of the window in the list of annotations
  #           if b[j][1] in a:
	      #print >> sys.stderr, j, b[j]
	      o=bisect_left(a[b[j][1]],(ws,we))

	      if o == 0:
		  if we >  a[b[j][1]][0][0]:
		      fl[0]=0                       # it is before any annotation but
		      br=1                          # may overlap the first one
		      print >> sys.stderr, j, b[j]
		      break
		  else:
		      fl[0]=1
		      continue

	      if o == len(a[b[j][1]]):
		  if ws < a[b[j][1]][len(a[b[j][1]])-2][1]:
		      fl[0]=0
		      br=1
		      print >> sys.stderr, j, b[j]
		      break
		  else: 
		      fl[0]=1
		      continue

	      if o > 0 and o < len(a[b[j][1]]):

		  ast1=int(a[b[j][1]][o-1][0])      # annot1 start
		  aen1=int(a[b[j][1]][o-1][1])      # annot1 end
		  
		  ast2=int(a[b[j][1]][o][0])      # annot2 start
		  aen2=int(a[b[j][1]][o][1])      # annot2 end
		  

		  if (aen1-ws <= ((aen1-ast1) + (we-ws)) and aen1-ws > 0) or (aen2-ws <= ((aen2-ast2) + (we-ws)) and aen2-ws > 0):
		      if ws < aen1:
			  fl[1]=aen1-ws + s            # return to local coordinates for trimmed positions
		      if we > ast2:
			  fl[2]=e - (we-ast2) 
		      fl[0]=0
		      br=1
		      print >> sys.stderr, j, b[j]
  #                   print "and now ....", fl
		      break
		  else:
		      fl[0]=1

    return fl

# Function to walk. From a position (p) in a window (bl) the function will walk in either
# direction (d) that is, either left or right, and returns the last position where
# no more than X (given in options) consecutive mismatched columns have been seen 

def walk(bl,p,d):
    xc = 0
    cnt=0
    
    if d == "left":
        st = p-1
        en = -1
        lastmatch = p
        step=-1
        flushends=+1
        
    if d == "right":
        st = p+1
        en = len(bl)
        lastmatch = p
        step=1
        flushends=-1
#        print "walking right from ", p
    k=p

    for k in range(st,en,step):
        lastmatch = k+1

        if bl[k][len(ESP)]==0:
            xc=xc+1
        else: 
            xc=0
        cnt=cnt+1
        if xc > X:
            lastmatch=k+(flushends * xc)
            break

    if bl[lastmatch-1][len(ESP)]==0:
            lastmatch=k+(flushends * xc)
#    print "walked ", cnt, " steps"
    return lastmatch

# Function to construct the profile of identical bases in a window
## Attention, profil defini par rapport a l'homme...

def profile(w):
    p=str("")
    for i in range(len(w)):
        if w[i][len(ESP)+1]==1:
            c=str(w[i][0])
        elif w[i][len(ESP)]==1:
            #c="~"
	    c=str(w[i][0]).lower()
        else:
            c="."
        p=p+c
#    print "window size = ", len(w), "profile length = ", len(p)
    return  p

# Function to pretty print selected windows in output
# takes a block, start, end, a list of infos as input
# prints sequence coordinate positions in 1-based counting

def prettyprint(b,s,e,idx,p,f,t,t2,lab):
    gapps={}
    gappe={}

    print "%s id= %s Length= %s Trimmed= %s" % (lab,idx[0], e-s+1, t), scoreval, f, t2
    for i in ESP:
        sub=[]
        gapps[i]=b[i][6].count("-",0,s)
        gappe[i]=b[i][6].count("-",0,e)

        if b[i][4] == "+":
            pos = int(b[i][2])+s-gapps[i]+1
            end = int(b[i][2])+s-gappe[i]+1 + (e-s)
        else:
            pos = int(b[i][2]) + ( int(b[i][3]) - ( e - gappe[i]) )
            end = int(b[i][2]) + ( int(b[i][3]) - ( s - gapps[i]) )

        print "%-22s %10s %10s %-2s %s" % (b[i][1], pos, end, b[i][4], b[i][6][s:e+1]), gapps[i]
    print "profile\t\t\t\t\t\t%s" % (p)

    print ""


########
# Main #
########

#for ligne in open(sys.argv[1], 'r').xreadlines():
for ligne in gzip.open(sys.argv[1], 'r').readlines(): 
	champs=ligne.split()
	n=n+1

	# at each blank line, scan previous block if block criterion OK
	if len(champs) == 0 :
		if esp0 != '':
			if int(block[esp0][3]) >= W and len(ESP) >= N :
				#print >> sys.stderr, ESP
				currwin=[]
				scorid=[]
				i=0
				
				#collect statistics on blocks
				for j in range(int(block[esp0][3])):
					col=getcol(block,j)
					currwin.append(colstat(col))

				#collect statistics on anchor windows of size W
				for j in range(int(block[esp0][3])-W):
					lstat=winstat(currwin[j:j+W])
					scorid.append([int(lstat[0]),j,currwin[j][len(ESP)],currwin[j+W][len(ESP)]])
				
				#Compute left & right extensions as long as unused anchor windows remain
	    			while len(scorid) > 0:
					scorid.sort()#Sort windows from lowest to highest score. 
					
					last=len(scorid)-1 #Window with highest score
					
					if testwin(scorid[last],W,ID) == 1:
						
						if len(lis)>0:                  #if option -a and/or -r are set, skip & delete window if overlaps annotations
							FL=overlap(annot, block, scorid[last][1], scorid[last][1]+W-1)
			 				if FL[0]==0: #overlap an annotation
#                            					print "anchor window overlap", scorid[last]
								del scorid[last] #suppress this best-scored window
  #                           		 			print "should get smaller and smaller", scorid
			      					continue
						
						st=int(scorid[last][1])
						en=int(scorid[last][1]+W-1)
						currwin2=currwin[st:en+1]
						
						left  = walk(currwin,st,"left")
						right = walk(currwin,en,"right")
						
						t=t+1
						
						#if option -a and/or -r are set, trim edges of window to not overlap annotations
						if len(lis)>0:
							tr=""
							FL=overlap(annot, block, left, right)
							#print "checking overlap of extension", scorid[last], FL, left, right
							if FL[0]==0:
								left=FL[1]
								right=FL[2]
								tr="YES"
						else:
							tr="NO"
							FL=""
						#finst=winstat(currwin[left:right])
						finst=winstat(currwin[left:right+1])
				     		pro=profile(currwin[left:right+1])
						
						if F != "": # if -f option is set, filters out the profile
							fi=myfilter(pro)
						
						if fi == 0:
							labnr=labnr+1 # numero du CNE
							zeros=6-len(str(labnr))
							label= "SB_4PM_" + "0"*zeros + str(labnr)
							prettyprint(block,left,right,finst,pro,fi,tr, FL, label)
						else :
							#print "do not pass filter, deleting window"
							del scorid[last]
							continue
						
						# Purge block data that overlaps selected window
						for m in range(left-1,right-1):
							currwin[m][len(ESP)]=0
						
						scorid=[m for m in scorid if m[1] < left - W or m[1] > right]
						left=right=0
					
					else :
						scorid=[]
		block={}
        	continue

	# only select blocks above a minimum alignment score MAS
	#print >> sys.stderr, champs
	if champs[0]=="a":
		y=1
		cs=0
		score=champs[1].split("=")
		scoreval=int(score[1][0:len(score[1])-7])
		#scoreval=int(score[1][0:len(score[1])-2])
		ESP = []
		esp0 = ''
		if scoreval < MAS:
            		y=0

	# only include sequence lignes of species in given catalogue 
	# and counts if ALL species have been seen.

	#if y==1 and champs[0]=="s" and "_" not in champs[1]:
	if y==1 and champs[0]=="s":
        	species=champs[1].split(".")
		if species[0] in scat:
			if esp0 == '' :
				esp0 = species[0]
			champs=minusflip(champs)
			block[species[0]]=champs
			cs=cs+1
			ESP.append(species[0]) # Vector with species (in the order as they appear in the block)
			#print >> sys.stderr, esp0
