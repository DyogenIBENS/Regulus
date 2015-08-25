#!/usr/bin/python

import sys
from optparse import OptionParser



annot={}

n=0
N=0

##############
# Parameters #
##############

usage = "usage: %prog [redundant UCSC annotation file]"

##########
# Module #
##########

# Module to read options
parser = OptionParser(usage=usage)

parser.add_option("-w", dest="W", default="0", help="Minimum size of anchor window to initiate a search, as specified in the scanmaf.py script. If two EST positions are separated by less than this distance, then they will be merged into one interval. Useful to increase the compaction of the annotation file, since no alignment can be selected in a region less than the window size. [default = 0 bp]")

(options, args) = parser.parse_args()

W=int(options.W)



########
# MAIN #
########


fle= sys.argv[1]
print >> sys.stderr, "Reading annotation file", fle

for ligne in open(fle,'r').xreadlines():
    champs=ligne.split()
    exstart=champs[8].split(",")
    exend=champs[9].split(",")
    cle=champs[1]
    for i in range(int(champs[7])):
        if (cle) not in annot:
            annot[cle]=[(int(exstart[i]),int(exend[i]))]
        else:
            annot[cle].append((int(exstart[i]),int(exend[i])))

for i in annot:

    n=0
    N=(len(annot[i])/10)+1

    print >> sys.stderr,"Removing redundancy in", i, "with ", len(annot[i]) , "exons"

    print >> sys.stderr,"Sorting...",i,    
    annot[i].sort(reverse=True)
    print >> sys.stderr,"Done"
    temp2=[]
    (p1,p2)=annot[i].pop()
    print >> sys.stderr,"Removing redundancy in...",i, "...",  len(annot[i])    
    while len(annot[i])>0:
        (p3,p4)=annot[i].pop()
        n=n+1
        if p2+W >= p3:
            p2=max(p2,p4)
        else:
            temp2.append((p1,p2))
            p1=p3
            p2=p4
        if n % N == 0:
            print >> sys.stderr, (10*n)/N, "%...",
    print >> sys.stderr, "Done"
    temp2.append((p1,p2))
    temp2.sort()
    annot[i]=temp2

    print "%s\t%s" % (i,  annot[i])
    
    print >> sys.stderr,"Removed redundancy with", len(annot[i]), "=", int(100*len(annot[i])/(10*N)), " % of initial intervals left"
        
