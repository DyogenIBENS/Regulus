This README.txt file describes three scripts that together form the REGULUS package:

### scanmaf.py ###
Scans a maf multiple alignement to search regions of particular conservation.
 
Usage example:
scanmaf.py example_hg19_chr22.maf.gz -w 10 -x 3 -i 90 -f 10 -t 88 -l 80 -n 3 -c hg19,gorGor3,otoGar3 > example.res
 
Options: 
- w : Minimum size of anchor window to initiate a search [default = 50 bp]
- x : Max cutoff for consecutive column mismatches [default = 3 columns]
- i : Minimum %ID of alignment in anchor window [default = 90 %]
- s : Minimum score of a block in the MAF file before it is considered for a search [default = 300]
- f : Filters out output profiles that contain low complexity sequences and/or that fall below the length provided as argument and/or that overlap small cap nucleotides (repeat masked) in one of the species given in the -c option. Low complexity means that at least one base has a frequency below 5%. [default = 0 bp]
- l : Threshold for low complexity sequences filtering; sequences containing more thant L% of the same nucleotide are discarded [default = 80%]
- t : Minimum similarity of a column to be considered as conserved [default = 60 %]
- n : Minimum number of species required in the alignment block to be considered [default = 6]
- c : Catalogue of species (at least 2) to consider when searching the MAF file for the regions of interest. Argument is a coma separated list of species, with species names as indicated in the MAF file. [compulsory ; default = None]
- a : If this option is set, a directory called SM_annot must be in the current directory with annotation files in UCSC format (>10 tab separated columns) called species_annot.ucsc, where species is replaced by the species names as indicated in the MAF file. Annotation coordinates are O-based, in UCSC fashion. The argument is a coma separated list of species for which annotations files are provided in the SM_annot directory. When this option is set, no output will be produced in regions of the multiple alignements of the MAF file that overlap the annotations in these genomes. Generally these annotations are non redundant (e.g. ensembl genes) and will be used as provided. This option can be combined with -r for other species.  [default = none].
- r : If this option is set, a directory called SM_annot must be in the current directory with annotation files in rr format. This format is produced by the remred.py script, which removes redundancy in large redundant annotation files such as mammalian EST alignment files. Annotation coordinates are O-based, in UCSC fashion. The argument is a coma separated list of species for which annotations files are provided in the SM_annot directory. When this option is set, no output will be produced in regions of the multiple alignements of the MAF file that overlap the annotations in these genomes. This option can be combined with -a for other species. In general, using this option instead of -a reduces computation time [default = none]


### genes_env.py ###
Lists all the genes that are neighbour of each CNE in Human, along with their 'status' in the other species possessing the considered CNE.
 
Usage example:
 genes_env.py chr22 fam_chr22.list 

Needs several input files:
- fam_chr22.list : list of elements we search targets for (example provided in Examples/). Each line corresponds to one element, with its number and its name in the different species where it is present.
- ancGenes.Boreoeutheria.list.bz2 and ancGenes.Coelomata.list.bz2 : lists of homologous gene families dating back to Boreoeutheria or Coelomata ancestors (retrieved from Ensembl Compara: http://www.ensembl.org/info/genome/compara/index.html)
- genesST.%species.list.bz2 : list of species genes from Ensembl.
- 46_species2.list : list of the different species analyzed (see Examples/).
- 46_species_sizes.list : file listing the parameters used for each species to calculate gene scores (weighting factors and maximum CNE-gene distance allowed, functions of genome size and assembly quality). See Examples/.
- %species_chr22.bed.gz : bedfiles with the coordinates of the CNEs in the different species.


### score_genes_targets.py ###
Calculate the association score of each gene found in the environment of a family of orthologous CNEs.

Usage example :
score_genes_targets.py genes_env_chr22.list.gz
with "genes_env_chr22.list.gz" the output of genes_env.py.

Needs several input files:
- species_synteny.list : a list of th weighting factors depending of the degree of rearrangement of each species genome compared to the human one (see Examples/).
- ancGenes.Boreoeutheria.list.bz2 : lists of homologous gene families dating back to Boreoeutheria ancestors (retrieved from Ensembl Compara: http://www.ensembl.org/info/genome/compara/index.html)
- names.Homo.sapiens.list.bz2 : human gene names.
- genesST.Homo.sapiens.list.bz2 : list of human genes from Ensembl.



