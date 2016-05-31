# phylogenomic_pipeline
## This repository contains a few steps for building a phylogenomic dataset from aTRAM output


## 1. Exon-stitching pipeline
###Will extract and stitch Exons together from aTRAM contigs using the program Exonerate
This code will take aTRAM outputs for many genes for many taxa (e.g. 'best.fasta' files) use a reference gene and the program exonerate to annotate exons fromt the aTRAM contigs and stich the exons together.  The final output will be one file for each gene that includes all the taxa with the exons annotated. Missing sections are filled in with NNN so the alignment step is simpler. 

All the scripts need to be in the same folder as the aTRAM contig files. Usually we use the .best.fasta files from the aTRAM output.  

dependencies:  Exonerate  http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate
	script was written with Exonerate v. 2.2

Edit the Loop_Exonerate_pipeline.sh  file:

need:  common file output name from aTRAM.  e.g.   Gene.library.atram.best.out    would use <atram.best.ou>
need:  path to reference genes
need:  List of library names

usage:  sh Loop_Exonerate_pipeline.sh

## 2.  Tree building 
## code used to build large phylogenomic trees with ASTRAL 
