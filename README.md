# enjoyMotif

Date : February 2017 
Last version : June 2017

**This tool based on Knuth-Morris-Pratt Algorithm searches motifs like binding sites**
**of RNA-binding protein on Drosophila 3'UTR. It searches too polyadenylation signal.**
**It can search microRNA binding site using Bowtie2.**	 

Installation required:

- Perl5
- BioPerl
- R (version 3.3.3 or higher)
- R package SeqLogo :
	source("https://bioconductor.org/biocLite.R")
	biocLite("seqLogo")
- Bowtie2

*Install this tool on PATH: export PATH="/path/to/tool:$PATH"*

Usage:

	enjoymotif.sh -d [DISTANCE] -g [GENELIST] -G [UTRGenome] -m [MOTIF] -R [T/F] -r [mirna.fa] -o [OUTPUT]

Options:

	-d: maximum distance between 2 binding sites
	-g: gene list
	-G: fasta file of 3'UTR genome sequences from FlyBase
	-m: fasta file of motifs
	-R: T or F: use Bowtie2 to align miRNA on 3'UTR
	-r: fasta file of miRNA sequences
  	-o: output directory

Author :
	Coralie Gimonnet (https://github.com/coraliegimonnet)
