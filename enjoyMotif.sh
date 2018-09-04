#!/bin/bash

#$1 = distance minimale fixée
#$2 = liste de gènes
#$3 = UTR génome
#$4 = fasta de motifs
#$5 = sortie alignement
#$6 = sortie position
#$7 = sortie distance entre motif


#export PATH = $HOME/bin

chemin=`dirname $0`

function print_usage()
{
	echo "Usage :"
	echo "$0 -d [DISTANCE] -g [GENELIST] -G [UTRGenome] -m [MOTIF] -R [T/F] -r [mirna.fa] -o [OUTPUT]"
	echo " "
	echo "$0 finds motifs in interest gene list and calculate distance between different motifs"
	echo " "
	echo "Default values :"
	echo "Distance : 20"
	echo "R : F"
	echo "OUTPUTDIR : result_utr"
}

if [ "$1" == "" ]
then
	print_usage
	exit
fi

while getopts d:g:G:m:R:r:o:h option
do
	case $option in
		d) DIST=$OPTARG ;;
		g) GENELIST=$OPTARG ;;
		G) UTRGenome=$OPTARG ;;
		m) MOTIF=$OPTARG ;;
		R) MIR=$OPTARG ;;
		r) MIRFILE=$OPTARG ;;
		o) OUTPUTDIR=$OPTARG ;;
		?) print_usage ;;
	esac
done

if [ -d $OUTPUTDIR ]
then
	echo "Please remove or rename the directory '$OUTPUTDIR'."
	exit
else
	mkdir $OUTPUTDIR
fi

# file for 3'UTR visualization
cp $chemin/dessin.js $OUTPUTDIR/visualization.js
cp $chemin/dessin.html $OUTPUTDIR/visualization.html

cd $OUTPUTDIR

mkdir $OUTPUTDIR/tmp_file
mkdir $OUTPUTDIR/tmp_stat

grep -e '$' $chemin/PAS.txt | cat - $MOTIF > $OUTPUTDIR/tmp_file/motifs.txt

perl $chemin/three_utr_sort.pl $GENELIST $UTRGenome

# miRNA alignment option
if [ "$MIR" = "T" ]
then
# alignment on UTR's query genes
	mkdir $OUTPUTDIR/tmp_file/index
	bowtie2-build -f $OUTPUTDIR/tmp_file/UTR_file $OUTPUTDIR/tmp_file/index/UTR
	bowtie2 -a -x $OUTPUTDIR/tmp_file/index/UTR -f -U $MIRFILE > $OUTPUTDIR/alignment.sam
# alignment on UTR's random genes
	mkdir $OUTPUTDIR/tmp_file/index_rand
	bowtie2-build -f $OUTPUTDIR/tmp_file/random_genes $OUTPUTDIR/tmp_file/index_rand/UTR_rand
	bowtie2 -a -x $OUTPUTDIR/tmp_file/index_rand/UTR_rand -f -U $MIRFILE > $OUTPUTDIR/tmp_file/alignment_rand.sam
	perl $chemin/sam_sort.pl
	Rscript $chemin/mirna_stat.r $OUTPUTDIR/tmp_file/stat_mir.txt $OUTPUTDIR/mirna_stat
else
	touch $OUTPUTDIR/mirna_alignment.txt
fi

perl $chemin/random_gene.pl $OUTPUTDIR/tmp_file/random_genes $OUTPUTDIR/tmp_file/motifs.txt

if [ -e tmp_file/UTR_file ]
then
    perl $chemin/alignment.pl $OUTPUTDIR/tmp_file/motifs.txt $OUTPUTDIR/tmp_file/UTR_file $OUTPUTDIR/aln $OUTPUTDIR/position
else
    echo "Can't find UTR_file create by three_utr.pl"
    exit
fi

if [ -e position.csv ]
then
	perl $chemin/distance_calculation.pl
	perl $chemin/distance_treat.pl -d$DIST
else
	echo "Can't find position file"
	exit
fi

for file in pwm_*
do
	Rscript $chemin/seqlogo.r $file $file.png
done

for file in tmp_stat/*
do
	Rscript $chemin/aln_statistic.r $file $file.stat
	cat $file.stat >> $OUTPUTDIR/tmp_stat/cat_stat.txt
done

if [ -e tmp_stat/cat_stat.txt ]
then
	perl $chemin/stat_creation.pl
fi

perl $chemin/gene_compare.pl

for file in $OUTPUTDIR/tmp_file/*_count
do 
	Rscript $chemin/gene_stat.r $file $file.stat
	cat $file.stat >> $OUTPUTDIR/tmp_file/gene_stat.stat
done

if [ -e tmp_file/gene_stat.stat ]
then
	perl $chemin/stat_gene_creation.pl 
fi

if [ -e tmp_file/Statistic_motif.csv ] && [ -e tmp_file/Statistic_gene.csv ]
then
	perl $chemin/final_stat.pl
fi

rm -Rf tmp_stat
rm -Rf tmp_file