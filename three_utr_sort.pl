#!/usr/bin/perl

#######################################
### File creation of 3'UTR sequence ###
#######################################

use strict;
use warnings;
use Data::Dumper;


#########################
### File recuperation ###
#########################

#command line decomposition
my $commandLine = join " ", @ARGV;

############
### HELP ###
############

if (scalar (@ARGV) <=1) 
{
	print <<EOF;
	USAGE : perl three_utr_interesting.pl <Genes file> <UTR genome file> 
		Genes file could be :
			- genes symbol list
			- FlyBase genes ID
			- FlyBase transcripts ID
EOF
	exit;
}

#input and output name recuperation
my $file_gene = shift;
my $file_UTR = shift;

############################
### Gene's interest file ###
############################

my $gene;
my %interest;

open (F1,"<",$file_gene) or die "can't open $_, cause $!";
while (my $li=<F1>)
{
	chomp $li;
	my $gene = $li;
	$interest{$gene}=1;
}
close F1;


##############################
### 3'UTR sequences genome ###
##############################

my $ID_UTR;
my %UTR;
my %tirage;
my $count = 0;
my $MAX;

open (F2,"<",$file_UTR) or die "can't open $_, cause $!";
while (my $li=<F2>)
{
	chomp $li;
	if ($li=~m/^>/)
	{
		$count+=1;
		$ID_UTR = $li;
		$tirage{$count} = $ID_UTR;
	}
	else
	{
		$UTR{$ID_UTR} .= $li;
	}
	$MAX = $count;
}
close F2;

#print Dumper (\%tirage);


#####################################
### ID_gene and ID_UTR comparison ###
#####################################

my $query = 0;

open (FE,">","tmp_file/UTR_file") or die "can't create $_, cause $!";
#open (FE,">","UTR_file") or die "can't create $_, cause $!";
foreach my $ID_gene (keys (%interest))
{
	foreach my $ID_UTR (keys (%UTR))
	{
		if ($ID_UTR=~m/name=$ID_gene\S+\;\s/)
		{
			$query+=1;
			print FE "$ID_UTR\n$UTR{$ID_UTR}\n";
		}
	}
}
close FE;

#################################
### Random gene file creation ###
#################################

my $s = 1;

open (FE,">","tmp_file/random_genes") or die "can't create $_, cause $!";
#open (FE,">","random_genes") or die "can't create $_, cause $!";
while ($s <= $query)
{
	my $random = int (rand($MAX));
	if (exists $tirage{$random})
	{
		my $SEQ = $tirage{$random};
		if (exists $UTR{$SEQ})
		{
			print FE "$SEQ\n$UTR{$SEQ}\n";
		}
	}
	$s++;
}
close FE;
exit;