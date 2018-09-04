#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

#####################
### Query's Genes ###
#####################

my %gene;
my $count = 0;

open (F1,"<","tmp_file/UTR_file") or die ("can't open UTR file, cause $!\n");
while (my $li = <F1>)
{
	chomp $li;
	if ($li=~m/^>/)
	{
		$count++;
		my $id;
		if ($li=~m/(FBtr\d+)/)
		{
			$id = $1;
		}
		if ($li=~m/name=(\S+);/)
		{
			$gene{$id} = $1;
		}
	}
}

my %start;
my %end;
my %seq;
my %aln;
my %mir_q;

open (F1,"<","alignment.sam") or die ("can't open SamFile, cause $!\n");
while (my $li = <F1>)
{
	chomp $li;
	if (!($li=~m/^@/))
	{
		my @T = split ("\t",$li);
		if ($T[2] ne "*")
		{
			my $mir = "$T[0]_$T[2]";
			$mir_q{$T[0]}++;
			$aln{$mir} = $1;
			$start{$mir} = $T[3];
			$seq{$mir} = $T[9];
			$end{$mir} = $T[3] + length($T[9]);
		}
	}
}
close F1;

######################
### Random's Genes ###
######################

my %gene_r;

open (F1,"<","tmp_file/random_genes") or die ("can't open UTR file, cause $!\n");
while (my $li = <F1>)
{
	chomp $li;
	if ($li=~m/^>/)
	{
		my $id;
		if ($li=~m/(FBtr\d+)/)
		{
			$id = $1;
		}
		if ($li=~m/name=(\S+);/)
		{
			$gene_r{$id} = $1;
		}
	}
}

my %start_r;
my %end_r;
my %seq_r;
my %aln_r;
my %mir_r;

open (F1,"<","tmp_file/alignment_rand.sam") or die ("can't open SamFile, cause $!\n");
while (my $li = <F1>)
{
	chomp $li;
	if (!($li=~m/^@/))
	{
		my @T = split ("\t",$li);
		if ($T[2] ne "*")
		{
			my $mir = "$T[0]_$T[2]";
			$mir_r{$T[0]}++;
			$aln_r{$mir} = $1;
			$start_r{$mir} = $T[3];
			$seq_r{$mir} = $T[9];
			$end_r{$mir} = $T[3] + length($T[9]);
		}
	}
}
close F1;


####################
### miRNAs stats ###
####################

open (FE1,">","tmp_file/stat_mir.txt") or die "can't create statFile, cause $!\n";
print FE1 "UTR,query,random\n";
foreach my $k (keys (%mir_q))
{
	if (defined $mir_r{$k})
	{
		print FE1 "$k,$mir_q{$k},$mir_r{$k}\n";
	}
	else
	{
		print FE1 "$k,$mir_q{$k},0\n";
	}
}
foreach my $k (keys (%mir_r))
{
	if (!defined $mir_q{$k})
	{
		print FE1 "$k,0,$mir_r{$k}\n";
	}
}
close FE1;



###################
### miRNAs file ###
###################

my $mir;
my $FBtr;
my %fbtr;
open (FE1,">","mirna_alignment.txt") or die ("can't create $_, cause : $!\n");
foreach my $k1 (keys (%aln))
{
	if ($k1=~m/(\S+)\_(\S+)/)
	{
		$mir = $1;
		$FBtr = $2;
	}
	if (exists $gene{$FBtr})
	{
		print FE1 "$gene{$FBtr},";
		$fbtr{$gene{$FBtr}}=1;
	}
	if (defined $start{$k1} && $end{$k1} && $seq{$k1})
	{
		print FE1 "$mir,$start{$k1},$end{$k1},$seq{$k1}\n";
	}
}
close FE1;


#########################
### miRNAs proportion ###
#########################

my $Count;

foreach my $k (keys (%fbtr))
{
	$Count++;
}

my $proportion = ($Count/$count) * 100;
print "The proportion of genes potentially targeted by miRNAs is: $proportion\n ";