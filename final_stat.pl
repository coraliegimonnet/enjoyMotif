#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my %motif;

open (F1,"<","tmp_file/Statistic_motif.csv") or die "can't open $_, cause : $! \n";
while (my $li = <F1>)
{
	chomp $li;
	if (!($li =~ m /^Transcript/))
	{
		$motif{$li}=1;
	}
}
close F1;

my %gene;

open (F1,"<","tmp_file/Statistic_gene.csv") or die "can't open $_, cause : $! \n";
while (my $li = <F1>)
{
	chomp $li;
	if (!($li=~m/^Protein/))
	{
		my @L = split (",",$li);
		$gene{$L[0]}="$L[1],$L[2]";
	}
}
close F1;


open (FE1,">","statistic.csv") or die "can't create $_, cause : $! \n";
print FE1 "Transcript ID,RBP or PAS,quantitative p-value,significant,qualitative p-value,significant\n";
foreach my $k1 (keys (%motif))
{
	foreach my $k2 (keys (%gene))
	{
		my @T = split (",",$k1);
		if ($T[1] eq $k2)
		{
			print FE1 "$k1,$gene{$k2}\n";
		}
	}
}
close FE1;
