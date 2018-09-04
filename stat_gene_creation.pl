#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;


my $ID_stat;
my %p_val;

open (F1,"<","tmp_file/gene_stat.stat") or die "can't open $_, cause $!\n";
while (my $li = <F1>)
{
	if ($li=~m/tmp_file\/(\S+)\_count/)
	{
		$ID_stat = $1;
	}
	if ($li=~m/p-value\s=\s(\S+)/)
	{
		$p_val{$ID_stat} = $1;
	}
}
close F1;


open (FE1,">","tmp_file/Statistic_gene.csv") or die "can't open $_, cause $!\n";
print FE1 "Protein Binding,p-value,Significant\n";
foreach my $ID (keys (%p_val))
{
	if ($p_val{$ID}=~m/[\d+\.]/)
	{
		if ($p_val{$ID} < 0.1)
		{
			print FE1 "$ID,$p_val{$ID},*\n";
		}
		else 
		{
			print FE1 "$ID,$p_val{$ID},NS\n";
		}
	}
	else
	{
		print FE1 "$ID,$p_val{$ID},NA\n";
	}
}
print FE1 "Poly(A)Signal,NA,*\n";
close FE1;