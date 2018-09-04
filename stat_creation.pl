#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;


my $ID_stat;
my %p_val;

open (F1,"<","tmp_stat/cat_stat.txt") or die "can't open $_, cause $!\n";
while (my $li = <F1>)
{
	if ($li=~m/tmp_stat\/(\S+)\"/)
	{
		$ID_stat = $1;
	}
	if ($li=~m/p-value\s=\s(\S+)/)
	{
		$p_val{$ID_stat} = $1;
	}
}
close F1;


my $SEQ;
my $PROT;
my %hash;


foreach my $ID (keys (%p_val))
{
	if ($ID=~m/(\S+)\_(\S+)/)
	{
		$SEQ = $1;
		$PROT = $2;
	}
	if (!exists $hash{$SEQ})
	{
		$hash{$SEQ}=[];
	}
	if ($p_val{$ID}=~m/[\d+\.]/)
	{
		if ($p_val{$ID} < 0.001)
		{
			push(@{$hash{$SEQ}},"$PROT,$p_val{$ID},***");
		}
		elsif ($p_val{$ID} < 0.01)
		{
			push(@{$hash{$SEQ}},"$PROT,$p_val{$ID},**");
		}
		elsif ($p_val{$ID} < 0.05)
		{
			push(@{$hash{$SEQ}},"$PROT,$p_val{$ID},*");
		}
		else 
		{
			push(@{$hash{$SEQ}},"$PROT,$p_val{$ID},NS");
		}
	}
	else
	{
		push(@{$hash{$SEQ}},"$PROT,$p_val{$ID},NA");
	}
}


open (FE1,">","tmp_file/Statistic_motif.csv") or die "can't open $_, cause $!\n";
print FE1 "Transcript ID,Protein Binding,p-value,Significant\n";
foreach my $k1 (keys(%hash))
{
	for (my $i = 0 ; $i < scalar(@{$hash{$k1}}) ; $i++)
	{
		print FE1 "$k1,${$hash{$k1}}[$i]\n";
		print FE1 "$k1,Poly(A)Signal,NA,*\n";
	}
}
close FE1;
