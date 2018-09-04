#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my %dist;
my %count_match;

open (F1,"<","position.csv") or die ("can't open $_, cause $!\n");
while (my $li = <F1>)
{
	chomp $li;
	if (!($li =~ m/^UTR/))
	{
		my @L = split(",",$li);
		$count_match{$L[0]}++;
		if (!exists $dist{$L[0]})
		{
			$dist{$L[0]}=[];
		}
		
		my $val = "$L[1],$L[2],$L[4],$L[5]";
		push(@{$dist{$L[0]}},$val);
	}
}
close F1;

#print Dumper (\%dist);

############################
### Distance Calculation ###
############################

my $prot1;
my $prot2;
my $motif1;
my $motif2;
my $start1;
my $start2;
my $end1;
my $end2;
my %count_dist;

open (FE2,">","tmp_file/all_distances.txt") or die "can't create file $_, cause $!\n";
foreach my $k1 (keys (%dist))
{
	my $y =1;
	for (my $x = 0 ; $x < scalar (@{$dist{$k1}}) ; $x++)
	{
		my @T = split (",", ${$dist{$k1}}[$x]);
		$prot1 = $T[0];
		$motif1 = $T[1];
		$start1 = $T[2];
		$end1 = $T[3];
		for (my $y = 0 ; $y < scalar (@{$dist{$k1}}) ; $y++)
		{
			my @R = split (",", ${$dist{$k1}}[$y]);
			$prot2 = $R[0];
			$motif2 = $R[1];
			$start2 = $R[2];
			$end2 = $R[3];
			if ($start1 != $start2)
			{
				if ($start2 > $start1)
				{
					my $add = $start2 - $end1;
					$count_dist{$k1}++;
					print FE2 "$k1,$prot1,$motif1,$start1,$end1,$prot2,$motif2,$start2,$end2,$add\n";
				}
				else
				{
					my $add = $start1 - $end2;
					$count_dist{$k1}++;
					print FE2 "$k1,$prot1,$motif1,$start1,$end1,$prot2,$motif2,$start2,$end2,$add\n";
				}
			}
		}
	}
}

close FE2;

