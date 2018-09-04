#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;


my %query;
my $cle1;

open (F1,"<","tmp_file/count_file") or die ("can't open $_, cause $!\n");
while (my $li = <F1>)
{
	chomp $li;
	my @L = split (",",$li);
	$cle1 = "$L[0],$L[1]";
	$query{$cle1}=$L[2];

}
close F1;


my %rand;
my $cle2;

open (F1,"<","tmp_file/rand_file") or die ("can't open $_, cause $!\n");
while (my $li = <F1>)
{
	chomp $li;
	my @L = split (",",$li);
	$cle2 = "$L[0],$L[1]";
	$rand{$cle2}=$L[2];
	
}
close F1;

my %rbp;

foreach my $k (keys (%query))
{
	if (exists ($rand{$k}))
	{
		my @T = split (",",$k);
		my $RBP = $T[0];
		if (!exists $rbp{$RBP})
		{
			$rbp{$RBP} = [];
		}
		my $val = "$T[1],$query{$k},$rand{$k}";
		push (@{$rbp{$RBP}},$val);
	}
}

foreach my $k1 (keys (%rbp))
{
	if (scalar(@{$rbp{$k1}}) != 1)
	{
		open (FE1,">>","tmp_file/$k1\_count") or die "can't create $_, cause $! \n";
		for (my $i = 0 ; $i < scalar (@{$rbp{$k1}}) ; $i++)
		{
			print FE1 "${$rbp{$k1}}[$i]\n";
		}
		close FE1;
	}
}