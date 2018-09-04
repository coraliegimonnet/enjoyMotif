#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

open (F1,"<","tmp_file/random_genes") or die "can't open $_, cause $!\n";
while (my $li = <F1>)
{
	chomp $li;
}
close F1;