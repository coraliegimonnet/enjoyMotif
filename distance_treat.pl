#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;


my $commandLine = join " ", @ARGV;

############
### HELP ###
############

if (scalar (@ARGV) == 0) 
{
	print "Usage : perl distance_treat.pl -d<int>\n";
	exit;
}

############
### ARGV ###
############

my $distance_users = shift;
my $distance;
if ($distance_users =~m/\-d(\d+)/)
{
	$distance = $1;
}

my %count_match;

open (F1,"<","position.csv") or die "can't open $_, cause : $!";
while (my $li = <F1>)
{
	chomp $li;
	my @T = split (",",$li);
	$count_match{$T[0]} += 1;
}
close F1;

my %final_dist;
my %count_dist;

open (F1,"<","tmp_file/all_distances.txt") or die "can't open $_, cause : $!";
while (my $li = <F1>)
{
	chomp $li;
	my @T = split (",",$li);
	if ($T[1] ne $T[5])
	{
		if ($T[3] < $T[7])
		{
			if ($T[9] <= $distance)
			{
				$final_dist{$li} = 1;
				$count_dist{$T[0]}+=1;
			}
		}
	}
}
close F1;

###############################
### Distance table creation ###
###############################

open (FE1,">","distance.csv") or die ("can't create file $_, cause $!\n");
print FE1 "UTR,protein 1,motif 1,start 1,end 1,protein 2,motif 2,start 2,end 2,distance\n";
foreach my $k (keys (%final_dist))
{
	print FE1 "$k\n";
}
close FE1;

#########################
### Score calculation ###
#########################

open (FE1,">","ratio_distance.csv") or die ("can't create file $_, cause $! \n");
print FE1"UTR,Ratio\n";
foreach my $k1 (keys (%count_match))
{
	if ($count_match{$k1} != 0 && exists $count_dist{$k1})
	{
		my $c = $count_dist{$k1} ;
		my $score = ($c / $count_match{$k1}) * 100 ;
		my $sc = sprintf("%0.2f",$score);
		print FE1 "$k1,$sc\%\n";
	}
	else
	{
		print FE1 "$k1,NA\n";
	}
}

close FE1;