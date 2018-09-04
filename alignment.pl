#!/usr/bin/perl

#######################
### Motifs research ###
#######################

use strict;
use warnings;
use Data::Dumper;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::IUPAC;

my $commandLine = join " ", @ARGV;

############
### HELP ###
############

if (scalar (@ARGV) <=1) 
{
	print "Usage : perl alignment.pl <Motif> <UTR file> <output alignment> <output matching index>\n";
	exit;
}

######################
### Input & Output ###
######################

my $motif = shift;
my $file_UTR = shift;
my $outputfile_aln = shift;
my $outputfile_match = shift;


#####################
### Files Reading ###
#####################

my %sequence;
my %ID_sequence;
my $ID;

open (F1,"<",$file_UTR) or die "can't open $_, cause $!\n"; 
while (my $li = <F1>)
{
	chomp ($li);
	$li =~ s/\r//gs;
	if ($li =~m/^>/)
	{
		$ID = $li;
		$ID_sequence{$ID} = 1;
	}
	else 
	{
		$sequence{$ID} = $li;
	}
}
close F1;

my $ID_motif;
my %motifs;
my %distrib_motif;
my %distrib_ID_motif;

open (F1,"<",$motif) or die "can't open $_,cause $!\n"; 
while (my $li = <F1>)
{
	$li =~ s/ //gs;
	chomp $li;
	if ($li =~m/^>(\S+)/)
	{
		$ID_motif = $1;
		if (!exists ($motifs{$ID_motif}))
		{
			$motifs{$ID_motif}=[];
			$distrib_ID_motif{$ID_motif} = 0;
		}
	}
	if ($li=~m/^(\w+)/)
	{
		my $Seq = $1;
		my $seq_obj = Bio::Seq->new(-seq => "$Seq", -alphabet => 'dna');
		my $iupac = Bio::Tools::IUPAC->new(-seq => $seq_obj);
		while (my $uniqueseq = $iupac->next_seq())
		{
			
			push (@{$motifs{$ID_motif}},$uniqueseq->seq);
			$distrib_motif{$uniqueseq->seq} = 0;
		}
	}
}
close F1;

#print Dumper (\%motifs);
#print Dumper (\%distrib_motif);

#####################
### KMP execution ###
#####################

my %aln_sequence; 
my %match_sequence; 
my $seq_cut; 
my %seq_60bp;
my %distrib_stat;

foreach my $k (keys (%sequence))
{
	my $ID_seq = $k ;
	if (!exists $seq_60bp{$ID_seq})
	{
		$seq_60bp{$ID_seq} = [];
	}	
	my $fraction = 0;
	while ($fraction <length($sequence{$k}))
	{
		$seq_cut = substr($sequence{$k},$fraction,60);
		push (@{$seq_60bp{$ID_seq}},$seq_cut);
		$fraction+=60;
	}
	foreach my $c (keys (%motifs))
	{
		for (my $m = 0 ; $m<scalar(@{$motifs{$c}}); $m++)
		{
			$distrib_stat{$k}{${$motifs{$c}}[$m]} = 0;
			my $ID_match = $motifs{$c}[$m];
			my $fraction = 0;
			while ($fraction <length($sequence{$k}))
			{
				$seq_cut = substr($sequence{$k},$fraction,60);
				$fraction+=60;
				
				if (!exists $match_sequence{$seq_cut}{$ID_match})
                {
                    $match_sequence{$seq_cut}{$ID_match} = [];
                }
				
				$aln_sequence{$seq_cut}{${$motifs{$c}}[$m]} = [];
				
                my ($ref_resultats , $ref_match) = knuth_morris_pratt ($seq_cut,${$motifs{$c}}[$m]);
                my @result = @$ref_resultats;
                my @matching = @$ref_match;
                
               	#Match Preparation
                for (my $m = 0 ; $m < scalar (@matching) ; $m ++)
                {
                    push (@{$match_sequence{$seq_cut}{$ID_match}}, $matching[$m]);
                }                
                			
				# Alignment Preparation
				my @seq = split ("",$seq_cut);
				my @motif = split ("",${$motifs{$c}}[$m]);
				my @alignment;
				
				for (my $l = 0; $l <scalar(@seq) ; $l++)
				{
					$alignment[$l]="-";
				}
				
				for (my $i = 0; $i <scalar(@result) ; $i++)
				{	
					if ($result[$i] != -1)
					{
						for (my $j = 0; $j <scalar (@motif) ; $j++)
						{
							$alignment[$result[$i]+$j] =  $motif[$j];
						}
					}
				}
				my $alignment_seq = join("",@alignment);
				push (@{$aln_sequence{$seq_cut}{${$motifs{$c}}[$m]}},$alignment_seq);
			}
		}
	}
}

#print Dumper (\%match_sequence);


###################################
### Alignement View Preparation ###
###################################

my $MAX=0;

foreach my $c (keys (%motifs))
{
	my @T = split ("",$c);
	if ($MAX<scalar (@T))
	{
		$MAX = scalar(@T);
	}
}


#space in front of sequence
my @vide;
for (my $y = 0;$y<($MAX);$y++)
{
    $vide[$y]=" ";
}
my $vide = join ("",@vide);


my $ID_seq;
my $prot;
my %test_aln;

foreach my $k1 (sort keys (%aln_sequence))
{
	foreach my $k3 (sort keys (%seq_60bp))
	{
		for (my $i = 0 ; $i < scalar (@{$seq_60bp{$k3}}) ; $i ++)
		{
			if ($k1 eq ${$seq_60bp{$k3}}[$i])
			{
				$ID_seq = $k3;
				if (!exists $test_aln{$k1})
				{
					$test_aln{$k1} = [];
				}
			}
		}
	}
	foreach my $k2 (sort (keys (%{$aln_sequence{$k1}})))
	{
		foreach my $k4 (keys (%motifs))
		{
			for (my $j = 0 ; $j < scalar (@{$motifs{$k4}}) ; $j++)
			{
				if ($k2 eq ${$motifs{$k4}}[$j])
				{
					#creation entete
			        my @T = split ("",$k4);
			        while (scalar(@T)<$MAX)
			        {
			            push (@T," ");
			        }
			        $prot = join ("",@T);
				}
			}
		}
		my @aln_sort = sort(@{$aln_sequence{$k1}{$k2}});
		for (my $l = 0 ; $l < scalar (@aln_sort) ; $l++)
		{
			if ($aln_sort[$l] =~ m/[\w+]/)
			{
				push (@{$test_aln{$k1}},"$prot $aln_sort[$l]\n");
			}
		}
	}
}



my %ID_seq_60bp;

foreach my $k (keys (%seq_60bp))
{
	if (exists $ID_sequence{$k})
	{
		for (my $i = 0 ; $i < scalar (@{$seq_60bp{$k}}) ; $i++)
		{
			if (!exists $ID_seq_60bp{$seq_60bp{$k}[$i]})
			{
				$ID_seq_60bp{$seq_60bp{$k}[$i]} = [];
			}
			push(@{$ID_seq_60bp{$seq_60bp{$k}[$i]}},$k);
		}
	}
}

#print Dumper (\%ID_seq_60bp);
 


##################################
###  Matching View Preparation ###
##################################

my %match;

foreach my $k1 (keys (%sequence))
{
	my $count = 0;
	
	my @TABLE;
	my $fraction = 0;
	while ($fraction<length($sequence{$k1}))
	{
		my $seq = substr($sequence{$k1},$fraction,60);
		push (@TABLE,$seq);
		$fraction+=60;
	}
	
	for (my $i = 0 ; $i < scalar (@TABLE) ; $i ++)
	{
		if (exists $match_sequence{$TABLE[$i]})
		{
			
			foreach my $k2 (keys (%{$match_sequence{$TABLE[$i]}}))
			{
				my $keys2 = "$k2.$count";
				if (!exists $match{$k1}{$keys2})
				{
					$match{$k1}{$keys2} = [];
				}
				@{$match{$k1}{$keys2}}= @{$match_sequence{$TABLE[$i]}{$k2}} ;
			}
			$count++;
		}
	}
}
	 


my %test_match;
my $Motif ;
my $Count ;
my $value;

foreach my $k1 (keys (%match))
{
	foreach my $k2 (keys (%{$match{$k1}}))
	{
		if ($k2=~m/(\w+)\.(\d+)/)
		{
			$Motif = $1;
			$Count = $2;
			if (!exists $test_match{$k1}{$Motif})
			{
				$test_match{$k1}{$Motif} = [];			
			}
		}
		for (my $i = 0 ; $i < scalar (@{$match{$k1}{$k2}}) ; $i++)
		{
			$value = ${$match{$k1}{$k2}}[$i] + ($Count * 60);
			push (@{$test_match{$k1}{$Motif}},$value);
		}
	}
}


######################
### Alignment View ###
######################


open (FE1,">",$outputfile_aln) or die "can't create file $_, cause $!\n";
 
foreach my $k (keys (%sequence))
{
	my @TABLE;
	my $fraction = 0;
	while ($fraction<length($sequence{$k}))
	{
		my $seq = substr($sequence{$k},$fraction,60);
		push (@TABLE,$seq);
		$fraction+=60;
	}
	print FE1 "$k\n";
	for (my $i = 0 ; $i < scalar (@TABLE) ; $i ++)
	{
		if (exists $test_aln{$TABLE[$i]})
		{
			my $ALN = join ("",@{$test_aln{$TABLE[$i]}});
			print FE1 $ALN;
			print FE1 "$vide $TABLE[$i]\n";
		}

	}
}

close FE1;		 		


#####################
### Matching View ###
#####################

my $PROT;
my $ID_motif1;
my $ID_motif2;
my $ID_SEQ;
my %total_length;
my $L;
my %Count;

open (FE1,">",$outputfile_match.".csv") or die "can't create file $_, cause $!\n";

print FE1 "UTR,protein binding,motif,UTR length,start,end\n";
foreach my $k1 ( sort (keys (%test_match)))
{
	$Count{$k1} = 0;
	if ($k1=~m/length=(\d+);/)
	{
		$L = $1;
	}
	if ($k1=~m/name=(\S+);/)
	{
		$ID_SEQ = $1;	
		foreach my $k2 (keys (%{$test_match{$k1}}))
		{
			foreach my $k3 (sort (keys (%motifs)))
			{
				for (my $i = 0 ; $i < scalar (@{$motifs{$k3}}) ; $i++)
				{
					if (${$motifs{$k3}}[$i] eq $k2)
					{
						$PROT = $k3;
					}
				}
			}
			
			my @begin = Array_Unique (@{$test_match{$k1}{$k2}});
			for (my $j = 0 ; $j < scalar (@begin) ; $j++)
			{
				my $long_motif = length ($k2);
				$total_length{$k1}+=$long_motif;
				$Count{$k1}+=1;
				
				if (exists $begin[$j])
				{	
					if (exists $distrib_motif{$k2})
					{
						$distrib_motif{$k2} += 1;
						$distrib_ID_motif{$PROT} += 1;
						$distrib_stat{$k1}{$k2} += 1;
					}
					my $end = $begin[$j] + $long_motif ;
					print FE1 "$ID_SEQ,$PROT,$k2,$L,$begin[$j],$end\n";
				}
			}
		}
	}	
}
close FE1;	


##############################
### Statistic on gene data ###
##############################

open (FE1,">","tmp_file/count_file") or die ("can't create $_, cause $!\n");
foreach my $k1 (keys (%motifs))
{
	for (my $i = 0 ; $i < scalar(@{$motifs{$k1}}) ; $i++)
	{
		if (exists $distrib_motif{${$motifs{$k1}}[$i]})
		{
			print FE1 "$k1,${$motifs{$k1}}[$i],$distrib_motif{${$motifs{$k1}}[$i]}\n";
		}
	}	
}
close FE1;


################################
### Distribution Calculation ###
################################

my %motif_matrix;
my %distribution_motif;

foreach my $k2 (sort (keys(%motifs)))
{
	if (!exists $motif_matrix{$k2})
	{
		$motif_matrix{$k2} = [];
	}
	for (my $i = 0 ; $i < scalar (@{$motifs{$k2}}) ; $i ++)
	{
		for (my $row = 0 ; $row < 4 ; $row ++)
		{
			for (my $col = 0 ; $col < length (${$motifs{$k2}}[$i]) ; $col ++)
			{
				${$motif_matrix{$k2}}[$col][$row] = 0;
			}
		}
	}

	foreach my $k1 (keys (%distrib_motif))
	{
		for (my $i = 0 ; $i < scalar (@{$motifs{$k2}}) ; $i ++)
		{
			if (${$motifs{$k2}}[$i] eq $k1)
			{
				my $distribution;
				my @motif = split("",$k1);
				for (my $j = 0 ; $j < scalar(@motif) ; $j ++)
				{
					$distribution = ($distrib_motif{$k1} / $distrib_ID_motif{$k2});
					$distribution_motif{$k1} = $distribution;
					if ($motif[$j] eq 'A')
					{
						${$motif_matrix{$k2}}[$j][0] += (1 * $distribution) ;
					}
					if ($motif[$j] eq 'C')
					{
						${$motif_matrix{$k2}}[$j][1] += (1 * $distribution);
					}
					if ($motif[$j] eq 'G')
					{
						${$motif_matrix{$k2}}[$j][2] += (1 * $distribution);
					}
					if ($motif[$j] eq 'T')
					{
						${$motif_matrix{$k2}}[$j][3] += (1 * $distribution);
					}
				}
				my $pourcent = $distribution * 100;
				#print "$k2\t$k1\t$pourcent\n";
			}
		}
	}
}




#################
### Print PWM ###
#################

foreach my $k (keys (%motif_matrix))
{
	open (FE4,">","pwm_".$k.".txt") or die "can't create file $_, cause $! \n"; 
	my @table = @{$motif_matrix{$k}};

	for (my $row = 0 ; $row < 4 ; $row ++)
	{
		for (my $col = 0 ; $col < scalar (@table) ; $col ++)
		{
			print FE4 "$table[$col][$row]\t";
		}
		print FE4 "\n";
	}
	close FE4;
}



###############################
### Random motif generation ###
###############################

my %nucleotid;
$nucleotid{0} = "A";
$nucleotid{1} = "C";
$nucleotid{2} = "G";
$nucleotid{3} = "T";


my %Random_motif;

foreach my $k1 (keys (%motifs))
{
	my $L_motif = 0;
	my $nb_motif= scalar (@{$motifs{$k1}});
	
	for (my $i = 0 ; $i < scalar (@{$motifs{$k1}}) ; $i ++)
	{
		$L_motif = length (${$motifs{$k1}}[$i]); 
		#print $k1."\t".$L_motif ."\n";
		my $m =1;
		my @rand_motif;
		while ($m <= $L_motif)
		{
			my $nb_nt = int (rand (3));
			#print "$k1\t$nb_nt\n";
			if (exists ($nucleotid{$nb_nt}))
			{
				my $nt = $nucleotid{$nb_nt};
				push (@rand_motif,$nt);
			}
			$m++;
		}
		my $motif_R = join ("",@rand_motif);
		$Random_motif{${$motifs{$k1}}[$i]} = $motif_R;
	}
}	


my %distrib_random;

foreach my $k1 (keys (%sequence))
{
	foreach my $k2 (keys (%Random_motif))
	{
		$distrib_random{$k1}{$k2} = 0;
		my ($random_result , $random_match) = knuth_morris_pratt ($sequence{$k1},$Random_motif{$k2});
		my @rand_match = @$random_match;
		$distrib_random{$k1}{$k2} += scalar (@rand_match); 
	}
}




#######################
### Statistic Table ###
#######################

my $id;

foreach my $k1 (keys (%distrib_random))
{
	if ($k1=~m/name=(\S+);/)
	{
		$id = $1;
	}
	foreach my $k2 (keys(%motifs))
	{
		if (scalar(@{$motifs{$k2}}) != 1)
		{
			open (FE2,">","tmp_stat/$id\_$k2") or die "can't create file $_, cause $!\n";
			#open (FE2,">","$id\_$k2") or die "can't create file $_, cause $!\n";
			print FE2 "motif,match,random\n";
			for (my $i = 0 ; $i < scalar (@{$motifs{$k2}}) ; $i++)
			{
				my $m = $motifs{$k2}[$i];
				if (defined $distribution_motif{$m})
				{
					my $stat = $distrib_stat{$k1}{$m} * $distribution_motif{$m};
					my $random = $distrib_random{$k1}{$m} * $distribution_motif{$m};
					print FE2 "$m,$stat,$random\n";
				}
			}
		}
	}
}

close FE2;

#######################################
###  Knuth-Morris-Pratt's Algorithm ###
#######################################

#computation of the prefix subroutine
sub knuth_morris_pratt_next
{
   my($P) = @_; #pattern
   use integer;
   my ( $m, $i, $j ) = ( length $P, 0, -1 );
   my @next;
   for ($next[0] = -1; $i < $m; ) 
   {
      # Note that this while() is skipped during the first for() pass.
      while ( $j > -1 && substr( $P, $i, 1 ) ne substr( $P, $j, 1 ) ) 
      {
         $j = $next[$j];
      }
      $i++;
      $j++;
      $next[$i] = substr( $P, $j, 1 ) eq substr( $P, $i, 1 ) ? $next[$j] : $j;
   }
   return ( $m, @next ); # Length of pattern and prefix function.
}

#matcher subroutine
sub knuth_morris_pratt
{
    my ( $T, $P , $M) = @_; # Text and pattern.
    use integer;
    my ($m,@next) = knuth_morris_pratt_next( $P );
    my ( $n, $i, $j ) = ( length($T), 0, 0 );
    
    my @match;
    
    #my @next;
    my @val;
    my $k=0;
    while ( $i < $n ) 
    {
        while ( $j > -1 && substr( $P, $j, 1 ) ne substr( $T, $i, 1 ) ) 
        {
            $j = $next[$j];
        }
        $i++;
        $j++;
        if($j>=$m)
        {
            $val[$k]= $i - $j; # Match
            push (@match,$val[$k]);
        }
        else
        {
            $val[$k]=-1; # Mismatch
        }
        $k++;
    }
    return ( \@val , \@match );
}



####################################
### Removing Duplicates in Table ###
####################################

sub Array_Unique
{
    my @List = @_;
    my %FutureList;
    foreach(@List)
    {
        $FutureList{$_} = 1; # supprime les doublons
    }
    return (keys(%FutureList));
}