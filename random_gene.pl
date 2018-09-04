#/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::IUPAC;

my $commandLine = join " ", @ARGV;

my $random_file = shift;
my $input_motif = shift;


my %random;
my $ID_rand;

open (F1,"<",$random_file) or die "can't open $_, cause $!";
while (my $li=<F1>)
{
	chomp ($li);
	if ($li=~m/^>/)
	{
		$ID_rand = $li;
	}
	else
	{
		$random{$ID_rand} = $li;
	}
}
close F1;

my %motifs;
my $ID_motif;

open (F1,"<",$input_motif) or die "can't open $_, cause $!";
while (my $li = <F1>)
{
  $li =~ s/ //gs;
	chomp $li;
	if ($li=~m/^>(\S+)/)
	{
		$ID_motif = $1;
	}
	if ($li=~m/^(\w+)/)
	{
		my $Seq = $1;
		my $seq_obj = Bio::Seq->new(-seq => "$Seq", -alphabet => 'dna');
		my $iupac = Bio::Tools::IUPAC->new(-seq => $seq_obj);
		while (my $uniqueseq = $iupac->next_seq())
		{
			push (@{$motifs{$ID_motif}},$uniqueseq->seq);
		}
	}
}
close F1;


my %rand;
foreach my $k1 (keys (%random))
{
	foreach my $k2 (keys (%motifs))
	{
		for (my $i = 0 ; $i < scalar (@{$motifs{$k2}}) ; $i++)
		{ 
			my ($rand_results , $rand_match) = knuth_morris_pratt ($random{$k1},${$motifs{$k2}}[$i]);
			my @rand_match = @$rand_match;
      $rand{${$motifs{$k2}}[$i]}+=scalar(@rand_match);  
		}
	}
}


#####################
### File creation ###
#####################

open (FE1,">","tmp_file/rand_file") or die ("can't create $_, cause $!\n");
foreach my $k1 (keys (%motifs))
{
  for (my $i = 0 ; $i < scalar(@{$motifs{$k1}}) ; $i++)
  {
    if (exists $rand{${$motifs{$k1}}[$i]})
    {
      print FE1 "$k1,${$motifs{$k1}}[$i],$rand{${$motifs{$k1}}[$i]}\n";
    }
  }
}
close FE1;


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