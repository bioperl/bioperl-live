#!/bin/perl
use strict;
use vars qw($USAGE);

# random sequence generator #
# -c=1 option will cause prot sequences to be built 
# using vertebrate aa frequencies, 
# with option -a putting a 1st methionine residues on. Frequencies are
# calculated from the NCBI human RefSeq protein sequences 
# -c and -a only affect protein sequences
# -a only works in conjunction with -c
# -n number of random sequences, default = 1

use Bio::PrimarySeq;
use Bio::SeqIO;
use Getopt::Long;
my ($length,$type,$filename,$comp,$met);

$USAGE = 'usage: generate_random_seq.pl --length=1000 --type=dna --filename=/tmp/test.seq --number=50';

my %alphabets = ( 'dna' => [qw(C A G T)],
                  'rna' => [qw(C A G U)],
                  'prot'=> [qw( A C D E F G H I K L M N P Q R S T V W Y)],
              );
# make random num from 1-10000. numbers in this array reflect the frequency,
# e.g., a random number from 1.744 = A, 745-991 = C etc;
my @aa_frequencies = qw(744 991 1398 2017 2378 3104 3349 3726 4239 5273 5443 
                        5749 6410 6848 7455 8263 8760 9340 9488 9713 10000);
my $number = 1;

&GetOptions
  (
   'l|length:s'          => \$length,
   't|type|m|alphabet:s' => \$type,
   'f|file|filename:s'   => \$filename,
   'c|composition:s'     => \$comp,
   'a|methionine:s'      => \$met,
   'n|number:s'          => \$number
  );

assert ( $type && defined ($alphabets{lc $type}),
         $USAGE);
assert ( $length && $length =~ /^\d+$/, $USAGE );

foreach my $num (1..$number) {
   my $sequence = "";
   my $alphabet = $alphabets{lc $type};
   my $sspace = scalar @$alphabet;
   if (!$comp || $type ne 'prot') {
      foreach ( 1..$length ) {
	 $sequence .= $alphabet->[ int rand($sspace) ];
      }
   }elsif ($type eq 'prot') {
      $sequence = build_seq($length, \@aa_frequencies);
   }
   my $seq =  Bio::PrimarySeq->new(-seq        => $sequence, 
				   -display_id => 'randomseq'.$num);
   my %args = (-format => 'fasta');
   if( $filename ) { $args{-file} = ">>$filename" }
   my $seqio = Bio::SeqIO->new(%args);
   $seqio->write_seq($seq);
}

sub assert { die $_[1] unless( $_[0] ); }

sub build_seq {
   #takes seqlen and ref to frequency data as parameters
   my ($len, $pf)  = @_;
   my $str = ($met)?'M':'';
   my $i = ($met)?1:0;
   for ($i..$len-1) {
      my $aa = int(rand (10000)) ;
      my $j = 0;
      while ($pf->[$j] < $aa && $j <19) {
	 $j++;
      }
      $str .= $alphabets{'prot'}[$j];
   }
   print "str is $str\n";
   return $str;
}
