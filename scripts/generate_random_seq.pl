#!/usr/local/bin/perl -w
use strict;
use vars qw($USAGE);

use Bio::SeqIO;
use Bio::PrimarySeq;
use Getopt::Long;
my ($length,$type,$filename);

$USAGE = 'usage: generate_random_seq.pl --length=1000 --type=dna --filename=/tmp/test.seq';

my %alphabets = ( 'dna' => [qw(C A G T)],
		  'rna' => [qw(C A G U)],
		  'prot'=> [qw( A C D E F G H I K L M N P Q R S T V W Y)],
	      );
&GetOptions 
    (
     'l|length:s'         => \$length,
     't|type|m|alphabet:s' => \$type,
     'f|file|filename:s'  => \$filename
     );

assert ( $type && defined ($alphabets{lc $type}),
	 $USAGE);
assert ( $length && $length =~ /^\d+$/, $USAGE );

my $sequence = "";
my $alphabet = $alphabets{lc $type};
my $sspace = scalar @$alphabet;
foreach ( 1..$length ) {
    $sequence .= $alphabet->[ int rand($sspace) ];
}

my $seq = new Bio::PrimarySeq(-seq => $sequence, -display_id => 'randomseq');
my %args = ( '-format' => 'fasta');
if( $filename ) { $args{'-file'} = ">$filename" }

my $seqio =  new Bio::SeqIO(%args);
$seqio->write_seq($seq);
sub assert { die $_[1] unless( $_[0] ); }
