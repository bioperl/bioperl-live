#!/usr/bin/perl

=head1 NAME

bp_mutate.pl - randomly mutagenize a single protein or DNA sequence

=head1 SYNOPSIS

  ./bp_mutate.pl -p 25 -i test.fa -n 5 -f swiss -o muts.swiss

  #or

  ./bp_mutate.pl --percent=25 --input=test.fa --number=5 -output=x.fa

=head1 DESCRIPTION

Randomly mutagenize a single protein or DNA sequence one or more times.
Specify percentage mutated and number of resulting mutant sequences.
Print mutagenized sequences to STDOUT or write to an output file.

  -h|--help    Help
  -p|--percent Percent mutagenized
  -n|--number  Number of mutant sequences created
  -o|--output  Output file (optional)
  -f|--format  Output format (default: fasta)
  -i|--input   Input file

=head1 FEEDBACK

User feedback is an integral part of the evolution of this and other
Bioperl scripts. Send your comments and suggestions to the Bioperl 
mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Brian Osborne, bosborne at alum.mit.edu

=cut

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;

my ( $help, $in_file, $percent, $out_file, $number );
my $format = "fasta";
my @dna    = qw(a g c t);
my @amino  = qw(a c d e f g h i k l m n p q r s t v w y);

GetOptions(
    "h|help"      => \$help,
    "p|percent:i" => \$percent,
    "n|number:i"  => \$number,
    "o|output:s"  => \$out_file,
    "f|format:s"  => \$format,
    "i|input:s"   => \$in_file
);

usage() if ( $help || !$percent || !$in_file || !$number || $percent > 100 );

# Seed the random number generator. "time|$$" combines the
# current time with the current process id
srand( time | $$ );

my $seqio   = Bio::SeqIO->new( -file => $in_file );
my $seqobj  = $seqio->next_seq;
my $num_mut = percent_to_num($percent);
my @seq_arr = ();

# don't keep a mutant that's already been made
while ( $number > $#seq_arr + 1 ) {
    my $mut_seq = mutate_all( $seqobj, $num_mut );
    push @seq_arr, $mut_seq unless ( grep /$mut_seq/, @seq_arr );
}

foreach my $mut_seq (@seq_arr) {
    my $name   = $seqobj->display_id . "-${percent}_percent-$number";
    my $outseq = Bio::Seq->new(
        -seq        => $mut_seq,
        -display_id => $name,
        -desc       => $seqobj->desc
    );
    my %args = ( -format => $format );
    $args{-file} = ">>$out_file" if $out_file;
    my $seqio = Bio::SeqIO->new(%args);
    $seqio->write_seq($outseq);
    $number--;
}

# mutagenize the sequence, one-by-one
sub mutate_all {
    my ( $seq_obj, $num ) = @_;
    my $type = $seq_obj->alphabet;
    my $str  = $seq_obj->seq;

    # store the mutagenized positions in $positions
    my $positions = "";
    for ( my $i = 0 ; $i < $num_mut ; ++$i ) {
        ( $str, $positions ) = mutate_one( $str, $type, $positions );
    }
    $str;
}

# mutagenize one position
sub mutate_one {
    my ( $str, $type, $positions ) = @_;
    my ( $position, $new_char );

    # pick a random position in the sequence, checking
    # that the position isn't already mutagenized
    do {
        $position = random_position($str);
    } until ( $positions !~ /\b$position\b/ );
    $positions .= "$position ";
    my $current_char = substr( $str, $position, 1 );

    # pick a random char that's not the existing char
    do {
        $new_char = random_char($type);
    } until ( $new_char ne $current_char );
    substr( $str, $position, 1, $new_char );
    ( $str, $positions );
}

# randomly select a position in the sequence
sub random_position {
    my $string = shift;
    int( rand( length($string) ) );
}

# randomly select one of the chars depending on alphabet
sub random_char {
    my $type = shift;
    $type eq "protein"
      ? return $amino[ rand @amino ]
      : return $dna[ rand @dna ];
}

sub percent_to_num {
    my $percent = shift;
    int( $percent / 100 * length( $seqobj->seq ) );
}

sub usage {
    exec( 'perldoc', $0 );
    exit(0);
}

__END__
