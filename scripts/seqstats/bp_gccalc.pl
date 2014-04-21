#!perl

use strict;
use warnings;

use Bio::SeqIO;
use Bio::Tools::SeqStats;
use Getopt::Long;
my $format = 'fasta';
my $file;
my $help =0;
my $aggregate;
GetOptions(
	    'f|format:s' => \$format,
	    'i|in:s'     => \$file,
	    'h|help|?'    => \$help,
	    'a|aggregate!'=> \$aggregate,
	    );


my $USAGE = "usage: gccalc.pl -f format -i filename\n";
if( $help ) {
    die $USAGE;
}

$file = shift unless $file;
my $seqin;
if( defined $file ) {
    print "Could not open file [$file]\n$USAGE" and exit unless -e $file;
    $seqin = new Bio::SeqIO(-format => $format,
			    -file   => $file);
} else {
    $seqin = new Bio::SeqIO(-format => $format,
			    -fh     => \*STDIN);
}
my ($total_base, $total_gc);
while( my $seq = $seqin->next_seq ) {
    next if( $seq->length == 0 );
    if( $seq->alphabet eq 'protein' ) {
	warn("gccalc does not work on amino acid sequences ...skipping this seq");
	next;
    }

    my $seq_stats  =  Bio::Tools::SeqStats->new('-seq'=>$seq);
    my $hash_ref = $seq_stats->count_monomers();  # for DNA sequence
    print "Seq: ", $seq->display_id, " ";
    print $seq->desc if $seq->desc;
    print " Len:", $seq->length, "\n";
    $total_base += $seq->length;
    $total_gc   += $hash_ref->{'G'} + $hash_ref->{'C'};
    printf "GC content is %.4f\n", ($hash_ref->{'G'} + $hash_ref->{'C'}) /
	$seq->length();

    foreach my $base (sort keys %{$hash_ref}) {
	print "Number of bases of type ", $base, "= ", $hash_ref->{$base},"\n";
    }
    print "--\n";
}
if( $aggregate ) {
  printf "Total GC content is %.4f out of %d bases\n", $total_gc / $total_base, $total_base;
}
# alternatively one could use code submitted by
# cckim@stanford.edu

sub calcgc {
    my $seq = $_[0];
    my @seqarray = split('',$seq);
    my $count = 0;
    foreach my $base (@seqarray) {
	$count++ if $base =~ /[G|C]/i;
    }

    my $len = $#seqarray+1;
    return $count / $len;
}


__END__

=head1 NAME

bp_gccalc - GC content of nucleotide sequences

=head1 SYNOPSIS

  bp_gccalc [-f/--format FORMAT] [-h/--help] filename
  or
  bp_gccalc [-f/--format FORMAT] < filename
  or
  bp_gccalc [-f/--format FORMAT] -i filename

=head1 DESCRIPTION

This scripts prints out the GC content for every nucleotide sequence
from the input file.

=head1 OPTIONS

The default sequence format is fasta.

The sequence input can be provided using any of the three methods:

=over 3

=item unnamed argument

  bp_gccalc filename

=item named argument

  bp_gccalc -i filename

=item standard input

  bp_gccalc < filename

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 HISTORY

Based on script code (see bottom) submitted by cckim@stanford.edu

Submitted as part of bioperl script project 2001/08/06

=cut
