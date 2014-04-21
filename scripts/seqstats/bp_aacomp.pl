#!perl
use strict;
use warnings;
use Carp;

use Bio::SeqIO;
use Getopt::Long;
use Bio::SeqUtils;
use Bio::Tools::IUPAC;

my $table = new Bio::SeqUtils;
my @BASES = $table->valid_aa(0);
my %all = $table->valid_aa(2);
my ($file,$format,$help) = ( undef, 'fasta');
GetOptions(
	   'i|in:s'  => \$file,
	   'f|format:s' => \$format,
	   'h|help|?'  => \$help,
	   );

my $USAGE = "usage: aacomp [-f format] filename\n\tdefault format is fasta\n";
$file = shift unless $file;

die $USAGE if $help;

my $seqin;
if( defined $file ) {
    print "Could not open file [$file]\n$USAGE" and exit unless -e $file;
    $seqin = new Bio::SeqIO(-format => $format,
			    -file   => $file);
} else {
    $seqin = new Bio::SeqIO(-format => $format,
			    -fh     => \*STDIN);
}

my %composition;
my $total;
foreach my $base ( @BASES ) {
    $composition{$base} = 0;
}
while ( my $seq = $seqin->next_seq ) {
    if( $seq->alphabet ne 'protein' ) {
	confess("Must only provide amino acid sequences to aacomp...skipping this seq");
	next;
    }
    foreach my $base ( split(//,$seq->seq()) ) {
	$composition{uc $base}++;
	$total++;
    }
}

printf("%d aa\n",$total); 
printf("%5s %4s\n", 'aa', '#' );
my $ct = 0;
foreach my $base ( @BASES ) {
    printf(" %s %s %3d\n", $base, $all{$base}, $composition{$base} );
    $ct += $composition{$base};
}
printf( "%6s %s\n", '','-'x5);
printf( "%6s %3d\n", '',$ct);


__END__


=head1 NAME

bp_aacomp - amino acid composition of protein sequences

=head1 SYNOPSIS

  bp_aacomp [-f/--format FORMAT] [-h/--help] filename
  or
  bp_aacomp [-f/--format FORMAT] < filename
  or
  bp_aacomp [-f/--format FORMAT] -i filename

=head1 DESCRIPTION

This scripts prints out the count of amino acids over all protein
sequences from the input file.

=head1 OPTIONS

The default sequence format is fasta.

The sequence input can be provided using any of the three methods:

=over 3

=item unnamed argument

  bp_aacomp filename

=item named argument

  bp_aacomp -i filename

=item standard input

  bp_aacomp < filename

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

Based on aacomp.c from an old version of EMBOSS

=cut
