#!perl

use strict;
use warnings;

use Bio::SeqIO;
use Getopt::Long;
use GD;

use vars qw( $USAGE %VALIDFORMATS);

%VALIDFORMATS = ( 'png'  => 1,
		  'jpeg' => 1,
		  'gd2'  => 1,
		  'gd'   => 1,
		  'gif'	 => 1,
		  'wbmp' => 1 );

$USAGE = "usage:\tchaos_plot -i/--input=INPUTFILE -f/--format=SEQFORMAT \n".
    "\t-o/--output=OUTPUTFILE -g/--graphics=GRAPHIC TYPE\n".
    "\t-w/--width=600 -h/--height=400\n";

$USAGE .= "\tValid graphics formats: (" . join(",", ( keys %VALIDFORMATS )) .")\n";
$USAGE .= "\tImage size defaults to 600x400, SEQFORMAT to fasta\n";
$USAGE .= "\tINPUTFILE can also be read from STDIN\n";

my ($format,$graph,$width,$height,$seqfile,$output) = ('fasta', 'png', 600, 400);
GetOptions( "i|input:s"           => \$seqfile,
	    "f|format:s"          => \$format,
	    "o|output:s"          => \$output,
	    "g|graph|graphics:s"  => \$graph,
            "width:i"             => \$width,
            "height:i"            => \$height
	    );

if( ! $output || ! $VALIDFORMATS{$graph} ) {
    die $USAGE ;
}
my $seqin;
$seqfile = shift unless $seqfile;
if( defined $seqfile ) {
    print "Could not open file [$seqfile]\n$USAGE" and exit unless -e $seqfile;
    $seqin = new Bio::SeqIO(-format => $format,
			    -file   => $seqfile);
} else {
    $seqin = new Bio::SeqIO(-format => $format,
			    -fh     => \*STDIN);
}

my $img = new GD::Image($width,$height);
my $white = $img->colorAllocate(255,255,255); 
my $black = $img->colorAllocate(0,0,0); 

my $seq = $seqin->next_seq;
die("Sequence type must be DNA not " . $seq->alphabet())
    unless $seq->alphabet ne 'dna' or $seq->alphabet ne 'rna';
my %nmerdata;
my $len = $seq->length();
my $max = 0;

my ($x,$y) = ( 0.5, 0.5);
$img->string(gdGiantFont, 1,1, 'A', $black);
$img->string(gdGiantFont, 0,$height - 15, 'C', $black);
$img->string(gdGiantFont, $width - 15,1, 'T', $black);
$img->string(gdGiantFont, $width - 15,$height -20, 'G', $black);

for( my $i = 1; $i <= $len; $i++ ) {

    my $base = lc $seq->subseq($i,$i);
    if( $base eq 'a' ) {
	$x *= 0.5;
	$y  *= 0.5;
    } elsif ( $base eq 'g' ) {
	$x = ( $x + 1.0 ) * 0.5;
	$y  = ( $y + 1.0  ) * 0.5; 
    } elsif ( $base eq 'c' ) {
	$x *= 0.5;
	$y  = ( $y + 1.0  ) * 0.5; 
    } elsif ( $base eq 't' or $base eq 'u' ) {
	$x = ( $x + 1.0 ) * 0.5;
	$y  *= 0.5;
    }

    $img->setPixel($x * $width,$y * $height, $black);
}
open my $OUT, '>', $output or die "Could not write file '$output': $!\n";
binmode $OUT;
$graph =~ s/jpg/jpeg/;

print $OUT $img->$graph();
close $OUT;


__END__

=head1 NAME

bp_chaos_plot - a chaos plot from DNA and RNA sequences

=head1 SYNOPSIS

  bp_chaos_plot.pl -i/--input=INPUTFILE -f/--format=SEQFORMAT
        -o/--output=OUTPUTFILE -g/--graphics=GRAPHIC FORMAT
        -w/--width=WIGHT -h/--height=HEIGHT

=head1 DESCRIPTION

This scripts generates image files using GD image library to visualize
nucleotide sequences using chaos plot.

=head1 OPTIONS

Valid graphics formats are currently gd, gd2, png, wbmp, jpeg and gif.

The default size of the image file is 600x400.

The sequence input can be provided using any of the three methods:

=over 3

=item unnamed argument

  bp_chaos_plot filename

=item named argument

  bp_chaos_plot -i filename

=item standard input

  bp_chaos_plot < filename

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

This code is based on EMBOSS C code for chaos.c by Ian Longden.
Included are documentation from EMBOSS code:

Chaos produces a chaos plot.  The original application is part of the
ACEDB genome database package, written by ** Richard Durbin (MRC LMB,
UK) rd@mrc-lmba.cam.ac.uk, and Jean Thierry-Mieg (CRBM du CNRS,
France) mieg@crbm1.cnusc.fr

=cut
