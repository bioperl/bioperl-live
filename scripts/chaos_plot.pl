#!/usr/local/bin/perl -w

# $Id$

# This script provides a chaos plot of DNA
# Cared for by Jason Stajich <jason@cgt.mc.duke.edu>
# Copyright Jason Stajich
# You may distribute this module under the same terms as perl itself

# This code is based on EMBOSS C code for chaos.c by Ian Longden
# included are documentation from EMBOSS code:
#
# Chaos produces a chaos plot.
# The original application is part of the ACEDB genome database
# package, written by ** Richard Durbin (MRC LMB, UK)
# rd@mrc-lmba.cam.ac.uk, and Jean Thierry-Mieg (CRBM du CNRS,
# France) mieg@crbm1.cnusc.fr

use strict;

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

$USAGE = "usage:\tchaos_plot.pl -i/--i=INPUTFILE -f/--format=SEQFORMAT \n".
    "\t-o/--output=OUTPUTFILE -g/--graphics=GRAPHIC TYPE\n";

$USAGE .= "\tvalid graphics formats: (" . join(",", ( keys %VALIDFORMATS )) .")\n";

my ($format,$graph ,$seqfile,$output) = ('fasta', 'png');
GetOptions( "i|input:s"           => \$seqfile,
	    "f|format:s"          => \$format,
	    "o|output:s"          => \$output,
	    "g|graph|graphics:s"  => \$graph,
	    );

if( ! $seqfile || ! -e $seqfile || ! $output ) {
    die( $USAGE );
}
my $seqin = new Bio::SeqIO(-format => $format,
			   -file   => $seqfile);

my ($width, $height) = (600,400);

my $img = new GD::Image($width,$height);
my $white = $img->colorAllocate(255,255,255); 
my $black = $img->colorAllocate(0,0,0); 

my $seq = $seqin->next_seq;
die("Sequence type must be DNA not " . $seq->moltype())
    unless $seq->moltype ne 'dna' or $seq->moltype ne 'rna';
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
open(OUT, ">$output");
binmode OUT;
$graph =~ s/jpg/jpeg/;

print OUT $img->$graph();
