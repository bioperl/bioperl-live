#!perl
# Author:  Jason Stajich <jason@bioperl.org>
# Purpose: Bioperl implementation of Sean Eddy's sreformat
#          We're not as clever as Sean's squid library though so
#          you have to specify the input format rather than letting
#          the application guess.

# TODO - support STDIN/STDOUT piping?
#      - finish POD
#      - stress test

use strict;
use Bio::SeqIO;
use Bio::AlignIO;
use Getopt::Long;

my $USAGE = "bpsreformat -if INFORMAT -of OUTFORMAT -i FILENAME -o output.FORMAT

-h/--help               Print this help
-if/--informat          Specify the input format
-of/--outformat         Specify the output format
-i/--input              Specify the input file name
                        (to pass in data on STDIN omit this flag)
-o/--output             Specify the output file name
                        (to pass data out on STDOUT omit this flag
--msa                   Specify this is multiple sequence alignment data
--special=specialparams Specify special params supported by some formats
                        Comma or space separated please.
                        These include:
                        nointerleaved   -- for phylip,non-interleaved format
                        idlinebreak       -- for phylip, makes it molphy format
                        percentages     -- for clustalw, show % id per line

";


my ($input,$output,$informat,$outformat,$msa,$special);

GetOptions(
	   'h|help:s'          => sub { die($USAGE) },
	   'i|input:s'         => \$input,
	   'o|output:s'        => \$output,
	   'if|informat:s'     => \$informat,
	   'of|outformat:s'    => \$outformat,
	   'msa'               => \$msa,
	   's|special:s'       => \$special,
	   );

unless( defined $informat && defined $outformat )
    { die("Cannot proceed without a defined input and output you gave ($informat,$outformat)\n") }

my ($in,$out);
my @extra;
if( $special ) {
    @extra = map { my @rc;
		   if( /nointerleaved/) {
		       @rc = ('-interleaved' => '0');
		   } else{ @rc = ("-$_" => 1) }
		   @rc;
	       } split(/[\s,]/,$special);
}
if( $msa ) {
    eval {
	$in = new Bio::AlignIO(-format => $informat, -file => $input);
    };
    if( $@ ) {
	die("Unknown MSA format to bioperl $informat\n");
    }
    eval {
	   if( $output ) {
	       $out = new Bio::AlignIO(-format => $outformat,
				       -file => ">$output", @extra);
	   } else {
	       # default to STDOUT for output
	       $out = new Bio::AlignIO(-format => $outformat,@extra);
	   }
       };
    if( $@ ) {
	die("Unknown MSA format to bioperl $outformat\n");
    }
    while( my $aln = $in->next_aln) { $out->write_aln($aln) }

} else {
    eval {
	$in = new Bio::SeqIO(-format => $informat, -file => $input);
    };
    if( $@ ) {
	die("Unknown sequence format to bioperl $informat\n");
    }
    eval {
	   if( $output ) {
	       $out = new Bio::SeqIO(-format => $outformat,
				       -file => ">$output");
	   } else {
	       # default to STDOUT for output
	       $out = new Bio::SeqIO(-format => $outformat);
	   }
       };
    if( $@ ) {
	die("Unknown sequence format to bioperl $outformat\n");
    }
    while( my $seq = $in->next_seq ) {
	$out->write_seq($seq);
    }
}

=head1 NAME

bpsreformat - convert sequence formats

=head1 DESCRIPTION

This module is an implementation that allows conversion of sequence
formats either sequence data or multiple sequence alignment data.  The
name comes from the fact that Sean Eddy's program sreformat (part of
the HMMER pkg) already does this.  Sean's program tries to guess the
input formats in our code we currently require your to specify what
the input and output formats are and if the data is from a multiple
sequence alignment or from straight sequence files.

Usage:

bpsreformat -if INFORMAT -of OUTFORMAT -i FILENAME -o output.FORMAT

  -h/--help        Print this help

  -if/--informat   Specify the input format

  -of/--outformat  Specify the output format

  -i/--input       Specify the input file name
                   (to pass in data on STDIN omit this flag)

  -o/--output      Specify the output file name
                   (to pass data out on STDOUT omit this flag

  --msa            Specify this is multiple sequence alignment data

  --special        Will pass on special parameters to the AlignIO/SeqIO
                   object -- most of these are for Bio::AlignIO objects
                   Comma separated list of the following
                   nointerleaved   -- for phylip,non-interleaved format
                   idlinebreak     -- for phylip, makes it molphy format
                   percentages     -- for clustalw, show % id per line

=cut
