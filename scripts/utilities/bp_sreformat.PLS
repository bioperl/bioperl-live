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
--msa                   Specify this is multiple sequence alignment data";

my ($input,$output,$informat,$outformat,$msa);

GetOptions(
	   'h|help:s'          => sub { die($USAGE) },
	   'i|input:s'         => \$input,
	   'o|output:s'        => \$output,
	   'if|informat:s'     => \$informat,
	   'of|outformat:s'    => \$outformat,
	   'msa'               => \$msa,
	   );

unless( defined $informat && defined $outformat ) { die("Cannot proceed without a defined input and output you gave ($informat,$outformat)\n") }

my ($in,$out);
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
				       -file => ">$output");
	   } else { 
	       # default to STDOUT for output
	       $out = new Bio::AlignIO(-format => $outformat);
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

=cut
