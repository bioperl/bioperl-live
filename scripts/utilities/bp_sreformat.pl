#!perl
# Author:  Jason Stajich <jason-at-bioperl-dot-org>
# Purpose: Bioperl implementation of Sean Eddy's sreformat
#          We're not as clever as Sean's squid library though so
#          you have to specify the input format rather than letting
#          the application guess.

use strict;
use warnings;
use Bio::SeqIO;
use Bio::AlignIO;
use Getopt::Long;

my $USAGE = "bp_sreformat -if INFORMAT -of OUTFORMAT -i FILENAME -o output.FORMAT

-h/--help               Print this help
-if/--informat          Specify the input format
-of/--outformat         Specify the output format
-i/--input              Specify the input file name
                        (to pass in data on STDIN use minus sign as filename)
-o/--output             Specify the output file name
                        (to pass data out on STDOUT use minus sign as filename)
--msa                   Specify this is multiple sequence alignment data
--special=specialparams Specify special params supported by some formats
                        Comma or space separated please.
                        These include:
                        nointerleaved   -- for phylip,non-interleaved format
                        idlinebreak     -- for phylip, makes it molphy format
                        percentages     -- for clustalw, show % id per line
                        flat            -- don't show start-end in seqid
                        linelength      -- line length for clustalw
                        mrbayes         -- for MrBayes proper NEXUS output
";


my ($input,$output,$informat,$outformat,$msa,$special);

GetOptions(
	   'h|help'          => sub { print STDERR ($USAGE); exit(0) },
	   'i|input:s'         => \$input,
	   'o|output:s'        => \$output,
	   'if|informat:s'     => \$informat,
	   'of|outformat:s'    => \$outformat,
	   'msa'               => \$msa,
	   's|special:s'       => \$special,
	   );

unless( defined $informat && defined $outformat ) { 
    die(sprintf("Cannot proceed without a defined input and output you gave (%s,%s)\n",
		defined $informat ? $informat : "''" ,
		defined $outformat ? $outformat : "''"));
}

my ($in,$out);
my @extra;
if( $special ) {
    @extra = map { my @rc;
		   if( /nointerleaved/) {
		       @rc = ('-interleaved' => '0');
		   } elsif( /mrbayes/ ) {
		       @rc = ('-show_symbols' => 0,
			      '-show_endblock' => 0);
		   } elsif( /(\S+)\=(\S+)/ ) { @rc = ( "-$1" => $2) } 
	           else{ @rc = ("-$_" => 1) }
		   @rc;
	       } split(/[\s,]/,$special);
}
# guess we're talking about MSA if any of the standard MSA names are used
if( $informat =~ /nexus|phylip|clustal|maf|stockholm|bl2seq|msf/ ||
    $outformat =~ /nexus|phylip|clustal|maf|stockholm|bl2seq|msf/ ) {
    $msa = 1;
}

if( $msa ) {
    eval {
	if( defined $input ) {
	    $in = new Bio::AlignIO(-format => $informat, -file => $input);
	} else {
	    $in = new Bio::AlignIO(-format => $informat, -fh => \*ARGV);
	}
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
    while( my $aln = $in->next_aln) { 
	 if( $special =~ /flat/ ) {$aln->set_displayname_flat(1); }
	 $out->write_aln($aln) }

} else {
    eval {
	if( defined $input ) {
	    $in = new Bio::SeqIO(-format => $informat, -file => $input);
	} else { 
	    $in = new Bio::SeqIO(-format => $informat, -fh => \*ARGV);
	}
    };
    if( $@ ) {
	if( $@ =~ /Could not open/ ) {
	    die("Could not open input file: $input\n");
	} else { 
	    die("Unknown sequence format to bioperl $informat\n");
	}	
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
	if( $@ =~ /Could not open/ ) {
	    die("Could not open output file: $output\n");
	} else { 
	    die("Unknown sequence format to bioperl $outformat: $@\n");
	}
    }
    while( my $seq = $in->next_seq ) {
	$out->write_seq($seq);
    }
}

=head1 NAME

bpsreformat - convert sequence formats

=head1 DESCRIPTION

This script uses the SeqIO system that allows conversion of sequence
formats either sequence data or multiple sequence alignment data.  The
name comes from the fact that Sean Eddy's program sreformat (part of
the HMMER pkg) already does this.  Sean's program tries to guess the
input formats while in our code we currently require your to specify what
the input and output formats are and if the data is from a multiple
sequence alignment or from straight sequence files.

Usage:

bpsreformat -if INFORMAT -of OUTFORMAT -i FILENAME -o output.FORMAT

  -h/--help        Print this help

  -if/--informat   Specify the input format

  -of/--outformat  Specify the output format

  -i/--input       Specify the input file name
                   (to pass in data on STDIN use minus sign as filename)
  -o/--output      Specify the output file name
                   (to pass data out on STDOUT use minus sign as filename)

  --msa            Specify this is multiple sequence alignment data

  --special        Will pass on special parameters to the AlignIO/SeqIO
                   object -- most of these are for Bio::AlignIO objects
                   Comma separated list of the following
                   nointerleaved   -- for phylip,non-interleaved format
                   idlinebreak     -- for phylip, makes it molphy format
                   percentages     -- for clustalw, show % id per line

=cut
