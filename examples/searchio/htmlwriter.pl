#!/usr/bin/perl 

# Demonstrates the use of a SearchIO Blast parser and a SearchWriterI object
# for producing HTML Blast output from a Blast report input stream.
#
# Usage:
#   STDIN:  none; supply filename of BLAST report on command-line
#   STDOUT: none; generates an output file "searchio.html"
#           containing HTML-formatted Blast Report
#   STDERR: Any errors that occurred.
#
# For more documentation about the writer, including
# a complete list of columns, see the docs for 
#   Bio::SearchIO::Writer::HTMLResultWriter.
#
# For more documentation about working with Blast result objects,
# see docs for these modules:
#   Bio::Search::Result::BlastResult
#   Bio::Search::Iteration::IterationI
#   Bio::Search::Hit::BlastHit
#   Bio::Search::HSP::BlastHSP
#
# For more documentation about the Blast parser, see docs for
#   Bio::SearchIO
#
# Author: Steve Chervitz <sac@bioperl.org>


use strict;
use lib '../../';

use Bio::SearchIO;
use Bio::SearchIO::Writer::HTMLResultWriter;

my $outfile = "searchio.html";
my $file = shift or die "Usage: $0 <BLAST-report-file>\n       HTML output is saved to $outfile\n";

my $in = Bio::SearchIO->new( -format => 'blast', 
                             -file => $file,  #comment this out to read STDIN
                             #-fh => \*ARGV,  #uncomment this to read from STDIN
                             -verbose => 0 );

my $writer = new Bio::SearchIO::Writer::HTMLResultWriter();
my $out = new Bio::SearchIO(-writer => $writer,
                            -file   => ">$outfile");


while ( my $result = $in->next_result() ) {
    eval {
        # printf STDERR "Report %d: $result\n", $in->result_count;
        $out->write_result($result, 1);
    };
    if($@) {
        warn "Warning: Blast parsing or writing exception caught for $result:\n$@\n";
    }
}

printf STDERR "\n%d Blast report(s) processed.\n", $in->result_count;
printf STDERR "Output sent to file: %s\n",  $out->file if $out->file;
