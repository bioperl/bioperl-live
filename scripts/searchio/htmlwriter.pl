#!/usr/bin/perl 

# Demonstrates the use of a SearchIO Blast parser and a SearchWriterI object
# for producing HTML Blast output from a Blast report input stream.
#
# Usage:
#   STDIN:  stream containing one or more BLAST or PSI-BLAST reports.
#   STDOUT: none, but generates an output file "searchio.html"
#           containing HTML-formatted Blast Report(s)
#   STDERR: Any errors that occurred.
#
# For more documentation about the writer, including
# a complete list of columns, see the docs for 
#   Bio::SearchIO::Writer::HTMLResultWriter.
#
# For more documentation about working with Blast result objects,
# see docs for these modules:
#   Bio::Search::Result::BlastResult
#   Bio::Search::Hit::BlastHit
#   Bio::Search::HSP::BlastHSP
#
# For more documentation about the PSI-Blast parser, see docs for
#   Bio::SearchIO::psiblast
#
# Author: Steve Chervitz <sac@bioperl.org>
# Revision: $Id$


use strict;
use lib '../../';

use Bio::SearchIO;
use Bio::SearchIO::Writer::HTMLResultWriter;

print "\nUsing SearchIO->new()\n";

my $in = Bio::SearchIO->new( -format => 'psiblast', 
                             -signif => 0.1, 
                             -verbose => 0 );

my $writer = new Bio::SearchIO::Writer::HTMLResultWriter();
my $out = new Bio::SearchIO(-writer => $writer,
                            -file   => ">searchio.html");


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
