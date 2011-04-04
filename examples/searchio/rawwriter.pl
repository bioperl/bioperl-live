#!/usr/bin/env perl 

# Demonstrates the use of a SearchIO Blast parser for producing
# output of raw HSP data from a Blast report input stream.
#
# Shows how to print out raw BLAST alignment data for each HSP.
#
# Usage:
#   STDIN:  stream containing one or more BLAST or PSI-BLAST reports.
#   STDOUT: Raw alignment data for each HSP of each hit (BLAST format)
#   STDERR: Progress info and any errors.
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
#
# TODO: 
#   * Implement a Bio::SearchIO::Writer::HSPTextWriter object
#     that can do this. Then this example can fit into the standard
#     model used by the other writer examples in which a writer
#     object is created and hooked up with a SearchIO output object.

use strict;

use lib '../../';

use Bio::SearchIO;

# In this case, we only want raw alignments, and we only need to screen
# on significance info (E- or P-value) so we don't need
# to do a full parse of the alignments. Thus, we're using a -shalow_parse
# flag to indicate that we don't need to parse alignments. This should
# result in faster processing.
# TODO: Convert this to use -format='blast'. Shallow-parse option not supported there.
my $in = Bio::SearchIO->new(-format => 'psiblast', 
                            -fh => \*ARGV,
			    -signif => 0.1,
			    -shallow_parse => 1,
			    -hold_raw_data => 1 );

while ( my $result = $in->next_result() ) {
  print STDERR "\nBLAST Results for $result\n\n";
  my $count = 0;
  foreach( $result->hits ) {
    print "Alignment for hit #", ++$count, "\n\n";
    print $_->raw_hit_data();
  }
  print "=" x 50 , "\n";
}

printf STDERR "\n%d Blast report(s) processed.\n", $in->result_count;

