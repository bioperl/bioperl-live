#!/usr/bin/env perl 

# Example usage of a Bio::SearchIO::psiblast parser 
# of traditional format Blast reports.
# Shows how to print out raw BLAST alignment data for each hit.
#
# Usage:
#   STDIN:  stream containing one or more BLAST or PSI-BLAST reports.
#   STDOUT: Raw alignment data for each HSP of each hit (BLAST format)
#   STDERR: Progress info.
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

use lib '../../../';

use Bio::SearchIO;

my $in = Bio::SearchIO->new(-format => 'psiblast', 
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

printf STDERR "\n%d Blast report(s) processed.\n", $in->report_count;

