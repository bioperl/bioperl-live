#!/usr/local/bin/perl 

# Example usage of a SearchIO::psiblast parser of traditional format Blast 
# and PSI-Blast reports.
# Illustrates how to grab a set of SeqFeatures from a Blast report.
# This parser represents a new and improved version of Bio/Tools/Blast.pm.
#
# Usage:
#   STDIN:  stream containing one or more BLAST or PSI-BLAST reports.
#   STDOUT: feature start, end data
#   STDERR: Processing info, such as the number of reports processed
#           and the number of hitless reports.
# 
# For more documentation about working with Blast result objects,
# see to documentation for these modules:
#   Bio::Search::Result::BlastResult
#   Bio::Search::Hit::BlastHit
#   Bio::Search::HSP::BlastHSP
#
# For more documentation about the PSI-Blast parser, see docs for
#   Bio::SearchIO::psiblast
#
# Author: Steve Chervitz <sac@bioperl.org>

use strict;
use lib '../../';
use Bio::SearchIO;

my $in = Bio::SearchIO->new( -format => 'psiblast',
                             -fh => \*ARGV,
			     -signif => 0.1, 
			     -verbose => 0 );
my @hitless_reports = ();

while ( my $blast = $in->next_result() ) {

    if( $blast->hits ) {
      while( my $feature = $blast->next_feature() ) {
	print "Feature from ", $feature->start, " to ", $feature->end, "\n";
      }
    }
    else {
      push @hitless_reports, $blast;
    }
}

printf STDERR "\n%d Blast report(s) processed.\n", $in->result_count;
printf STDERR "\n%d reports had no hits:\n", scalar(@hitless_reports);

foreach my $blast (@hitless_reports) {
    print STDERR "No hits for query ", $blast->query_name;
    print STDERR ($blast->no_hits_found ? "\n" : "(filtered)\n")
;
}
  

