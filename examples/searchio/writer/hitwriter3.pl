#!/usr/local/bin/perl 

# Example usage of a SearchIO::psiblast parser of traditional format Blast reports.
# Same as hitwriter.pl but showing how to use FileHandles instead of objects.
#
# Note that since we're working with file handles,
# we can't call $in->report_count as in hitwriter.pl
#
# Usage:
#   STDIN:  stream containing one or more BLAST or PSI-BLAST reports.
#   STDOUT: none, but generates an output file "hitwriter3.out"
#           containing tab-delimited data on a per-hit basis.
#   STDERR: Progress info.
#
# For more documentation about the writer, including
# a complete list of columns, execute:
#   perldoc Bio::SearchIO::Writer::HitTableWriter.
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
# Author: Steve Chervitz <steve_chervitz@affymetrix.com>
# Revision: $Id$


use strict;

use lib '../../../';

use Bio::SearchIO;
use Bio::SearchIO::Writer::HitTableWriter;

print STDERR "\nUsing BlastIO->newFh()\n";

my $in     = Bio::SearchIO->newFh( -format => 'psiblast',
				   -signif => 0.1 );
my $writer = Bio::SearchIO::Writer::HitTableWriter->new();
my $out    = Bio::SearchIO->newFh( -format => 'psiblast',
				   -writer => $writer,
				   -file   => ">hitwriter3.out" );
my $hit_count = 0;
while ( my $blast = <$in> ) {
  # Can't do this since $in isn't a Bio::SearchIO object reference:
  #printf STDERR "Report %d: $blast\n", $in->report_count;
  printf STDERR "\nReport #%d: $blast\n", ++$hit_count;
  
  if( $blast->hits ) {
    $hit_count++;
    print $out $blast;
  }
  else {
    print STDERR "Hitless Blast Report: $blast ";
    print STDERR ($blast->no_hits_found ? "\n" : "(filtered)\n");
  }
}
# Can't do this since $in isn't a Bio::SearchIO object reference:
#    printf STDERR "\n%d Blast report(s) processed.\n", $in->report_count;
#    printf STDERR "Output sent to file: %s\n",  $out->file if $out->file;
printf STDERR "\nOutput sent to file: hitwriter3.out\n";


