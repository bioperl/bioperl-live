#!/usr/local/bin/perl 

# Example usage of a SearchIO::psiblast parser of traditional format Blast reports.
# Same as hitwriter.pl but showing how to use the -hit_filter option.
#
# Usage:
#   STDIN:  stream containing one or more BLAST or PSI-BLAST reports.
#   STDOUT: none, but generates an output file "hitwriter2.out"
#           containing tab-delimited data on a per-hit basis.
#   STDERR: Progress info.
#
# For more documentation about the writer, including
# a complete list of columns, execute:
#   perldoc Bio::SearchIO::Writer::HitTableWriter.
#
# For more documentation about working with Blast result objects,
# see docs foor these modules:
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
use Bio::SearchIO::Writer::HitTableWriter;

# These are the columns that will be in the output table of BLAST results.
my @columns = qw(
		 query_name
		 query_length
                 hit_name
                 hit_length
		 num_hsps
                 expect
                 frac_identical_query
                 length_aln_query
                 gaps_total
                 strand_query
                 strand_hit
		);

# These labels just override the defaults for the specified columns.
# (first column is 1).
my %labels = (
	      1 => 'QUERY_GI',
	      3 => 'HIT_IDENTIFIER',
	     );


# BLAST hit filtering function. All hits of each BLAST report must satisfy 
# this criteria to be retained. If a hit fails this test, it is ignored.
# If all hits of a report fail, the report will be considered hitless,
# but we can output a "(filtered)" string to indicate that this is the reason.
my $filt_func = sub{ my $hit=shift; 
                     $hit->frac_identical('query') >= 0.5 
                         && $hit->frac_aligned_query >= 0.50
                     };

print STDERR "\nUsing SearchIO->new()\n";

# Note that all parameters for the $in, $out, and $writer objects are optional.
# Default in = STDIN; Default out = STDOUT; Default writer = all columns 
# In this example, we're reading from STDIN and  writing to a file 
# called "hitwriter.out"

my $in     = Bio::SearchIO->new( -format => 'psiblast',
				 -hit_filter => $filt_func,);
my $writer = Bio::SearchIO::Writer::HitTableWriter->new( -columns => \@columns, 
							 -labels  => \%labels 
						       );
my $out    = Bio::SearchIO->new( -format => 'psiblast',
				 -writer => $writer,
				 -file   => ">hitwriter2.out" );
my $hit_count = 0;
while ( my $blast = $in->next_result() ) {
  printf STDERR "\nReport %d: $blast\n", $in->report_count;
  
  if( $blast->hits ) {
    $hit_count++;
    $out->write_result($blast, $hit_count==1 );
  }
  else {
    print STDERR "Hitless Blast Report: ";
    print STDERR ($blast->no_hits_found ? "\n" : "(filtered)\n");
  }
  
  ## For a simple progress monitor, uncomment this line:
  #print STDERR "."; print STDERR "\n" if $in->report_count % 50 == 0;
}

printf STDERR "\n%d Blast report(s) processed.\n", $in->report_count;
printf STDERR "Output sent to file: %s\n",  $out->file if $out->file;

