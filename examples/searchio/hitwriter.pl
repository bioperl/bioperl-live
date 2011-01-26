#!/usr/bin/perl 

# Demonstrates the use of a SearchIO Blast parser and a SearchWriterI object
# for producing tab-delimited output of hit data from a Blast report 
# input stream.
#
# Each row in the output represents data for a single hit.
# For hits containing multiple HSPs, the output information represents a
# summary across all HSPs.
#
# This parser represents a new and improved version of Bio::Tools::Blast.
#
# Usage:
#   STDIN:  stream containing one or more BLAST or PSI-BLAST reports.
#   STDOUT: none, but generates an output file "hitwriter.out"
#           containing tab-delimited data on a per-hit basis.
#   STDERR: Progress info and any errors.
# 
# In this example, we create a SearchIO parser that screens out hits 
# based on expect (or P) scores and a default HitTableWriter. This writer
# provides the same functionality as the original Bio::Tools::Blast::table()
# function (i.e., a tab-delimited summary of each hit per row).
# HitTableWriter, however, is customizable so you can specify just the columns
# you want to have in the output table.
#
# For more documentation about the writer, including
# a complete list of columns, execute:
#   perldoc Bio::SearchIO::Writer::HitTableWriter.
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
use Bio::SearchIO::Writer::HitTableWriter;

# These are the columns that will be in the output table of BLAST results.
my @columns = qw(
		 query_name
		 query_length
                 hit_name
                 hit_length
		 num_hsps
                 expect
                 frac_aligned_query
                 frac_identical_query
                 length_aln_query
                 gaps_total
                 strand_query
                 strand_hit
		);

# The following columns require HSP alignment data:
# 		  num_hsps
#                 frac_identical_query
#                 length_aln_query
#                 gaps_total
#                 strand_query
#                 strand_hit

print STDERR "\nUsing SearchIO->new()\n";

# Note that all parameters for the $in, $out, and $writer objects are optional.
# Default in = STDIN; Default out = STDOUT; Default writer = all columns 
# In this example, we're reading from STDIN and  writing to a file 
# called "hitwriter.out"
# TODO: write hitless reports to STDERR and note if filtered.
my $in     = Bio::SearchIO->new( -format => 'blast', 
				 -fh => \*ARGV,
				 -signif => 0.1, 
				# -verbose=> 2
                               );
my $writer = Bio::SearchIO::Writer::HitTableWriter->new( -columns => \@columns
						       );
my $out    = Bio::SearchIO->new( -format => 'blast',
				 -writer => $writer,
				 -file   => ">hitwriter.out" );
# Need to keep a separate count of reports with hits
# to know when to include labels. The first report may be hitless, 
# so we can't use $in->result_count
my $hit_count = 0;
while ( my $blast = $in->next_result() ) {
  printf STDERR "\nReport %d: $blast\n", $in->result_count;
  
  printf STDERR "query=%s, length=%d\n", $blast->query_name, $blast->query_length;

  if( $blast->hits ) {
      print STDERR "# hits= ", $blast->num_hits, "\n";
      $hit_count++;
      my @hits= $blast->hits;
      print STDERR "frac_aligned_query= ", $hits[0]->frac_aligned_query, "\n";

      $out->write_result($blast, $hit_count==1 );
  }
  else {
    print STDERR "Hitless Blast Report ";
    print STDERR ($blast->no_hits_found ? "\n" : "(filtered)\n");
  }
  
  ## For a simple progress monitor, uncomment this line:
  #print STDERR "."; print STDERR "\n" if $in->result_count % 50 == 0;
}

printf STDERR "\n%d Blast report(s) processed.\n", $in->result_count;
printf STDERR "Output sent to file: %s\n",  $out->file if $out->file;

