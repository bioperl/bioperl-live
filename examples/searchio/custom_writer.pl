#!/usr/bin/env perl 

# Demonstrates the use of a SearchIO Blast parser and a SearchWriterI object
# for producing custom output of Blast hit data from a Blast report 
# input stream.
#
# Here we define a custom SearchWriterI object that ouputs just the data we want
# from each BLAST report.
# 
# NOTE: If you just want pick and choose which columns you want 
# in the output table, you don't need to create your own custom
# SearchWriterI implementation as we do here. HitTableWriter and HSPTableWriter
# are configurable as to what columns and order you want. 
# The hitwriter*.pl and hspwriter.pl examples in this directory
# illustrate this.
#
# For a complete list of columns, see the docs for these modules:
#   Bio::SearchIO::Writer::HitTableWriter
#   Bio::SearchIO::Writer::HSPTableWriter
#
# This example serves as an illustration of how to use the 
# SearchWriterI api and plug it in to a SearchIO parser,
# which you may want to do if you want to generate data column(s)
# not provided by the available writers.
#
# Usage:
#   STDIN:  stream containing one or more BLAST or PSI-BLAST reports.
#   STDOUT: none, but generates an output file "custom_writer.out"
#           containing tab-delimited data on a per-hit basis.
#   STDERR: Progress info.
#
# Author: Steve Chervitz <sac@bioperl.org>

package MyBlastWriter;

use strict;
use lib '../../';
use Bio::Root::Root;
use Bio::SearchIO::SearchWriterI;

use base qw( Bio::Root::Root Bio::SearchIO::SearchWriterI );

sub to_string {
    my ($self, $result, @args) = @_;
    my $str = '';

    my $hits_reported = 0;

    foreach my $hit($result->hits) {

      # If this is a PSI-BLAST report, only report novel hits
      if( $result->psiblast ) {
	# Note that we could have supplied this has a -HIT_FILTER function
	# when we defined our input SearchIO object. Then we wouldn't need 
	# to define a custom writer.
	next unless $hit->iteration > 1 and not $hit->found_again;
      }

      $hits_reported++;
      printf STDERR "$hit\n";

      $str .= sprintf "%s\t%d\t%s\t%d\t%.2f\t%d\t%.1e\t%d\t%d\t%d\t%d\t%s\n", 
	               $result->query_name, $result->query_length, $hit->name,
	               $hit->length, $hit->frac_identical('query'), $hit->length_aln,
	               $hit->expect, $hit->score, $hit->bits,
	               $hit->gaps('total'), $hit->num_hsps, $hit->iteration || '-';
    }

    printf STDERR "\n%d hits written\n", $hits_reported;

    $str;

}

package main;

#===================================================
# Start of script 
#===================================================

use strict;

use lib '../../../';
use Bio::SearchIO;

select STDOUT; $|=1; 

my $in     = Bio::SearchIO->new( -format => 'blast', 
				 -fh => \*ARGV,
                                 -signif => 0.1 );
my $writer = MyBlastWriter->new();
my $out    = Bio::SearchIO->new( -format => 'blast',
				 -writer => $writer,
				 -file   => ">custom_writer.out" );

while ( my $result = $in->next_result() ) {
    printf STDERR "Report %d: $result\n", $in->result_count;
    $out->write_result($result);
}

printf STDERR "\n%d Results processed.\n", $in->result_count;
printf STDERR "Output sent to file: %s\n",  $out->file if $out->file;

