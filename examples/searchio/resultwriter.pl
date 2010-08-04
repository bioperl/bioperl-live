#!/usr/bin/perl 

# Demonstrates the use of a SearchIO Blast parser and a SearchWriterI object
# for producing tab-delimited output of result data from a Blast report 
# input stream.
#
# This writer only outputs information at the level of the result object.
# This shows that you can work with a writer that only knows about 
# Bio::Search::Result objects and doesn't care about hit or HSP data. 
# Therefore, the output from this example doesn't contain any information 
# about hits or HSPs. 
# See the hitwriter.pl and hspwriter.pl examples for that.
#
# This parser represents a new and improved version of Bio::Tools::Blast.
#
# Usage:
#   STDIN:  stream containing one or more BLAST or PSI-BLAST reports.
#   STDOUT: none, but generates an output file "resultwriter.out"
#           containing tab-delimited data on a per-report basis.
#   STDERR: Any errors that occurred.
#
# For more documentation about the writer, including
# a complete list of columns, see the docs for 
#   Bio::SearchIO::Writer::ResultTableWriter.
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
use Bio::SearchIO::Writer::ResultTableWriter;
use Bio::SearchIO::Writer::HTMLResultWriter;

print "\nUsing SearchIO->new()\n";


# Note that all parameters for the $in, $out, and $writer objects are optional.
# Default in = STDIN; Default out = STDOUT; Default writer = all columns
# In this example, we're reading from STDIN and  writing to STDOUT
# and using the default columns for the writer.
# We're also telling the script to timeout if input isn't received
# within 10 sec. (Note the clock is still ticking when you background the job.)
# Setting verbose to 1 is useful for debugging.
my $in = Bio::SearchIO->new( -format => 'blast', 
                             -fh => \*ARGV,
                             -signif => 0.1, 
                             -verbose => 0, 
                             -timeout_sec => 10 );
# not specifying any columns to get the default.
my $writer = Bio::SearchIO::Writer::ResultTableWriter->new();
my $out    = Bio::SearchIO->new( -format => 'blast', 
                                 -writer => $writer,
                                 -file => ">resultwriter.out");

my $writerhtml = new Bio::SearchIO::Writer::HTMLResultWriter();
my $outhtml = new Bio::SearchIO(-writer => $writerhtml,
                                -file   => ">searchio.html");


while ( my $result = $in->next_result() ) {
    eval {
        # printf STDERR "Report %d: $result\n", $in->result_count;
        $out->write_result($result, ($in->result_count - 1 ? 0 : 1) );
        
        $outhtml->write_result($result, 1);

        # To get at the statistical parameters:
        # Calling raw_statistics() returns a list containing the
        # unparsed lines of the parameters section of the report.
        # Here we're only interested in parameters beginning with "effective".
        #         print "Report Stats, effective data:\n";
        #         foreach( $result->raw_statistics) {
        #             print "$_" if /^effective/i;
        #         }
        
	## For a simple progress monitor, uncomment this line:
	#print STDERR "."; print STDERR "\n" if $in->result_count % 50 == 0;
    };
    if($@) {
        warn "Warning: Blast parsing or writing exception caught for $result:\n$@\n";
    }
}

printf STDERR "\n%d Blast report(s) processed.\n", $in->result_count;
printf STDERR "Output sent to file: %s\n",  $out->file if $out->file;
