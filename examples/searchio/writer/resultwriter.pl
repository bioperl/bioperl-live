#!/usr/bin/perl 

# Example usage of a SearchIO-based parser of existing Blast reports.
# This parser represents a new and improved version of Bio/Tools/Blast.pm.
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
use Bio::SearchIO::Writer::ResultTableWriter;
use Error qw(:try);

print "\nUsing SearchIO->new()\n";

try {

    # Note that all parameters for the $in, $out, and $writer objects are optional.
    # Default in = STDIN; Default out = STDOUT; Default writer = all columns
    # In this example, we're reading from STDIN and  writing to STDOUT
    # and using the default columns for the writer.
    # We're also telling the script to timeout if input isn't received
    # within 10 sec. (Note the clock is still ticking when you background the job.)
    # Setting verbose to 1 is useful for debugging.
    my $in = Bio::SearchIO->new( -format => 'psiblast', 
				 -signif => 0.1, 
				 -verbose => 0, 
				 -timeout_sec => 10 );
    # not specifying any columns to get the default.
    my $writer = Bio::SearchIO::Writer::ResultTableWriter->new();
    my $out    = Bio::SearchIO->new( -format => 'psiblast', 
				     -writer => $writer,
                                     -file => ">resultwriter.out");

    while ( my $result = $in->next_result() ) {
#       printf STDERR "Report %d: $result\n", $in->report_count;
        $out->write_result($result, ($in->report_count - 1 ? 0 : 1) );

	# To get at the statistical parameters:
        # Calling raw_statistics() returns a list containing the
        # unparsed lines of the parameters section of the report.
        # Here we're only interested in parameters beginning with "effective".
	 print "Report Stats, effective data:\n";
	 foreach( $result->raw_statistics) {
	     print "$_" if /^effective/i;
	 }

	## For a simple progress monitor, uncomment this line:
	#print STDERR "."; print STDERR "\n" if $in->report_count % 50 == 0;
    }

    printf STDERR "\n%d Blast report(s) processed.\n", $in->report_count;
    printf STDERR "Output sent to file: %s\n",  $out->file if $out->file;

}
catch Bio::Root::Exception with {
    my $err = shift;
    print STDERR "\nCaught Bio::Root::Exception:\n\n$err\n";
}
otherwise {
    my $err = shift;
    print STDERR "\nAn unanticipated exception occurred of type ", ref $err, "\n\n";
    print STDERR "$err\n";
};
