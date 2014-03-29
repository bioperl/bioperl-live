#!perl

# Author:      Jason Stajich <jason-at-bioperl-dot-org>
# Description: Turn SearchIO parseable report(s) into a TRIBE matrix
#
=head1 NAME

bp_search2tribe - Turn SearchIO parseable reports(s) into TRIBE matrix

=head1 SYNOPSIS

Usage:
  bp_search2tribe [-o outputfile] [-f reportformat] [-w/--weight] file1 file2 ..

=head1 DESCRIPTION

This script is probably too slow for most people's uses.  It is better
to use something like scripts/searchio/fastam9_to_table, -m 9 output
from BLAST, or the blast2table from the BLAST O'Reilly book to get a
tabular output from these programs and then feed the table into MCL
with the mcxdeblast script and the --m9 option.

This script will turn a protein Search report (BLASTP, FASTP, SSEARCH)
into a Markov Matrix for TribeMCL clustering.

The options are:

   -o filename          - the output filename [default STDOUT]
   -f format            - search result format (blast, fasta)
                          (ssearch is fasta format). default is blast.
   -w or --weight VALUE - Change the default weight for E(0.0) hits
                          to VALUE (default=200 (i.e. 1e-200) )
   -h                   - this help menu

Additionally specify the filenames you want to process on the
command-line.  If no files are specified then STDIN input is assumed.
You specify this by doing: bp_search2tribe E<lt> file1 file2 file3

=head1 AUTHOR

Jason Stajich, jason-at-bioperl-dot-org

=cut

use strict;
use warnings;
use Bio::SearchIO;
use Bio::SearchIO::FastHitEventBuilder; # employ a speedup
use Getopt::Long;
use constant DEFAULT_WEIGHT => 200;
use constant DEFAULT_FORMAT => 'blast';

my ($format,@files,$output,$weight);
$weight = DEFAULT_WEIGHT; # default weight value
$format = DEFAULT_FORMAT;

my ($help);

GetOptions(
	   'f|format:s'    => \$format,
	   'o|output:s'    => \$output,
	   'w|weight:i'  => \$weight,
	   'h|help'        => sub{ exec('perldoc',$0);
				   exit(0)
				   },
	   );

my $outfh;
if( $output ) { 
    open $outfh, '>', $output or die "Could not write output file '$output': $!\n";
} else {
    $outfh = *STDOUT;
}

my $parser = new Bio::SearchIO(-format => $format, -fh => \*ARGV);

# Let's throw away HSP events
$parser->attach_EventHandler(new Bio::SearchIO::FastHitEventBuilder);
while( my $report = $parser->next_result ) {
    my $q = $report->query_name;
    while( my $hit = $report->next_hit ) {
	my $evalue = $hit->significance;
	$evalue =~ s/^e/1e/i;

	if( $evalue == 0 ) {	    
	    $evalue = "1e-$weight";
	} else { 
	    $evalue = sprintf("%e",$evalue);
	}

	print $outfh join("\t",$q,$hit->name, split('e-',$evalue)), "\n"; 
    }
}
