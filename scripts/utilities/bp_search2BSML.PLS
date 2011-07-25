#!perl

# Author:      Jason Stajich <jason-at-bioperl-dot-org>
# Description: Turn SearchIO parseable report(s) into a GFF report
#
=head1 NAME

search2bsml - Turn SearchIO parseable reports(s) into a BSML report

=head1 SYNOPSIS

Usage:
  search2bsml [-o outputfile] [-f reportformat] [-i inputfilename]  OR file1 file2 ..

=head1 DESCRIPTION

This script will turn a protein Search report (BLASTP, FASTP, SSEARCH, 
AXT, WABA, SIM4) into a BSML File.

The options are:

   -i infilename        - (optional) inputfilename, will read
                          either ARGV files or from STDIN
   -o filename          - the output filename [default STDOUT]
   -f format            - search result format (blast, fasta,waba,axt)
                          (ssearch is fasta format). default is blast.
   -h                   - this help menu

Additionally specify the filenames you want to process on the
command-line.  If no files are specified then STDIN input is assumed.
You specify this by doing: search2gff E<lt> file1 file2 file3

=head1 AUTHOR

Jason Stajich, jason-at-bioperl-dot-org

=cut

use strict;
use warnings;
use Getopt::Long;
use Bio::SearchIO;

my ($output,$input,$format,$type,$help,$cutoff);
$format = 'blast'; # by default
GetOptions(
	   'i|input:s'  => \$input,
	   'o|output:s' => \$output,
	   'f|format:s' => \$format,
	   'c|cutoff:s' => \$cutoff,
	   'h|help'     => sub{ exec('perldoc',$0);
				exit(0)
				},
	   );
# if no input is provided STDIN will be used
my $parser = new Bio::SearchIO(-format => $format, 
			       -file   => $input);

my $out;
if( defined $output ) {
    $out = new Bio::SearchIO(-file => ">$output",
			     -output_format => 'BSMLResultWriter');
} else { 
    $out = new Bio::SearchIO(-output_format => 'BSMLResultWriter'); # STDOUT
}

while( my $result = $parser->next_result ) {
    $out->write_result($result);
}

