#!/usr/bin/perl -w
# $Id: search2gff.PLS,v 1.4 2003-05-19 17:27:09 jason Exp $

# Author:      Jason Stajich <jason@bioperl.org>
# Description: Turn SearchIO parseable report(s) into a GFF report
#
=head1 NAME

search2gff - Turn SearchIO parseable reports(s) into a GFF report

=head1 SYNOPSIS

Usage:
  search2gff [-o outputfile] [-f reportformat] [-i inputfilename]  OR file1 file2 ..

=head1 DESCRIPTION

This script will turn a protein Search report (BLASTP, FASTP, SSEARCH, 
AXT, WABA) into a GFF File.

The options are:

   -i infilename        - (optional) inputfilename, will read
                          either ARGV files or from STDIN
   -o filename          - the output filename [default STDOUT]
   -f format            - search result format (blast, fasta,waba,axt)
                          (ssearch is fasta format). default is blast.
   -t seqtype           - if you want to see query or hit information
                          in the GFF report
   -h                   - this help menu

Additionally specify the filenames you want to process on the
command-line.  If no files are specified then STDIN input is assumed.
You specify this by doing: search2gff E<lt> file1 file2 file3

=head1 AUTHOR

Jason Stajich, jason-at-bioperl-dot-org

=cut

use strict;
use Bio::Tools::GFF;
use Getopt::Long;
use Bio::SearchIO;

my ($output,$input,$format,$type,$help,$cutoff);
$format = 'blast'; # by default
$type   = 'query';
GetOptions(
	   'i|input:s'  => \$input,
	   'o|output:s' => \$output,
	   'f|format:s' => \$format,
	   't|type:s'   => \$type,
	   'c|cutoff:s' => \$cutoff,
	   'h|help'     => sub{ exec('perldoc',$0);
				exit(0)
				},
	   );
$type = lc($type);
if( $type ne 'query' && $type ne 'hit' ) {
    die("seqtype must be either 'query' or 'hit'");
} 
# if no input is provided STDIN will be used
my $parser = new Bio::SearchIO(-format => $format, 
			       -file   => $input);

my $out;
if( defined $output ) {
    $out = new Bio::Tools::GFF(-file => ">$output");
} else { 
    $out = new Bio::Tools::GFF(); # STDOUT
}
while( my $result = $parser->next_result ) {
    while( my $hit = $result->next_hit ) {
	next if( defined $cutoff && $hit->significance > $cutoff);
	while( my $hsp = $hit->next_hsp ) {
	    my $feature = new Bio::SeqFeature::Generic;
	    my $proxyfor;
	    if( $type eq 'query' ) {
		$proxyfor = $hsp->query;
		$feature->add_tag_value('Target', 'Sequence:'.$hit->name);
		if( $hsp->hit->strand < 0 ) {
		    $feature->add_tag_value('Target', $hsp->hit->end);
		    $feature->add_tag_value('Target', $hsp->hit->start);
		} else { 
		    $feature->add_tag_value('Target', $hsp->hit->start);
		    $feature->add_tag_value('Target', $hsp->hit->end);
		}
	    } else {
		$proxyfor = $hsp->hit;
		$feature->add_tag_value('Target', 'Sequence:'.
					$result->query_name);
		$proxyfor->score($hit->bits) unless( $proxyfor->score );

		if( $hsp->query->strand < 0 ) {
		    $feature->add_tag_value('Target', $hsp->query->end);
		    $feature->add_tag_value('Target', $hsp->query->start);
		} else { 
		    $feature->add_tag_value('Target', $hsp->query->start);
		    $feature->add_tag_value('Target', $hsp->query->end);
		}
	    }
	    $feature->location($proxyfor->location);
	    $feature->source_tag($proxyfor->source_tag);
	    $feature->score($proxyfor->score);
	    $feature->frame($proxyfor->frame);
	    $feature->seq_id($proxyfor->seq_id);
	    $feature->primary_tag('similarity');

	    $out->write_feature($feature);
	}
    }
}
