#!/usr/bin/perl -w
# $Id: search2gff.PLS,v 1.11 2004-07-19 01:11:02 jason Exp $

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
   -t/--type seqtype    - if you want to see query or hit information
                          in the GFF report
   -s/--source          - specify the source (will be algorithm name
                          otherwise like BLASTN)
   -h                   - this help menu
   --version            - GFF version to use (put a 3 here to use gff 3)
   --component          - generate GFF component fields (chromosome)
   -m/--match           - generate a 'match' line which is a container
                          of all the similarity HSPs
   -c/--cutoff          - specify an evalue cutoff

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

my ($output,$input,$format,$type,$help,$cutoff,$sourcetag,$comp,
    $gffver,$match,$quiet);
$format = 'blast'; # by default
$type   = 'query';
$gffver = 2;
GetOptions(
	   'i|input:s'  => \$input,
	   'component'  => \$comp,
	   'm|match'    => \$match,
	   'o|output:s' => \$output,
	   'f|format:s' => \$format,
	   's|source:s' => \$sourcetag,
	   't|type:s'   => \$type,
	   'c|cutoff:s' => \$cutoff,
	   'v|version:i'=> \$gffver,
           'q|quiet'    => \$quiet,
	   'h|help'     => sub{ exec('perldoc',$0);
				exit(0)
				},
	   );
$type = lc($type);
if( $type =~ /target/ ) { $type = 'hit' }
elsif( $type ne 'query' && $type ne 'hit' ) {
    die("seqtype must be either 'query' or 'hit'");
} 
# if no input is provided STDIN will be used
my $parser = new Bio::SearchIO(-format => $format,
			       -verbose => $quiet ? -1 : 0,
			       -file   => $input);

my $out;
if( defined $output ) {
    $out = new Bio::Tools::GFF(-gff_version => $gffver,
			       -file => ">$output");
} else { 
    $out = new Bio::Tools::GFF(-gff_version => $gffver); # STDOUT
}
my (%seen_hit,%seen);
my $other = $type eq 'query' ? 'hit' : 'query';

while( my $result = $parser->next_result ) {
    my $qname = $result->query_name;
    if ( $comp && $type eq 'query' && 
	 $result->query_length ) {
	$out->write_feature(Bio::SeqFeature::Generic->new
			    (-start       => 1,
			     -end         => $result->query_length,
			     -seq_id      => $qname,
			     -source_tag  => 'chromosome',
			     -primary_tag => 'Component',
			     -tag         => {
				 'Sequence' => $qname
				 }));
    }
    while( my $hit = $result->next_hit ) {
	next if( defined $cutoff && $hit->significance > $cutoff);
	my $acc = $qname;
	if( $seen{$qname."-".$hit->name}++ ) {
	    $acc = $qname."-".$seen{$qname.'-'.$hit->name};
	}
	
	if( $comp && $type eq 'hit' && $hit->length &&
	    ! $seen_hit{$hit->name}++ ) {
	    $out->write_feature(Bio::SeqFeature::Generic->new
				(-start       => 1,
				 -end         => $hit->length,
				 -seq_id      => $hit->name,
				 -source_tag  => 'chromosome',
				 -primary_tag => 'Component',
				 -tag         => {
				     'Sequence' => $hit->name
				     }));
	}
	my (%min,%max,$seqid,$name,$st);
	while( my $hsp = $hit->next_hsp ) {
	    my $feature = new Bio::SeqFeature::Generic;
	    my ($proxyfor,$otherf);
	    if( $type eq 'query' ) {
		($proxyfor,$otherf) = ($hsp->query,
				      $hsp->hit);
		$name  ||= $hit->name;
	    } else {
		($other,$proxyfor) = ($hsp->query,
				      $hsp->hit);
		$name ||= $acc;
	    }
	    $proxyfor->score($hit->bits) unless( $proxyfor->score );
	    if( $match && $gffver == 3) {
		$feature->add_tag_value('Parent', $name);
	    }
	    
	    $min{$type} = $proxyfor->start unless defined $min{$type} && $min{$type} < $proxyfor->start;
	    $max{$type} = $proxyfor->start unless defined $max{$type} && $max{$type} > $proxyfor->end;
	    $min{$other} = $otherf->start unless defined $min{$type} && $min{$type} < $otherf->start;
	    $max{$other} = $otherf->start unless defined $max{$type} && $max{$type} > $otherf->end;
	    
	    $feature->add_tag_value('Target', 'Sequence:'.$name);
	    $feature->add_tag_value('Target', $otherf->start);
	    $feature->add_tag_value('Target', $otherf->end);

	    $feature->location($proxyfor->location);
	    #  strand for feature is always going to be product of
	    #  query & hit strands so that target can always be just
	    #  '+'
	    $feature->strand ( $proxyfor->strand * $otherf->strand);
	    if( $sourcetag ) { 
		$feature->source_tag($sourcetag);
	    } else {
		$feature->source_tag($proxyfor->source_tag);
	    }
	    $feature->score($proxyfor->score);
	    $feature->frame($proxyfor->frame);
	    $feature->seq_id($proxyfor->seq_id );
	    $feature->primary_tag('similarity');
	    $seqid ||= $proxyfor->seq_id;
	    $out->write_feature($feature);
	    $st ||= $sourcetag || $proxyfor->source_tag;
	}
	if( $match ) { 
	    
	    my $matchf = Bio::SeqFeature::Generic->new
		(-start => $min{$type},
		 -end   => $max{$type},
		 -strand=> $hit->strand($type)*$hit->strand($other),
		 -primary_tag => 'match',
		 -source_tag  => $st,
		 -score => $hit->bits,
		 -seq_id => $seqid);
	    if( $gffver == 3 ) { 
		$matchf->add_tag_value('ID', $name);
	    }
	    $matchf->add_tag_value('Target', "Sequence:$name");
	    $out->write_feature($matchf);
	}
    }
}
