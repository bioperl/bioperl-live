#!/usr/bin/perl
use strict;
use Bio::SearchIO;
use Getopt::Long;
use Benchmark;

my ($infile,$outfile,$verbose);

GetOptions( 
	    'i|input:s'  => \$infile,
	    'o|output:s' => \$outfile,
	    'v|verbose'  => \$verbose,
	    );
$infile = shift unless $infile;

my $in = new Bio::SearchIO(-format => 'waba',
			   -file   => $infile, #comment out to read from STDIN
                           #-fh => \*ARGV,  # uncomment to read from STDIN
			   -verbose => $verbose);

my $out;
if( defined $outfile) {
    $out = new Bio::Tools::GFF(-file => ">$outfile");
} else {
    $out = new Bio::Tools::GFF(-verbose => $verbose);
}

while( my $r = $in->next_result ) {
    while( my $hit = $r->next_hit ) {
	while( my $hsp = $hit->next_hsp ) {
	    $out->write_feature($hsp);
	}
    }
}
