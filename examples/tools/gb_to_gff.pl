#!/usr/local/bin/perl
use strict;

use Bio::Tools::GFF;
use Bio::SeqIO;

my ($seqfile) = @ARGV;
die("must define a valid seqfile to read") unless ( defined $seqfile && -r $seqfile);

my $seqio = new Bio::SeqIO(-format => 'genbank',
			   -file   => $seqfile);
my $count = 0;
while( my $seq = $seqio->next_seq ) {
    $count++;
    # defined a default name
    my $fname = sprintf("%s.gff", $seq->display_id || "seq-$count");
    my $gffout = new Bio::Tools::GFF(-file => ">$fname" ,
				     -gff_version => 1);
    
    foreach my $feature ( $seq->top_SeqFeatures() ) {
	$gffout->write_feature($feature);
    }
}
