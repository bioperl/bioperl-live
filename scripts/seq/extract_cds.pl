#!/usr/bin/perl -w
# Contributed by Jason Stajich <jason@bioperl.org>

# simple extract the CDS features from a genbank file and 
# write out the CDS and Peptide sequences

use strict;
use Bio::SeqIO;
my $filename = shift || die("pass in a genbank filename on the cmd line");
my $in = new Bio::SeqIO(-file => $filename, -format => 'genbank');
my $out = new Bio::SeqIO(-file => ">$filename.cds");
my $outpep = new Bio::SeqIO(-file => ">$filename.pep");

while( my $seq = $in->next_seq ) {
  my @cds = grep { $_->primary_tag eq 'CDS' } $seq->get_SeqFeatures();
  foreach my $feature ( @cds ) {
    my $featureseq = $feature->spliced_seq;
    $out->write_seq($featureseq);
    $outpep->write_seq($featureseq->translate);
  }
}
