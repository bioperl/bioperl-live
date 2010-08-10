#!/usr/bin/perl

use lib "../src/";
use Bio::Seq::RefLocatableSeq;
use Bio::AlignIO::mummer;
use Bio::AlignIO;
use Data::Dumper;

#alignio test
my $in = Bio::AlignIO->new(-file => shift @ARGV || "../data/mummer_ex1.mum",
                           -format => 'mummer');
my $out = Bio::AlignIO->new(-file => shift @ARGV || ">output.mum",
                            -format => 'mummer');

while(my $aln = $in->next_aln()) {
   $out->write_aln($aln);
}

#direct instantiation test
#my $test2 = Bio::AlignIO::mummer->new();
#$test2->write_aln();

#compat test
#my $test3 = new Bio::Seq::RefLocatableSeq(-id => 'test');
#print Dumper($test3) . "\n";
#print "compatible!\n" if $test3->isa('Bio::LocatableSeq');

