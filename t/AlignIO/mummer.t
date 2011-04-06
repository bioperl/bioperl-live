#!/usr/bin/perl

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin();
	
    use_ok('Bio::Seq::RefLocatableSeq');
	use_ok('Bio::AlignIO::mummer');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

use Bio::AlignIO::mummer;
use Bio::AlignIO;
use Data::Dumper;

#alignio test
my $in = Bio::AlignIO->new(-file => test_input_file('mummer_ex1.mum'),
                           -format => 'mummer');
#my $out = Bio::AlignIO->new(-file => shift @ARGV || ">output.mum",
#                            -format => 'mummer');

while(my $aln = $in->next_aln()) {
    ok($aln);
}

#direct instantiation test
#my $test2 = Bio::AlignIO::mummer->new();
#$test2->write_aln();

#compat test
#my $test3 = new Bio::Seq::RefLocatableSeq(-id => 'test');
#print Dumper($test3) . "\n";
#print "compatible!\n" if $test3->isa('Bio::LocatableSeq');

done_testing();
