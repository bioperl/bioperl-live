#!/usr/bin/perl
use strict;
use warnings;
use Bio::Root::Test;
use File::Spec::Functions qw(catfile);
use FindBin;
test_begin(-tests => 1);

use Bio::SeqIO;

my $seqio = Bio::SeqIO->new(
    -file   => catfile($FindBin::Bin, 'test_clear_range.fastq'),
    -format => 'fastq'
);

while ( my $seq = $seqio->next_seq() ) {
    $seq->threshold(15);
    lives_ok { my $newqualobj = $seq->get_clear_range };
}
