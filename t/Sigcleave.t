# -*-Perl-*-
## Bioperl Test Harness Script for various modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    plan tests => 16;
}
use Bio::PrimarySeq;
use Bio::Tools::Sigcleave;

#load n-terminus of MGR5_HUMAN as test seq
my $protein = "MVLLLILSVLLLKEDVRGSAQSSERRVVAHMPGDIIIGALFSVHHQPTVDKVHERKCGAVREQYGI";

ok my $seq= Bio::PrimarySeq->new(-seq => $protein);

ok my $sig = Bio::Tools::Sigcleave->new();
ok $sig->seq($seq);
ok my $sout = $sig->seq;
ok $sout->seq eq $protein;
ok $sig->threshold, 3.5;
ok $sig->threshold(5), 5;
ok $sig->matrix, 'eucaryotic';
ok $sig->matrix('procaryotic'), 'procaryotic';
ok $sig->matrix('eucaryotic'), 'eucaryotic';

ok $sig->pretty_print =~ /Maximum score 7/;
ok my %results = $sig->signals;

ok $results{9}, 5.2, "unable to get raw sigcleave results";


$sig = Bio::Tools::Sigcleave->new(-seq=>$protein,
				 -threshold=>5);
ok %results = $sig->signals;
ok $results{9}, 5.2, "unable to get raw sigcleave results";
ok $sig->result_count, 5;

