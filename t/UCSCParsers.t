# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##

my $error;
use strict;
BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    $error = 0; 
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    use vars qw($TESTCOUNT);
    $TESTCOUNT = 48;
    plan tests => $TESTCOUNT;
}

use Bio::SearchIO;
use Bio::Root::IO;

my $pslparser = new Bio::SearchIO(-format => 'psl',
				  -file   => Bio::Root::IO->catfile
				  (qw(t data sbay_c545-yeast.BLASTZ.PSL)));

my $result = $pslparser->next_result;
ok($result->query_name, 'I');
ok($result->query_length, 230203);

my $hit    = $result->next_hit;
ok($hit->name, 'sbay_c545');
ok($hit->length, 28791);
my $hsp    = $hit->next_hsp;
ok($hsp->query->start,139871);
ok($hsp->query->end,141472);
ok($hsp->query->length, 1602);
ok($hsp->query->strand, 1);
ok($hsp->hit->strand, 1);
my $q_gapblocks = $hsp->gap_blocks('query');
ok(scalar @$q_gapblocks, 24);
ok($q_gapblocks->[0]->[1],45);
ok($q_gapblocks->[1]->[1],10);
ok($q_gapblocks->[1]->[0],139921);


$hsp       = $hit->next_hsp;
$hsp       = $hit->next_hsp;
ok($hsp->hit->start,27302);
ok($hsp->hit->end,27468);
ok($hsp->hit->length,167);
ok($hsp->query->start, 123814);
ok($hsp->query->end, 123972);
ok($hsp->query->length, 159);
ok($hsp->query->strand,-1);

$q_gapblocks = $hsp->gap_blocks('query');
ok(scalar @$q_gapblocks, 4);
ok($q_gapblocks->[0]->[1],116);
ok($q_gapblocks->[1]->[1],4);
ok($q_gapblocks->[1]->[0],123856);



#-----------------------------------


$pslparser = new Bio::SearchIO(-format => 'psl',
			       -file   => Bio::Root::IO->catfile
			       (qw(t data blat.psLayout3)));

$result = $pslparser->next_result;
ok($result->query_name, 'sequence_10');
ok($result->query_length, 1775);

$hit    = $result->next_hit;
ok($hit->name, 'sequence_10');
ok($hit->length, 1775);
$hsp    = $hit->next_hsp;
ok($hsp->query->start,1);
ok($hsp->query->end,1775);
ok($hsp->query->length,1775);
ok($hsp->query->strand,1);
ok($hsp->hit->strand,1);
$q_gapblocks = $hsp->gap_blocks('query');
ok(scalar @$q_gapblocks, 1);
ok($q_gapblocks->[0]->[1],1775);
ok($q_gapblocks->[1]->[1],undef);
ok($q_gapblocks->[1]->[0],undef);


$hsp       = $hit->next_hsp;
ok($hsp->hit->start,841);
ok($hsp->hit->end,1244);
ok($hsp->query->start, 841);
ok($hsp->query->end, 1244);
ok($hsp->query->length, 404);
ok($hsp->query->strand,-1);
ok($hsp->hit->strand,1);

$q_gapblocks = $hsp->gap_blocks('query');
ok(scalar @$q_gapblocks, 4);
ok($q_gapblocks->[0]->[1],14);
ok($q_gapblocks->[1]->[1],21);
ok($q_gapblocks->[1]->[0],1152);


