# -*-Perl-*- Test Harness script for Bioperl
# $Id: SearchIO_cross_match.t 11788 2007-12-03 23:37:59Z jason $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 15);
	
	use_ok('Bio::SearchIO');
}

my ($searchio, $result,$iter,$hit,$hsp);

# The cross_match SearchIO parser is not event-based; it directly adds
# information to the relevant Bio::Search objects as the report is parsed.
# The parser currently misses much information present in the report. Also,
# methods expected to work somehow don't (hsp->length('hsp'), for instance).
# Unsure if this parses non-alignment-based cross-match reports accurately
# (see bioperl-live/t/data/consed_project/edit_dir/test_project.screen.out for
# an example).

# Note lots of ResultI/HitI/HSPI methods not tested yet!

$searchio = Bio::SearchIO->new('-format' => 'cross_match',
				  '-file'   => test_input_file('testdata.crossmatch'));

$result = $searchio->next_result;

is($result->algorithm, 'cross_match');
is($result->algorithm_version, '0.990329');

my @valid = ( [ 'msx1_ens2', 0]);
my $count = 0;
while( $hit = $result->next_hit ) {
    my $d = shift @valid;

    is($hit->name, shift @$d);
    is($hit->length, shift @$d);

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 19);
            is($hsp->query->end, 603);
            is($hsp->hit->start, 2824);
            is($hsp->hit->end, 3409);
            #is($hsp->length('hsp'), 820);  # shouldn't this work?
            is($hsp->start('hit'), $hsp->hit->start);
            is($hsp->end('query'), $hsp->query->end);
            is($hsp->strand('sbjct'), $hsp->subject->strand);# alias for hit
            is($hsp->gaps, 0);
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}

is(@valid, 0);
