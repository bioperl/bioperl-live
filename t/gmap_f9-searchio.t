#!perl

use strict;
use warnings;

use Test::More tests => 46;

use Bio::SearchIO;
use Data::Dumper;

my $searchio =
    Bio::SearchIO->new(-format => 'gmap_f9',
		       -file   => 't/data/gmap_f9.txt');

my $result = $searchio->next_result;
isa_ok($result, 'Bio::Search::Result::GenericResult', 'Did we get a Result?');
is($result->num_hits(), 1, 'Did we get the expected number of hits?');
is($result->algorithm(), 'gmap', 'Did we get the expected algorithm?');
is($result->query_name(), 'NM_004448', 'Did we get the expected query_name?');

my $hit = $result->next_hit;
isa_ok($hit, 'Bio::Search::Hit::GenericHit', 'Did we get a Hit?');
_check_hit($hit, {name => '17',
		  length => 4624,
		  num_hsps => 27,
		 } );

my $hsp = $hit->next_hsp;
_check_hsp($hsp, {algorithm => 'GMAP',
		  query_gaps => 1,
		  hit_gaps => 0,
		  query_length => 310,
		  hit_length => 311,
		  qseq => 'GGAGGAGGTGGAGGAGGAGG', # first 20 bases
		  hseq => 'GGAGGAGGTGGAGGAGGAGG', # ditto
		  query => {start => 1,
			    end => 310,
			    strand => 1},
		  hit => {start => 35109780,
			  end =>  35110090,
			  strand => 1},
		 } );

my $searchio_rev =
    Bio::SearchIO->new(-format => 'gmap_f9',
		       -file   => 't/data/gmap_f9-reverse-strand.txt');
my $result_rev = $searchio_rev->next_result;
isa_ok($result_rev,
       'Bio::Search::Result::GenericResult', 'Did we get a Result?');
is($result_rev->num_hits(), 1, 'Did we get the expected number of hits?');
is($result_rev->algorithm(), 'gmap', 'Did we get the expected algorithm?');
is($result_rev->query_name(),
   'NM_004448', 'Did we get the expected query_name?');

$hit = $result_rev->next_hit;
_check_hit($hit, {name => '17',
		  length => 4624,
		  num_hsps => 27,
		 } );

$hsp = $hit->next_hsp;
_check_hsp($hsp, {algorithm => 'GMAP',
		  query_gaps => 0,
		  hit_gaps => 0,
		  query_length => 974,
		  hit_length => 974,
		  qseq => 'TAGCTGTTTTCCAAAATATA', # first 20 bases
		  hseq => 'TAGCTGTTTTCCAAAATATA', # ditto
		  query => {start => 1,
			    end => 974,
			    strand => 1},
		  hit => {start => 35137468,
			  end =>  35138441,
			  strand => -1},
		 } );


$searchio =  Bio::SearchIO->new(-format => 'gmap_f9',
				-file   => 't/data/gmap_f9-multiple_results.txt');

my $result_count = 0;
while (my $result = $searchio->next_result) {
    $result_count++;
#    while(my $hit = $result->next_hit){
#	my @hsps = $hit->hsps();
#	print $hit->name, " ", $hit->length, scalar(@hsps);
#    }
}

is($result_count, 58, "Can we loop over multiple results properly (expecting 58)?");

exit(0);


sub _check_hit {
    my ($hit, $info) = @_;

    isa_ok($hit, 'Bio::Search::Hit::HitI');
    is($hit->name, $info->{name}, 'Check the name');
    is($hit->length, $info->{length}, 'Check the hit length');
    is($hit->num_hsps, $info->{num_hsps}, 'Check the number of hsps');

}

sub _check_hsp {
    my($hsp, $info) = @_;
    isa_ok($hsp, 'Bio::Search::HSP::HSPI');
    is($hsp->algorithm, $info->{algorithm}, 'Check the algorithm');
    is($hsp->gaps('query'), $info->{query_gaps}, 'Count gaps in the query');
    is($hsp->gaps('hit'), $info->{hit_gaps}, 'Count gaps in the hit');
    is($hsp->length('query'), $info->{query_length}, 'Length of the query');
    is($hsp->length('hit'), $info->{hit_length}, 'Length of the hit');
    is(substr($hsp->query_string, 0, 20), $info->{qseq}, 'Query sequence');
    is(substr($hsp->hit_string, 0, 20), $info->{hseq}, 'Hit sequence');
    is($hsp->query->start, $info->{query}->{start}, "Check query start");
    is($hsp->query->end, $info->{query}->{end}, "Check query end");
    is($hsp->query->strand, $info->{query}->{strand}, "Check query end");
    is($hsp->hit->start, $info->{hit}->{start}, "Check hit start");
    is($hsp->hit->end, $info->{hit}->{end}, "Check hit end");
    is($hsp->hit->strand, $info->{hit}->{strand}, "Check hit end");
}
