use strict;

BEGIN {     
    
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;

    plan tests => 14;
	
	use_ok('Bio::Root::IO');
	use_ok('Bio::Tools::Tmhmm');
}



ok my $infile = Bio::Root::IO->catfile(qw(t data tmhmm.out)), 'catfile()';

ok my $parser =  Bio::Tools::Tmhmm->new(-file=>$infile), 'new()';

my @feat;
while ( my $feat = $parser->next_result ) {
  push @feat, $feat;
}

is @feat, 3, 'got 3 feat';

ok $feat[0]->seq_id,      'my_sequence_id';
ok $feat[0]->source_tag,  'TMHMM2.0';
ok $feat[0]->primary_tag, 'transmembrane';


my $raa_test_data = [
  [ 54, 76],
  [116, 138],
  [151, 173],
];

for (0..(scalar(@feat)-1)) {
	is $feat[$_]->start, $raa_test_data->[$_]->[0];
	is $feat[$_]->end, $raa_test_data->[$_]->[1];
}

