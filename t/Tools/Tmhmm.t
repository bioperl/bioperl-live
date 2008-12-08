# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 12);
	
	use_ok('Bio::Tools::Tmhmm');
}

my $infile = test_input_file('tmhmm.out');

ok my $parser = Bio::Tools::Tmhmm->new(-file=>$infile), 'new()';

my @feat;
while ( my $feat = $parser->next_result ) {
  push @feat, $feat;
}

is @feat, 3, 'got 3 feat';

is $feat[0]->seq_id,      'my_sequence_id';
is $feat[0]->source_tag,  'TMHMM2.0';
is $feat[0]->primary_tag, 'transmembrane';


my $raa_test_data = [
  [ 54, 76],
  [116, 138],
  [151, 173],
];

for (0..(scalar(@feat)-1)) {
	is $feat[$_]->start, $raa_test_data->[$_]->[0];
	is $feat[$_]->end, $raa_test_data->[$_]->[1];
}
