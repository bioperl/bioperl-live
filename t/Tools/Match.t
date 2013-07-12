# -*-Perl-*- Test Harness script for Bioperl
# $Id: Match.t,v 1.15 2007/06/27 10:16:38 sendu Exp $

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 38,
			   -requires_module => 'URI::Escape');
	
    use_ok('Bio::Tools::Match');
}

ok my $parser = Bio::Tools::Match->new(-file => test_input_file('match.output'));

my $count = 0;
my @expected = ([qw(4338 4349 1.000 0.945 V$MYOD_01)],
                [qw(7390 7401 1.000 0.932 V$MYOD_01)],
                [qw(8503 8514 1.000 0.941 V$MYOD_01)],
                [qw(8767 8778 1.000 0.937 V$MYOD_01)],
                [qw(33 47 0.693 0.779 V$E47_01)]);
while (my $feat = $parser->next_result) {
    $count++;
    my @exp = @{shift(@expected)};
    
    isa_ok $feat, 'Bio::SeqFeature::Generic';
    is $feat->source_tag, 'transfac_match', 'correct source';
    is $feat->start, shift(@exp), 'feature start correct';
    is $feat->end, shift(@exp), 'feature end correct';
    
    my $core_score = $feat->score;
    my $matrix_score = ($feat->annotation->get_Annotations('matrix_score'))[0]->value;
    my $matrix_id = ($feat->annotation->get_Annotations('matrix_id'))[0]->value;
    
    is $core_score, shift(@exp), 'feature core score correct';
    is $matrix_score, shift(@exp), 'feature matrix score correct';
    is $matrix_id, shift(@exp), 'feature matrix id correct';
    
    last if $count == 5;
}

is $count, 5, "correct number of results managed to get tested";
