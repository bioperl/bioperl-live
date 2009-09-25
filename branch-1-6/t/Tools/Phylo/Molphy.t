# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 18,
			   -requires_module => 'IO::String');
	
	use_ok('Bio::Tools::Phylo::Molphy');
}

my $verbose = test_debug();

my $inmolphy = Bio::Tools::Phylo::Molphy->new(-file => test_input_file('lysozyme6.simple.protml'));
ok($inmolphy);
my $r = $inmolphy->next_result;
ok($r);
is($r->model, 'JTT');
is($r->search_space,50);
my @trees;
while( my $t = $r->next_tree ) { 
    push @trees, $t;
}
is(@trees,5);
 $inmolphy = Bio::Tools::Phylo::Molphy->new(-file => test_input_file('lysozyme6.protml'));
ok($inmolphy);
$r = $inmolphy->next_result;
is($r->model, 'JTT');
is($r->search_space,50);
@trees = ();
while( my $t = $r->next_tree ) { 
    push @trees, $t;
}
is(@trees,5);

is($trees[0]->score, -1047.8);
is($trees[-1]->id, 9);

my $tpm = $r->transition_probability_matrix;
is($tpm->{'Val'}->{'Val'}, -122884);
is($tpm->{'Ala'}->{'Arg'}, 2710);

my $sub_mat = $r->substitution_matrix;
is($sub_mat->{'Val'}->{'Tyr'}, 50);
is($sub_mat->{'Arg'}->{'Ile'}, 72);
is($sub_mat->{'Met'}->{'Met'}, '');

my %fmat = $r->residue_frequencies();
is($fmat{'D'}->[0], 0.052);
