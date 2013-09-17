# -*-Perl-*- Test Harness script for Bioperl
# $Id: gerp.t,v 1.15 2007/06/27 10:16:38 sendu Exp $

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 33,
			   -requires_module => 'URI::Escape');
	
    use_ok('Bio::Tools::Phylo::Gerp');
}

ok my $parser = Bio::Tools::Phylo::Gerp->new(-file => test_input_file('ENr111.mfa.example.elems'));

my $count = 0;
my @expected = ([qw(334180 334352 449 1.03744e-165)],
                [qw(337735 337915 458.2 5.02405e-164)],
                [qw(262604 262861 473.1 3.64789e-117)],
                [qw(285427 285608 386.1 8.42494e-113)],
                [qw(309563 309744 383.6 2.88895e-111)]);
while (my $feat = $parser->next_result) {
    $count++;
    my @exp = @{shift(@expected)};
    
    isa_ok $feat, 'Bio::SeqFeature::Generic';
    is $feat->source_tag, 'GERP', 'correct source';
    is $feat->start, shift(@exp), 'feature start correct';
    is $feat->end, shift(@exp), 'feature end correct';
    is $feat->score, shift(@exp), 'feature score correct';
    my ($p_value) = $feat->annotation->get_Annotations('pvalue');
    is ref $p_value ? $p_value->value : $p_value, shift(@exp), 'feature pvalue correct';
}

is $count, 5, "correct number of results parsed out";
