# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 7);
	
    use_ok('Bio::Tools::Promoterwise');
}

my $file = test_input_file('promoterwise.out');
my $parser = Bio::Tools::Promoterwise->new(-file=>$file);
isa_ok $parser,'Bio::Tools::Promoterwise';
my @fp;
while (my $fp = $parser->next_result){
  push @fp,$fp;
}
my $first = $fp[0]->feature1;
my $second = $fp[0]->feature2;

my @sub = $first->sub_SeqFeature;
my @sub2 = $second->sub_SeqFeature;

is $sub[0]->start,4;
is $sub2[0]->start,29;
is $sub[0]->end,18;
is $sub2[0]->end,43;
is $sub[0]->score,1596.49
