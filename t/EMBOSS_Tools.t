# -*-Perl-*-
# $Id$
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
    plan tests => 12;
}

use Bio::Tools::EMBOSS::Palindrome;
use Bio::Tools::GFF;
ok(1);

my $parser = new Bio::Tools::EMBOSS::Palindrome
    (-file => Bio::Root::IO->catfile(qw( t data humts1.pal)));

my $seq = $parser->next_seq;
ok($seq);
ok($seq->display_id, 'HUMTS1');
ok($seq->length, 18596);
my @features = $seq->get_SeqFeatures();
ok(scalar @features, 23);

ok($features[0]->feature1->start, 126);
ok($features[0]->feature1->end, 142);
ok($features[0]->feature1->strand, 1);
ok($features[0]->feature1->seq_id, 'HUMTS1');


ok($features[0]->feature2->start, 201);
ok($features[0]->feature2->end, 217);
ok($features[0]->feature2->strand, -1);

if( $DEBUG ) {
    my $out = new Bio::Tools::GFF(-gff_version => 2);
    $out->write_feature($features[0]);
}
