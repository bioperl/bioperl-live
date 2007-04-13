# -*-Perl-*-
# $Id$
## Bioperl Test Harness Script for Modules
##

my $error;
use strict;
my $DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
BEGIN { 
    $error = 0; 
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;
    plan tests => 13;
	use_ok('Bio::Tools::EMBOSS::Palindrome');
	use_ok('Bio::Tools::GFF');
}

my $parser = new Bio::Tools::EMBOSS::Palindrome
    (-file => Bio::Root::IO->catfile(qw( t data humts1.pal)));

my $seq = $parser->next_seq;
ok($seq);
is($seq->display_id, 'HUMTS1');
is($seq->length, 18596);
my @features = $seq->get_SeqFeatures();
is(scalar @features, 23);

is($features[0]->feature1->start, 126);
is($features[0]->feature1->end, 142);
is($features[0]->feature1->strand, 1);
is($features[0]->feature1->seq_id, 'HUMTS1');


is($features[0]->feature2->start, 201);
is($features[0]->feature2->end, 217);
is($features[0]->feature2->strand, -1);

if( $DEBUG ) {
    my $out = new Bio::Tools::GFF(-gff_version => 2);
    $out->write_feature($features[0]);
}
