#-*-Perl-*-
## Bioperl Test Harness Script for Modules

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;

BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan test => 16;
}

use Bio::Seq;
use Bio::Tools::GFF;
use Bio::SeqFeatureI;
use Bio::SeqFeature::Generic;
my $feat = new Bio::SeqFeature::Generic( -start => 10, -end => 100,
				-strand => -1, -primary => 'repeat',
				-source => 'repeatmasker',
				-score  => 1000,
				-tag    => {
				    new => 1,
				    author => 'someone',
				    sillytag => 'this is silly!;breakfast' } );
ok($feat);
my $gff1out = Bio::Tools::GFF->new(-gff_version => 1, -file => ">out1.gff");
ok($gff1out);
my $gff2out = Bio::Tools::GFF->new(-gff_version => 2, -file => ">out2.gff");
ok($gff2out);

$gff1out->write_feature($feat);
$gff2out->write_feature($feat);

$gff1out->close();
$gff2out->close();

my $gff1in = Bio::Tools::GFF->new(-gff_version => 1,  -file => "out1.gff");
ok($gff1in);
my $gff2in = Bio::Tools::GFF->new(-gff_version => 2, -file => "out2.gff");
ok($gff2in);

my $feat1 = $gff1in->next_feature();
ok($feat1);
ok($feat1->start, $feat->start);
ok($feat1->end, $feat->end);
ok($feat1->primary_tag, $feat->primary_tag);
ok($feat1->score, $feat->score);

my $feat2 = $gff2in->next_feature();
ok($feat2);
ok($feat2->start, $feat->start);
ok($feat2->end, $feat->end);
ok($feat2->primary_tag, $feat->primary_tag);
ok($feat2->score, $feat->score);
ok(($feat2->each_tag_value('sillytag'))[0], 'this is silly!;breakfast');
END {
    unlink("out1.gff", "out2.gff");
}
