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
    plan test => 1;
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
				    sillytag => 'this is silly!' } );

my $gff1out = Bio::Tools::GFF->new(-gff_version => 1);
my $gff2out = Bio::Tools::GFF->new(-gff_version => 2);

$gff1out->write_feature($feat);
$gff2out->write_feature($feat);

my $gff1in = Bio::Tools::GFF->new(-gff_version => 1);
my $gff2in = Bio::Tools::GFF->new(-gff_version => 2);

my $feat1 = $gff1in->next_feature();
my $feat2 = $gff2in->next_feature();
