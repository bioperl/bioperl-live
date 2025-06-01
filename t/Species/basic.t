#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 8;

use_ok('Bio::Species');

my $sps = Bio::Species->new();
$sps->classification(qw(sapiens Homo Hominidae Catarrhini Primates Eutheria Mammalia Vertebrata Chordata Metazoa Eukaryota));

is $sps->binomial, 'Homo sapiens';

$sps->sub_species('sapiensis');
is $sps->binomial, 'Homo sapiens';
is $sps->binomial('FULL'), 'Homo sapiens sapiensis';

is $sps->sub_species, 'sapiensis';

my $species = Bio::Species->new(
    -classification => [qw(sapiens Homo Hominidae Catarrhini Primates Eutheria Mammalia Vertebrata Chordata Metazoa Eukaryota)],
    -common_name    => 'human'
);
is $species->binomial, 'Homo sapiens';
is $species->genus, 'Homo';

