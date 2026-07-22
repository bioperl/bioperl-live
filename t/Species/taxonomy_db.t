#!/usr/bin/env perl

use strict;
use warnings;

use Bio::Root::Test;
use Bio::Species;
use Bio::DB::Taxonomy;

test_begin(-tests => 5,
           -requires_module     => 'Bio::DB::Taxonomy::entrez',
           -requires_networking => 1);

my $species = Bio::Species->new(-id => 51351);
my $taxdb   = Bio::DB::Taxonomy->new(-source => 'entrez');

eval { $species->db_handle($taxdb); };
skip("Unable to connect to Entrez database; no network or server busy?", 5) if $@;

is $species->binomial, 'Brassica rapa subsp.';
is $species->binomial('FULL'), 'Brassica rapa subsp. pekinensis';
is $species->genus, 'Brassica';
is $species->species, 'rapa subsp.';
is $species->sub_species, 'pekinensis';

