# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
# Test script for Bio::Species.pm
# Hilmar Lapp <Hilmar.Lapp@pharma.novartis.com>, <hlapp@gmx.net>
# Fairly rudimentary
# Header code for this test was borrowed from Blast.t

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

    plan tests => 21;
}

use Bio::Species;
use Bio::Taxonomy::Node;

ok(1);

my $sps = Bio::Species->new();
ok defined $sps;
my $msg = 'bug in either Species->classification() or elsewhere in Bio::Species';
$sps->classification(qw( sapiens Homo Hominidae
			 Catarrhini Primates Eutheria Mammalia Vertebrata
			 Chordata Metazoa Eukaryota));

ok( $sps->binomial(), 'Homo sapiens');

$sps->sub_species('sapiensis');
ok($sps->binomial(), 'Homo sapiens');
ok($sps->binomial('FULL'), 'Homo sapiens sapiensis');
ok($sps->sub_species(), 'sapiensis');

$sps->classification(qw( sapiens Homo Hominidae
			 Catarrhini Primates Eutheria Mammalia Vertebrata
			 Chordata Metazoa Eukaryota));
ok $sps->binomial(), 'Homo sapiens', $msg;


# test cmd line initializtion
my $species = new Bio::Species( -classification => 
				[ qw( sapiens Homo Hominidae
				      Catarrhini Primates Eutheria 
				      Mammalia Vertebrata
				      Chordata Metazoa Eukaryota) ] );
ok( $species);
ok $species->binomial(), 'Homo sapiens';

# test if Taxonomy::Node can masquerade as a Bio::Species object
my $node = Bio::Taxonomy::Node->new(-rank => 'species');
ok( $node );
$node->classification(qw( sapiens Homo Hominidae
			 Catarrhini Primates Eutheria Mammalia Vertebrata
			 Chordata Metazoa Eukaryota));

ok($node->binomial(),$sps->binomial);

$node->sub_species('sapiensis');
ok($node->binomial(), 'Homo sapiens');
ok($node->binomial('FULL'), 'Homo sapiens sapiensis');
ok($node->sub_species(), 'sapiensis');

# test object level initializtion
my $newnode = new Bio::Taxonomy::Node( 
				       -rank           => 'species',
				       -organelle      => 'mito',
				       -sub_species    => 'sapiensis',
				       -variant        => 'varianthere',
				       -ncbi_taxid     => 9606,
				       -common_name    => 'human',
				       -classification => 
    [ qw( sapiens Homo Hominidae Catarrhini Primates Eutheria 
	  Mammalia Vertebrata Chordata Metazoa Eukaryota) ] );
ok( $newnode );
ok( $newnode->binomial(), 'Homo sapiens');
ok( $newnode->binomial('full'), 'Homo sapiens sapiensis');
ok( $newnode->sub_species, 'sapiensis');
ok( $newnode->organelle, 'mito');
ok( $newnode->variant, 'varianthere');
ok( $newnode->ncbi_taxid, 9606);
