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

use Test;
use strict;

BEGIN { plan tests => 9 }
use Bio::Species;

ok(1);

my $sps = Bio::Species->new();
ok defined $sps;
my $msg = 'bug in either Species->classification() or elsewhere in Bio::Species';
$sps->classification('sapiens', 'Homo', 'Hominidae',
		     'Catarrhini', 'Primates', 'Eutheria', 'Mammalia',
		     'Vertebrata', 'Chordata', 'Metazoa', 'Eukaryota');
ok $sps->binomial(), 'Homo sapiens', $msg;

$sps->classification('sapiens', 'Homo', 'Hominidae',
		     'Catarrhini', 'Primates', 'Eutheria', 'Mammalia',
		     'Vertebrata', 'Chordata', 'Metazoa', 'Eukaryota');
$sps->sub_species('sapiensis');
ok $sps->binomial(), 'Homo sapiens', $msg;
ok $sps->binomial('FULL'), 'Homo sapiens sapiensis', $msg;
ok $sps->sub_species(), 'sapiensis', $msg;

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
