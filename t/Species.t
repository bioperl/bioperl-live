## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
# Test script for Bio::Species.pm
# Hilmar Lapp <Hilmar.Lapp@pharma.novartis.com>, <hlapp@gmx.net>
# Fairly rudimentary
# Header code for this test was borrowed from Blast.t

## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------

BEGIN {	
    ## We start with some black magic to print on failure.
    $| = 1; print "1..5\n"; 
    # $^W = 0;  # I'm not sure what this means.
}
END {
    print "not ok 1\n" unless $loaded;
}

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

# use lib '../';
use Bio::Species;

$loaded = 1;

sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

my $sps = Bio::Species->new();
my $msg = 'bug in either Species->classification() or elsewhere in Bio::Species';
$sps->classification('sapiens', 'Homo', 'Hominidae',
		     'Catarrhini', 'Primates', 'Eutheria', 'Mammalia',
		     'Vertebrata', 'Chordata', 'Metazoa', 'Eukaryota');
test 1, $sps->binomial() eq 'Homo sapiens', $msg;

$sps->classification('sapiens', 'Homo', 'Hominidae',
		     'Catarrhini', 'Primates', 'Eutheria', 'Mammalia',
		     'Vertebrata', 'Chordata', 'Metazoa', 'Eukaryota');
$sps->sub_species('sapiensis');
test 2, $sps->binomial() eq 'Homo sapiens', $msg;
test 3, $sps->binomial('FULL') eq 'Homo sapiens sapiensis', $msg;
test 4, $sps->sub_species() eq 'sapiensis', $msg;

$sps->classification(qw( sapiens Homo Hominidae
			 Catarrhini Primates Eutheria Mammalia Vertebrata
			 Chordata Metazoa Eukaryota));
test 5, $sps->binomial() eq 'Homo sapiens', $msg;



