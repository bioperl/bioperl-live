# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

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

    plan tests => 4 }

use Bio::Factory::EMBOSS;
use Bio::Root::IO;

my $compseqoutfile = '/tmp/dna1.4.compseq';
END { unlink($compseqoutfile) }
my $verbose = $ENV{'BIOPERLDEBUG'};
ok(1);

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my $factory = new Bio::Factory::EMBOSS(-verbose => $verbose);

ok($factory);
my $compseqapp = $factory->program('compseq');
ok($compseqapp);
my %input = ( '-word' => 4,
	      '-sequence' => Bio::Root::IO->catfile('t',
						   'data',
						   'dna1.fa'),
	      '-outfile' => $compseqoutfile);
$compseqapp->run(\%input);
ok(-e $compseqoutfile);
#if( open(IN, $compseqoutfile) ) {
#    while(<IN>) { print }
#}

		       
