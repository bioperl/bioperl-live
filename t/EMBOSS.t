# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN {
    use vars qw($NTESTS);
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    $NTESTS = 4;
    plan tests => $NTESTS }

use Bio::Factory::EMBOSS;
use Bio::Root::IO;

my $compseqoutfile = '/tmp/dna1.4.compseq';
END { 
    foreach ( $Test::ntest..$NTESTS ) { skip(1,1,"EMBOSS not installed locally");}
    unlink($compseqoutfile) }
my $verbose = $ENV{'BIOPERLDEBUG'} || -1;
ok(1);

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my $factory = new Bio::Factory::EMBOSS(-verbose => $verbose);

ok($factory);
my $compseqapp = $factory->program('compseq');
if( ! $compseqapp ) { 
    # no EMBOSS installed
    exit();
}

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

		       
