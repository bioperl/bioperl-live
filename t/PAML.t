# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $error);


BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    $error = 0;
    if( $@ ) {
	use lib 't';
    }
    use Test;

    $NUMTESTS = 38;
    plan tests => $NUMTESTS;
    eval { require IO::String; 
	   require Bio::Tools::Phylo::PAML;}; 
    if( $@ ) { print STDERR "no IO string installed\n"; 
	$error = 1;
	}
}

END { 
    foreach ( $Test::ntest .. $NUMTESTS ) {
	skip("unable to run all of the PAML tests",1);
    }
}


exit(0) if( $error );

my $testnum;
my $verbose = 0;

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 
use Bio::Root::IO;

my $inpaml = new Bio::Tools::Phylo::PAML(-file => Bio::Root::IO->catfile(qw(t data codeml.mlc)));
ok($inpaml);
my $result = $inpaml->next_result;
ok($result);
ok($result->version, qr'3\.12');
my $MLmat = $result->get_MLmatrix;
my $NGmat = $result->get_NGmatrix;

ok($NGmat->[0]->[1]->{'omega'}, 0.2507);
ok($NGmat->[0]->[1]->{'dN'}, 0.0863);
ok($NGmat->[0]->[1]->{'dS'}, 0.3443);
ok($NGmat->[2]->[3]->{'omega'}, 0.2178);
ok($NGmat->[2]->[3]->{'dN'}, 0.1348);
ok($NGmat->[2]->[3]->{'dS'}, 0.6187);

ok($MLmat->[0]->[1]->{'omega'}, 0.1948);
ok($MLmat->[0]->[1]->{'dN'}, 0.0839);
ok($MLmat->[0]->[1]->{'dS'}, 0.4309);
ok($MLmat->[0]->[1]->{'lnL'}, -1508.607268);
ok($MLmat->[2]->[3]->{'omega'}, 0.1611);
ok($MLmat->[2]->[3]->{'dN'}, 0.1306);
ok($MLmat->[2]->[3]->{'dS'}, 0.8105);
ok($MLmat->[2]->[3]->{'lnL'},-1666.440696);

my @codonposfreq = $result->get_codon_pos_basefreq();
ok($codonposfreq[0]->{'A'}, 0.23579);
ok($codonposfreq[0]->{'T'}, 0.14737);
ok($codonposfreq[1]->{'C'}, 0.25123);
ok($codonposfreq[2]->{'G'}, 0.32842);

$inpaml = new Bio::Tools::Phylo::PAML(-file => Bio::Root::IO->catfile(qw(t data aaml.mlc)));

ok($inpaml);
#$result = $inpaml->next_result;
#ok($result);

$inpaml = new Bio::Tools::Phylo::PAML(-file => Bio::Root::IO->catfile(qw(t data yn00.mlc)));

ok($inpaml);
$result = $inpaml->next_result;

ok($result);
$MLmat = $result->get_MLmatrix;
$NGmat = $result->get_NGmatrix;


ok($NGmat->[0]->[1]->{'omega'}, 0.251);
ok($NGmat->[0]->[1]->{'dN'}, 0.0863);
ok($NGmat->[0]->[1]->{'dS'}, 0.3443);
ok($NGmat->[2]->[3]->{'omega'}, 0.218);
ok($NGmat->[2]->[3]->{'dN'}, 0.1348);
ok($NGmat->[2]->[3]->{'dS'}, 0.6187);

ok($MLmat->[0]->[1]->{'omega'}, 0.1625);
ok($MLmat->[0]->[1]->{'dN'}, 0.0818);
ok($MLmat->[0]->[1]->{'dS'}, 0.5031);
ok($MLmat->[2]->[3]->{'omega'}, 0.1262);
ok($MLmat->[2]->[3]->{'dN'}, 0.1298);
ok($MLmat->[2]->[3]->{'dN_SE'}, 0.0149);
ok($MLmat->[2]->[3]->{'dS'}, 1.0286);
ok($MLmat->[2]->[3]->{'dS_SE'}, 0.2614);

