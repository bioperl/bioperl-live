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

    $NUMTESTS = 72;
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

# AAML parsing - Empirical model
$inpaml = new Bio::Tools::Phylo::PAML(-file => Bio::Root::IO->catfile
				      (qw(t data aaml.mlc)));

ok($inpaml);
$result = $inpaml->next_result;
ok($result);
my @trees = $result->get_trees;
ok(@trees, 1);
ok($trees[0]->score, -1042.768973);

ok((scalar grep { $_->is_Leaf } $trees[0]->get_nodes), $result->get_seqs);

my $aadistmat = $result->get_AADistMatrix();
ok($aadistmat);
ok($aadistmat->get_entry('Cow', 'Horse'), 0.5462);
ok($aadistmat->get_entry('Baboon', 'Langur'), 0.1077);

my %aafreq = %{$result->get_AAFreqs()};
ok(%aafreq);
ok($aafreq{'Human'}->{'N'}, 0.0769);
ok($aafreq{'Human'}->{'R'}, 0.1077);

my %ratfreqs = %{$result->get_AAFreqs('Rat')};
ok($ratfreqs{'R'},0.0923);
ok($ratfreqs{'F'},0.0154);
my %avgfreqs = %{$result->get_AAFreqs('Average')};
ok($avgfreqs{'Q'},0.0411);

ok($result->get_AAFreqs('Average')->{'I'},0.0424);

my $patterns = $result->patterns;
my @pat = @{$patterns->{'-patterns'}};
ok(scalar @pat, 98);
ok($patterns->{'-ns'}, 6);
ok($patterns->{'-ls'}, 130);

ok(($result->get_stat_names)[0], 'constant_sites');
ok($result->get_stat('constant_sites'), 46);
ok($result->get_stat('constant_sites_percentage'), 35.38);

# AAML parsing - pairwise model
$inpaml = new Bio::Tools::Phylo::PAML(-file => Bio::Root::IO->catfile
				      (qw(t data aaml_pairwise.mlc)));

ok($inpaml);
$result = $inpaml->next_result;
ok($result);
ok($result->get_stat('loglikelihood'),-1189.106658);
ok($result->get_stat('constant_sites'), 170);
ok($result->get_stat('constant_sites_percentage'), 59.65);

ok($result->get_AAFreqs('Average')->{'R'},0.0211);
ok($result->get_AAFreqs('rabbit')->{'L'},0.1228);

$aadistmat = $result->get_AADistMatrix();
ok($aadistmat);
ok($aadistmat->get_entry('rabbit', 'marsupial'), 0.2877);
ok($aadistmat->get_entry('human', 'goat-cow'), 0.1439);

$aadistmat = $result->get_AAMLDistMatrix();
ok($aadistmat);
ok($aadistmat->get_entry('rabbit', 'marsupial'), 0.3392);
ok($aadistmat->get_entry('human', 'goat-cow'), 0.1551);

my @seqs = $result->get_seqs;
ok($seqs[0]->display_id, 'human');

# YN00 parsing, pairwise Ka/Ks from Yang & Nielsen 2000
$inpaml = new Bio::Tools::Phylo::PAML(-file => Bio::Root::IO->catfile
				      (qw(t data yn00.mlc)));

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


