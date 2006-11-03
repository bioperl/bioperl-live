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

    $NUMTESTS = 190;
    plan tests => $NUMTESTS;
    eval { require IO::String; 
	   require Bio::Tools::Phylo::PAML;}; 
    if( $@ ) {
	print STDERR "no IO::String installed\n"; 
	$error = 1;
    }
}

END {
	foreach ( $Test::ntest .. $NUMTESTS ) {
		skip("Unable to run all of the PAML tests",1);
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
ok($result->model, 'several dN/dS ratios for branches');
ok($result->version, qr'3\.12');
my $MLmat = $result->get_MLmatrix;
my $NGmat = $result->get_NGmatrix;

ok($NGmat->[0]->[1]->{'omega'}, 0.2507);
ok($NGmat->[0]->[1]->{'dN'}, 0.0863);
ok($NGmat->[0]->[1]->{'dS'}, 0.3443);

ok($NGmat->[2]->[3]->{'omega'}, 0.2178);
ok($NGmat->[2]->[3]->{'dN'}, 0.1348);
ok($NGmat->[2]->[3]->{'dS'}, 0.6187);

ok($MLmat->[0]->[1]->{'omega'}, 0.19479);
ok($MLmat->[0]->[1]->{'dN'}, 0.0839);
ok($MLmat->[0]->[1]->{'dS'}, 0.4309);
ok($MLmat->[0]->[1]->{'lnL'}, -1508.607268);
ok($MLmat->[0]->[1]->{'t'}, 0.47825);
ok($MLmat->[0]->[1]->{'kappa'}, 2.29137);

ok($MLmat->[2]->[3]->{'omega'}, 0.16114);
ok($MLmat->[2]->[3]->{'dN'}, 0.1306);
ok($MLmat->[2]->[3]->{'dS'}, 0.8105);
ok($MLmat->[2]->[3]->{'lnL'},-1666.440696);
ok($MLmat->[2]->[3]->{'t'}, 0.85281);
ok($MLmat->[2]->[3]->{'kappa'}, 2.21652);

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
ok($result->model, 'Empirical (wag.dat)');
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

ok((sort $result->get_stat_names)[0], 'constant_sites');
ok($result->get_stat('constant_sites'), 46);
ok($result->get_stat('constant_sites_percentage'), 35.38);

# AAML parsing - pairwise model
$inpaml = new Bio::Tools::Phylo::PAML(-file => Bio::Root::IO->catfile
				      (qw(t data aaml_pairwise.mlc)));

ok($inpaml);
$result = $inpaml->next_result;
ok($result);
ok($result->model, 'Empirical_F (wag.dat)');
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

# codeml NSSites parsing

$inpaml = new Bio::Tools::Phylo::PAML
    (-file => Bio::Root::IO->catfile(qw(t data codeml_nssites.mlc)));

ok($inpaml);
$result = $inpaml->next_result;

ok($result);
ok($result->model, 'One dN/dS ratio dGamma (ncatG=11)');
ok($result->version, 'paml 3.13, August 2002');
$NGmat = $result->get_NGmatrix;
ok($NGmat);

ok($NGmat->[0]->[1]->{'omega'}, 0.2782);
ok($NGmat->[0]->[1]->{'dN'}, 0.0133);
ok($NGmat->[0]->[1]->{'dS'}, 0.0478);
ok($NGmat->[1]->[2]->{'omega'}, 1.1055);
ok($NGmat->[1]->[2]->{'dN'}, 0.0742);
ok($NGmat->[1]->[2]->{'dS'}, 0.0671);
          # this is
          #   model num  description
          #   kappa   log-likelihood tree length time used
          #   shape   alpha/gamma r          f
my @tstr = ([qw(0 one-ratio 0
		4.54006 -906.017440    0.55764
		)],
	    [qw(1 neutral 2
		4.29790 -902.503869    0.56529
		)],
	    [qw(2 selection 3 
		5.12250 -900.076500    0.6032
		)],
	     );
my $iter = 0;
my $lastmodel;
foreach my $model ( $result->get_NSSite_results ) {    
    my $i = 0;
    my $r = shift @tstr;
    ok($model->model_num, $r->[$i++]);
    ok($model->model_description, qr/$r->[$i++]/);
    ok($model->num_site_classes,$r->[$i++]);
    my $tree = $model->next_tree;
    ok($model->kappa, $r->[$i++]);
    ok($model->likelihood,$r->[$i]);
    ok($tree->score, $r->[$i++]);
    ok($tree->total_branch_length, $r->[$i++]);
    if( $iter == 0 ) {
	my $params = $model->shape_params;
	ok($params->{'shape'}, 'alpha');
	ok($params->{'gamma'},   '0.50000');
	ok($params->{'r'}->[0], '1.00000');
	ok($params->{'f'}->[0], '1.00000');
    } elsif( $iter == 2 ) {
	my $class = $model->dnds_site_classes;
	ok($class->{'p'}->[0], '0.38160');
	ok($class->{'w'}->[1], '1.00000');
    }
    $iter++;
    $lastmodel = $model;
}

my ($firstsite) = $lastmodel->get_pos_selected_sites;
ok($firstsite->[0], 15);
ok($firstsite->[1], 'L');
ok($firstsite->[2], 0.6588);

# codeml NSSites parsing
# for M0 model

my $codeml_m0 = new Bio::Tools::Phylo::PAML
    (-file => Bio::Root::IO->catfile(qw/t data M0.mlc/));
ok($codeml_m0);
my $result_m0 = $codeml_m0->next_result;
my ($nssite_m0,$nssite_m1) = $result_m0->get_NSSite_results;
ok($nssite_m0->num_site_classes,1);
my $class_m0 = $nssite_m0->dnds_site_classes;
ok($class_m0->{q/p/}->[0],q/1.00000/);
ok($class_m0->{q/w/}->[0],0.09213);

@trees= $nssite_m0->get_trees;
ok (@trees , 1 ); 
# model 0
ok($trees[0]->score, -30.819156);

@trees= $nssite_m1->get_trees;
ok($trees[0]->score, -30.819157);

# test BASEML
# pairwise first

my $baseml_p = Bio::Tools::Phylo::PAML->new
    (-file => Bio::Root::IO->catfile(qw(t data baseml.pairwise)));
ok($baseml_p);
my $baseml = $baseml_p->next_result;
my @b_seqs =  $baseml->get_seqs;
ok($b_seqs[0]->seq, 'GTAGAGTACTTT');
ok($b_seqs[1]->seq, 'GTAAGAGACGAT');

my @otus = map { $_->display_id } @b_seqs;
ok(scalar @otus, 3);
my $ntfreq = $baseml->get_NTFreqs;
ok($ntfreq);
ok($ntfreq->{$otus[0]}->{'A'}, '0.3333');
ok($ntfreq->{$otus[1]}->{'G'}, '0.2105');
my $kappaM = $baseml->get_KappaMatrix;
ok($kappaM);
ok($kappaM->get_entry($otus[1],$otus[0]), '0.3240');
ok($kappaM->get_entry($otus[0],$otus[1]), 
   $kappaM->get_entry($otus[1],$otus[0]));
ok($kappaM->get_entry($otus[1],$otus[2]), '0.1343');
my $alphaM = $baseml->get_AlphaMatrix;
ok($alphaM);
ok($alphaM->get_entry($otus[1],$otus[0]), '9.3595');
ok($alphaM->get_entry($otus[0],$otus[1]), 
   $alphaM->get_entry($otus[1],$otus[0]));
ok($alphaM->get_entry($otus[1],$otus[2]), '1.1101');
ok($alphaM->get_entry($otus[0],$otus[2]), '33.1197');

# codeml NSSites parsing
# for only 1 model

my $codeml_single = new Bio::Tools::Phylo::PAML
    (-file => Bio::Root::IO->catfile(qw/t data singleNSsite.mlc/));
ok($codeml_single);
my $result_single = $codeml_single->next_result;
my ($nssite_single) = $result_single->get_NSSite_results;
ok($nssite_single->num_site_classes,q/3/);
ok($nssite_single->kappa, q/5.28487/);
ok($nssite_single->likelihood,q/-30.819156/);

ok($baseml->get_stat('loglikelihood'),-110.532715);
ok($baseml->get_stat('constant_sites'),46);
ok($baseml->get_stat('constant_sites_percentage'),'80.70');
ok($baseml->model,'HKY85 dGamma (ncatG=5)');

# user trees
$baseml_p = Bio::Tools::Phylo::PAML->new
    (-file => Bio::Root::IO->catfile(qw(t data baseml.usertree)));
$baseml = $baseml_p->next_result;

@trees = $baseml->get_trees;
ok(@trees, 1);
ok($trees[0]->score, -129.328757);

# codeml NSSites parsing
# for branch site model/clade model

my $codeml_bs = new Bio::Tools::Phylo::PAML
    (-file => Bio::Root::IO->catfile(qw/t data branchSite.mlc/));
ok($codeml_bs);
my $result_bs = $codeml_bs->next_result;
my ($nssite_bs) = $result_bs->get_NSSite_results;
ok($nssite_bs->num_site_classes,q/4/);
my $class_bs = $nssite_bs->dnds_site_classes;
ok($class_bs->{q/p/}->[1],q/0.65968/);
ok($class_bs->{q/w/}->[1]->{q/background/},q/0.00000/);
ok($class_bs->{q/w/}->[2]->{q/foreground/},q/999.00000/);

# Let's parse the RST file

my $paml = Bio::Tools::Phylo::PAML->new
    (-file => Bio::Root::IO->catfile(qw(t data codeml_lysozyme mlc)),
     -dir  => Bio::Root::IO->catfile(qw(t data codeml_lysozyme)));

$result = $paml->next_result;

my ($rst) = grep {$_->id eq 'node#8'} $result->get_rst_seqs;
ok($rst);
ok($rst->seq, join('',qw(
AAGGTCTTTGAAAGGTGTGAGTTGGCCAGAACTCTGAAAAGATTGGGACTGGATGGCTAC
AGGGGAATCAGCCTAGCAAACTGGATGTGTTTGGCCAAATGGGAGAGTGATTATAACACA
CGAGCTACAAACTACAATCCTGGAGACCAAAGCACTGATTATGGGATATTTCAGATCAAT
AGCCACTACTGGTGTAATAATGGCAAAACCCCAGGAGCAGTTAATGCCTGTCATATATCC
TGCAATGCTTTGCTGCAAGATAACATCGCTGATGCTGTAGCTTGTGCAAAGAGGGTTGTC
CGTGATCCACAAGGCATTAGAGCATGGGTGGCATGGAGAAATCATTGTCAAAACAGAGAT
GTCAGTCAGTATGTTCAAGGTTGTGGAGTG)),
   'node#8 reconstructed seq');

my ($first_tree) = $result->get_rst_trees;
my ($node) = $first_tree->find_node(-id => '5_Mmu_rhesus');
my @changes = $node->get_tag_values('changes');
my ($site) = grep { $_->{'site'} == 94 } @changes;
ok($site->{'anc_aa'}, 'A','ancestral AA');
ok($site->{'anc_prob'}, '0.947','ancestral AA prob');
ok($site->{'derived_aa'}, 'T','derived AA');

($node) = $first_tree->find_node(-id => '12');
@changes = $node->get_tag_values('changes');
($site) = grep { $_->{'site'} == 88 } @changes;
ok($site->{'anc_aa'}, 'N','ancestral AA');
ok($site->{'anc_prob'}, '0.993','ancestral AA prob');
ok($site->{'derived_aa'}, 'D','derived AA');
ok($site->{'derived_prob'}, '0.998','derived AA prob');

my $persite = $result->get_rst_persite;
# minus 1 because we have shifted so that array index matches site number
# there are 130 sites in this seq file
ok(scalar @$persite -1, $result->patterns->{'-ls'}); 
# let's score site 1
$site = $persite->[2]; 
# so site 2, node 2 (extant)
ok($site->[2]->{'codon'}, 'GTC');
ok($site->[2]->{'aa'}, 'V');
# site 2, node 3
ok($site->[3]->{'codon'}, 'ATC');
ok($site->[3]->{'aa'}, 'I');

# ancestral node 9
ok($site->[9]->{'codon'}, 'GTC');
ok($site->[9]->{'aa'},    'V');
ok($site->[9]->{'prob'},  '1.000');
ok($site->[9]->{'Yang95_aa'},'V');
ok($site->[9]->{'Yang95_aa_prob'},'1.000');

# ancestral node 10
ok($site->[10]->{'codon'}, 'ATC');
ok($site->[10]->{'aa'},    'I');
ok($site->[10]->{'prob'},  '0.992');
ok($site->[10]->{'Yang95_aa'},'I');
ok($site->[10]->{'Yang95_aa_prob'},'0.992');


## PAML 3.15
$paml = Bio::Tools::Phylo::PAML->new(-file => Bio::Root::IO->catfile(qw(t data codeml315.mlc)) );
$result = $paml->next_result;

ok($result->model, 'One dN/dS ratio');
ok($result->version, qr'3\.15');
$MLmat = $result->get_MLmatrix;
$NGmat = $result->get_NGmatrix;

ok($NGmat->[0]->[1]->{'omega'}, 0.2264);
ok($NGmat->[0]->[1]->{'dN'}, 0.0186);
ok($NGmat->[0]->[1]->{'dS'}, 0.0821);

ok($MLmat->[0]->[1]->{'omega'}, 0.32693);
ok($MLmat->[0]->[1]->{'dN'}, '0.0210');
ok($MLmat->[0]->[1]->{'dS'}, 0.0644);
