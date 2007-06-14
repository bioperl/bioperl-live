# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS);

BEGIN { 
    eval { require Test::More; };
    if( $@ ) {
		use lib 't/lib';
    }
    use Test::More;

    $NUMTESTS = 194;
    eval { 	require IO::String; }; 
    if( $@ ) {
		plan skip_all => "no IO::String installed\n";
	} else {
	    plan tests => $NUMTESTS;
    }
	use_ok('Bio::Tools::Phylo::PAML');
	use_ok('Bio::Root::IO');
}

my $testnum;
my $verbose = 0;

## End of black magic.

my $inpaml = Bio::Tools::Phylo::PAML->new(-file => Bio::Root::IO->catfile(qw(t data codeml.mlc)));
ok($inpaml);
my $result = $inpaml->next_result;
ok($result);
is($result->model, 'several dN/dS ratios for branches');
like($result->version, qr'3\.12');
my $MLmat = $result->get_MLmatrix;
my $NGmat = $result->get_NGmatrix;

is($NGmat->[0]->[1]->{'omega'}, 0.2507);
is($NGmat->[0]->[1]->{'dN'}, 0.0863);
is($NGmat->[0]->[1]->{'dS'}, 0.3443);

is($NGmat->[2]->[3]->{'omega'}, 0.2178);
is($NGmat->[2]->[3]->{'dN'}, 0.1348);
is($NGmat->[2]->[3]->{'dS'}, 0.6187);

is($MLmat->[0]->[1]->{'omega'}, 0.19479);
is($MLmat->[0]->[1]->{'dN'}, 0.0839);
is($MLmat->[0]->[1]->{'dS'}, 0.4309);
is($MLmat->[0]->[1]->{'lnL'}, -1508.607268);
is($MLmat->[0]->[1]->{'t'}, 0.47825);
is($MLmat->[0]->[1]->{'kappa'}, 2.29137);

is($MLmat->[2]->[3]->{'omega'}, 0.16114);
is($MLmat->[2]->[3]->{'dN'}, 0.1306);
is($MLmat->[2]->[3]->{'dS'}, 0.8105);
is($MLmat->[2]->[3]->{'lnL'},-1666.440696);
is($MLmat->[2]->[3]->{'t'}, 0.85281);
is($MLmat->[2]->[3]->{'kappa'}, 2.21652);

my @codonposfreq = $result->get_codon_pos_basefreq();
is($codonposfreq[0]->{'A'}, 0.23579);
is($codonposfreq[0]->{'T'}, 0.14737);
is($codonposfreq[1]->{'C'}, 0.25123);
is($codonposfreq[2]->{'G'}, 0.32842);

# AAML parsing - Empirical model
$inpaml = Bio::Tools::Phylo::PAML->new(-file => Bio::Root::IO->catfile
				      (qw(t data aaml.mlc)));

ok($inpaml);
$result = $inpaml->next_result;
ok($result);
is($result->model, 'Empirical (wag.dat)');
my @trees = $result->get_trees;
is(@trees, 1);
is($trees[0]->score, -1042.768973);

is((scalar grep { $_->is_Leaf } $trees[0]->get_nodes), $result->get_seqs);

my $aadistmat = $result->get_AADistMatrix();
ok($aadistmat);
is($aadistmat->get_entry('Cow', 'Horse'), 0.5462);
is($aadistmat->get_entry('Baboon', 'Langur'), 0.1077);

my %aafreq = %{$result->get_AAFreqs()};
ok(%aafreq);
is($aafreq{'Human'}->{'N'}, 0.0769);
is($aafreq{'Human'}->{'R'}, 0.1077);

my %ratfreqs = %{$result->get_AAFreqs('Rat')};
is($ratfreqs{'R'},0.0923);
is($ratfreqs{'F'},0.0154);
my %avgfreqs = %{$result->get_AAFreqs('Average')};
is($avgfreqs{'Q'},0.0411);

is($result->get_AAFreqs('Average')->{'I'},0.0424);

my $patterns = $result->patterns;
my @pat = @{$patterns->{'-patterns'}};
is(scalar @pat, 98);
is($patterns->{'-ns'}, 6);
is($patterns->{'-ls'}, 130);

is((sort $result->get_stat_names)[0], 'constant_sites');
is($result->get_stat('constant_sites'), 46);
is($result->get_stat('constant_sites_percentage'), 35.38);

# AAML parsing - pairwise model
$inpaml = Bio::Tools::Phylo::PAML->new(-file => Bio::Root::IO->catfile
				      (qw(t data aaml_pairwise.mlc)));

ok($inpaml);
$result = $inpaml->next_result;
ok($result);
is($result->model, 'Empirical_F (wag.dat)');
is($result->get_stat('loglikelihood'),-1189.106658);
is($result->get_stat('constant_sites'), 170);
is($result->get_stat('constant_sites_percentage'), 59.65);

is($result->get_AAFreqs('Average')->{'R'},0.0211);
is($result->get_AAFreqs('rabbit')->{'L'},0.1228);

$aadistmat = $result->get_AADistMatrix();
ok($aadistmat);
is($aadistmat->get_entry('rabbit', 'marsupial'), 0.2877);
is($aadistmat->get_entry('human', 'goat-cow'), 0.1439);

$aadistmat = $result->get_AAMLDistMatrix();
ok($aadistmat);
is($aadistmat->get_entry('rabbit', 'marsupial'), 0.3392);
is($aadistmat->get_entry('human', 'goat-cow'), 0.1551);

my @seqs = $result->get_seqs;
is($seqs[0]->display_id, 'human');

# YN00 parsing, pairwise Ka/Ks from Yang & Nielsen 2000
$inpaml = Bio::Tools::Phylo::PAML->new(-file => Bio::Root::IO->catfile
				      (qw(t data yn00.mlc)));

ok($inpaml);
$result = $inpaml->next_result;

ok($result);
$MLmat = $result->get_MLmatrix;
$NGmat = $result->get_NGmatrix;

is($NGmat->[0]->[1]->{'omega'}, 0.251);
is($NGmat->[0]->[1]->{'dN'}, 0.0863);
is($NGmat->[0]->[1]->{'dS'}, 0.3443);
is($NGmat->[2]->[3]->{'omega'}, 0.218);
is($NGmat->[2]->[3]->{'dN'}, 0.1348);
is($NGmat->[2]->[3]->{'dS'}, 0.6187);

is($MLmat->[0]->[1]->{'omega'}, 0.1625);
is($MLmat->[0]->[1]->{'dN'}, 0.0818);
is($MLmat->[0]->[1]->{'dS'}, 0.5031);
is($MLmat->[2]->[3]->{'omega'}, 0.1262);
is($MLmat->[2]->[3]->{'dN'}, 0.1298);
is($MLmat->[2]->[3]->{'dN_SE'}, 0.0149);
is($MLmat->[2]->[3]->{'dS'}, 1.0286);
is($MLmat->[2]->[3]->{'dS_SE'}, 0.2614);

# codeml NSSites parsing

$inpaml = Bio::Tools::Phylo::PAML->new
    (-file => Bio::Root::IO->catfile(qw(t data codeml_nssites.mlc)));

ok($inpaml);
$result = $inpaml->next_result;

ok($result);
is($result->model, 'One dN/dS ratio dGamma (ncatG=11)');
is($result->version, 'paml 3.13, August 2002');
$NGmat = $result->get_NGmatrix;
ok($NGmat);

is($NGmat->[0]->[1]->{'omega'}, 0.2782);
is($NGmat->[0]->[1]->{'dN'}, 0.0133);
is($NGmat->[0]->[1]->{'dS'}, 0.0478);
is($NGmat->[1]->[2]->{'omega'}, 1.1055);
is($NGmat->[1]->[2]->{'dN'}, 0.0742);
is($NGmat->[1]->[2]->{'dS'}, 0.0671);
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
    is($model->model_num, $r->[$i++]);
    like($model->model_description, qr/$r->[$i++]/);
    is($model->num_site_classes,$r->[$i++]);
    my $tree = $model->next_tree;
    is($model->kappa, $r->[$i++]);
    is($model->likelihood,$r->[$i]);
    is($tree->score, $r->[$i++]);
    is($tree->total_branch_length, $r->[$i++]);
    if( $iter == 0 ) {
	my $params = $model->shape_params;
	is($params->{'shape'}, 'alpha');
	is($params->{'gamma'},   '0.50000');
	is($params->{'r'}->[0], '1.00000');
	is($params->{'f'}->[0], '1.00000');
    } elsif( $iter == 2 ) {
	my $class = $model->dnds_site_classes;
	is($class->{'p'}->[0], '0.38160');
	is($class->{'w'}->[1], '1.00000');
    }
    $iter++;
    $lastmodel = $model;
}

my ($firstsite) = $lastmodel->get_pos_selected_sites;
is($firstsite->[0], 15);
is($firstsite->[1], 'L');
is($firstsite->[2], 0.6588);

# codeml NSSites parsing
# for M0 model

my $codeml_m0 = Bio::Tools::Phylo::PAML->new
    (-file => Bio::Root::IO->catfile(qw/t data M0.mlc/));
ok($codeml_m0);
my $result_m0 = $codeml_m0->next_result;
my ($nssite_m0,$nssite_m1) = $result_m0->get_NSSite_results;
is($nssite_m0->num_site_classes,1);
my $class_m0 = $nssite_m0->dnds_site_classes;
is($class_m0->{q/p/}->[0],q/1.00000/);
is($class_m0->{q/w/}->[0],0.09213);

is($nssite_m0->model_num, "0");
@trees= $nssite_m0->get_trees;
is (@trees , 1 ); 
# model 0
is($trees[0]->score, -30.819156);
is($nssite_m1->model_num, "1");
@trees= $nssite_m1->get_trees;
is($trees[0]->score, -30.819157);

# test BASEML
# pairwise first

my $baseml_p = Bio::Tools::Phylo::PAML->new
    (-file => Bio::Root::IO->catfile(qw(t data baseml.pairwise)));
ok($baseml_p);
my $baseml = $baseml_p->next_result;
my @b_seqs =  $baseml->get_seqs;
is($b_seqs[0]->seq, 'GTAGAGTACTTT');
is($b_seqs[1]->seq, 'GTAAGAGACGAT');

my @otus = map { $_->display_id } @b_seqs;
is(scalar @otus, 3);
my $ntfreq = $baseml->get_NTFreqs;
ok($ntfreq);
is($ntfreq->{$otus[0]}->{'A'}, '0.3333');
is($ntfreq->{$otus[1]}->{'G'}, '0.2105');
my $kappaM = $baseml->get_KappaMatrix;
ok($kappaM);
is($kappaM->get_entry($otus[1],$otus[0]), '0.3240');
is($kappaM->get_entry($otus[0],$otus[1]), 
   $kappaM->get_entry($otus[1],$otus[0]));
is($kappaM->get_entry($otus[1],$otus[2]), '0.1343');
my $alphaM = $baseml->get_AlphaMatrix;
ok($alphaM);
is($alphaM->get_entry($otus[1],$otus[0]), '9.3595');
is($alphaM->get_entry($otus[0],$otus[1]), 
   $alphaM->get_entry($otus[1],$otus[0]));
is($alphaM->get_entry($otus[1],$otus[2]), '1.1101');
is($alphaM->get_entry($otus[0],$otus[2]), '33.1197');

# codeml NSSites parsing
# for only 1 model

my $codeml_single = Bio::Tools::Phylo::PAML->new
    (-file => Bio::Root::IO->catfile(qw/t data singleNSsite.mlc/));
ok($codeml_single);
my $result_single = $codeml_single->next_result;
my ($nssite_single) = $result_single->get_NSSite_results;
is($nssite_single->num_site_classes,q/3/);
is($nssite_single->kappa, q/5.28487/);
is($nssite_single->likelihood,q/-30.819156/);

is($baseml->get_stat('loglikelihood'),-110.532715);
is($baseml->get_stat('constant_sites'),46);
is($baseml->get_stat('constant_sites_percentage'),'80.70');
is($baseml->model,'HKY85 dGamma (ncatG=5)');

# user trees
$baseml_p = Bio::Tools::Phylo::PAML->new
    (-file => Bio::Root::IO->catfile(qw(t data baseml.usertree)));
$baseml = $baseml_p->next_result;

@trees = $baseml->get_trees;
is(@trees, 1);
is($trees[0]->score, -129.328757);

# codeml NSSites parsing
# for branch site model/clade model

my $codeml_bs = Bio::Tools::Phylo::PAML->new
    (-file => Bio::Root::IO->catfile(qw/t data branchSite.mlc/));
ok($codeml_bs);
my $result_bs = $codeml_bs->next_result;
my ($nssite_bs) = $result_bs->get_NSSite_results;
is($nssite_bs->num_site_classes,q/4/);
my $class_bs = $nssite_bs->dnds_site_classes;
is($class_bs->{q/p/}->[1],q/0.65968/);
is($class_bs->{q/w/}->[1]->{q/background/},q/0.00000/);
is($class_bs->{q/w/}->[2]->{q/foreground/},q/999.00000/);

# Let's parse the RST file

my $paml = Bio::Tools::Phylo::PAML->new
    (-file => Bio::Root::IO->catfile(qw(t data codeml_lysozyme mlc)),
     -dir  => Bio::Root::IO->catfile(qw(t data codeml_lysozyme)));

$result = $paml->next_result;

my ($rst) = grep {$_->id eq 'node#8'} $result->get_rst_seqs;
ok($rst);
is($rst->seq, join('',qw(
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
is($site->{'anc_aa'}, 'A','ancestral AA');
is($site->{'anc_prob'}, '0.947','ancestral AA');
is($site->{'derived_aa'}, 'T','derived AA');

($node) = $first_tree->find_node(-id => '12');
@changes = $node->get_tag_values('changes');
($site) = grep { $_->{'site'} == 88 } @changes;
is($site->{'anc_aa'}, 'N','ancestral AA');
is($site->{'anc_prob'}, '0.993','ancestral AA');
is($site->{'derived_aa'}, 'D','derived AA');
is($site->{'derived_prob'}, '0.998','derived AA');

my $persite = $result->get_rst_persite;
# minus 1 because we have shifted so that array index matches site number
# there are 130 sites in this seq file
is(scalar @$persite -1, $result->patterns->{'-ls'}); 
# let's score site 1
$site = $persite->[2]; 
# so site 2, node 2 (extant)
is($site->[2]->{'codon'}, 'GTC');
is($site->[2]->{'aa'}, 'V');
# site 2, node 3
is($site->[3]->{'codon'}, 'ATC');
is($site->[3]->{'aa'}, 'I');

# ancestral node 9
is($site->[9]->{'codon'}, 'GTC');
is($site->[9]->{'aa'},    'V');
is($site->[9]->{'prob'},  '1.000');
is($site->[9]->{'Yang95_aa'},'V');
is($site->[9]->{'Yang95_aa_prob'},'1.000');

# ancestral node 10
is($site->[10]->{'codon'}, 'ATC');
is($site->[10]->{'aa'},    'I');
is($site->[10]->{'prob'},  '0.992');
is($site->[10]->{'Yang95_aa'},'I');
is($site->[10]->{'Yang95_aa_prob'},'0.992');


## PAML 3.15
$paml = Bio::Tools::Phylo::PAML->new(-file => Bio::Root::IO->catfile(qw(t data codeml315.mlc)) );
$result = $paml->next_result;

is($result->model, 'One dN/dS ratio');
like($result->version, qr'3\.15');
$MLmat = $result->get_MLmatrix;
$NGmat = $result->get_NGmatrix;

is($NGmat->[0]->[1]->{'omega'}, 0.2264);
is($NGmat->[0]->[1]->{'dN'}, 0.0186);
is($NGmat->[0]->[1]->{'dS'}, 0.0821);

is($MLmat->[0]->[1]->{'omega'}, 0.32693);
is($MLmat->[0]->[1]->{'dN'}, '0.0210');
is($MLmat->[0]->[1]->{'dS'}, 0.0644);
