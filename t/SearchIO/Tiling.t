#-*-perl-*-
#$Id$
use strict;
use warnings;
    use vars qw($EXHAUSTIVE $VERBOSE);
BEGIN {
    use lib '.';
    use lib '../..';

    use Bio::Root::Test;
    $EXHAUSTIVE = $ENV{BIOPERL_TILING_EXHAUSTIVE_TESTS};
    $VERBOSE    = $ENV{BIOPERL_TILING_VERBOSE_TESTS};
    test_begin(-tests => ($EXHAUSTIVE ? 6519 : 1141) );
}

use_ok('Bio::Search::Tiling::MapTiling');
use_ok('Bio::Search::Tiling::MapTileUtils');
use_ok('Bio::SearchIO');
use_ok('Bio::Search::Hit::BlastHit');
use_ok('File::Spec');

my ($blio, $result, $hit, $tiling, $hsp);

my @normal_formats = qw( blast  wublast
                         blastn wublastn
                         blastp wublastp
                         multiblast 
                         megablast
                         rpsblast
                         psiblast );
my @xltd_formats  = qw(  blastx wublastx
                         tblastn wutblastn
                         tblastx wutblastx );
                         
                 
# an exhaustive listing of search reports in 
# t/data 
        
my %test_files = (
    'blast' => [qw(
               ecolitst.bls
               frac_problems.blast
               frac_problems2.blast
               frac_problems3.blast
               bl2seq.out
               )],
    'multiblast' => [qw(
               multi_blast.bls
               )],
    'blastn' => [qw(
               a_thaliana.blastn
               bl2seq.blastn
               new_blastn.txt
               hsinsulin.blastcl3.blastn
               )],
    'wublastn' =>[qw(
               brassica_ATH.WUBLASTN
               echofilter.wublastn
               )],
    'blastp' => [qw(
               blastp2215.blast
               no_hsps.blastp
               catalase-webblast.BLASTP
               )],
    'wublastp' => [qw(
               dcr1_sp.WUBLASTP
               ecolitst.wublastp
               contig-by-hand.wublastp
               ecolitst.noseqs.wublastp
               )],
    'blastx' => [qw(
               bl2seq.blastx.out
               )],
    'wublastx' => [qw(
               dnaEbsub_ecoli.wublastx
               )],
    'wublast' => [qw(
               tricky.wublast
               )],
    'tblastn' => [qw(
               tblastn.out
               1ZZ19XR301R-Alignment.tblastn
               )],
    'wutblastn' => [qw(
               dnaEbsub_ecoli.wutblastn
               )],
    'tblastx' => [qw(
               bl2seq.tblastx.out
               HUMBETGLOA.tblastx
               )],
    'wutblastx' => [qw(
               dnaEbsub_ecoli.wutblastx
               )],
    'megablast' => [qw(
               503384.MEGABLAST.2
               )],
    'rpsblast' => [qw(
               ecoli_domains.rpsblast
               )],
    'psiblast' => [qw(
               psiblastreport.out
               )],
    'bug2942'  => [qw(
               bug2942.blastx
               )]
    );

# a subset of search reports for 
# run-o-the-mill regression tests

my %example_files = (
    'blast' => [qw(
               ecolitst.bls
               )],
    'blastn' => [qw(
               a_thaliana.blastn
               )],
    'wublastn' =>[qw(
               brassica_ATH.WUBLASTN
               )],
    'blastp' => [qw(
               no_hsps.blastp
               catalase-webblast.BLASTP
               )],
    'wublastp' => [qw(
               dcr1_sp.WUBLASTP
               )],
    'blastx' => [qw(
               bl2seq.blastx.out
               )],
    'wublastx' => [qw(
               dnaEbsub_ecoli.wublastx
               )],
    'wublast' => [qw(
               tricky.wublast
               )],
    'tblastn' => [qw(
               tblastn.out
               )],
    'wutblastn' => [qw(
               dnaEbsub_ecoli.wutblastn
               )],
    'tblastx' => [qw(
               HUMBETGLOA.tblastx
               )],
    'wutblastx' => [qw(
               dnaEbsub_ecoli.wutblastx
               )],
    'megablast' => [qw(
               503384.MEGABLAST.2
               )]
    );

ok( $blio = new Bio::SearchIO( 
	-file=>test_input_file('dcr1_sp.WUBLASTP'),
	-format=>'blast'), 'parse data file');

$result = $blio->next_result;
while ( $_ = $result->next_hit ) {
    last if $_->name =~ /ASPTN/;
}
ok($hit = $_, 'got test hit');
ok($tiling = Bio::Search::Tiling::MapTiling->new($hit), 'create tiling');


# TilingI compliance

isa_ok($tiling, 'Bio::Search::Tiling::TilingI');
foreach ( qw( next_tiling rewind_tilings identities conserved length ) ) {
    ok( $tiling->$_, "implements '$_'" );
}

# regression test on original calculations

my @orig_id_results = ( 387,388,388,381,382,389 );
my @orig_cn_results = ( 622,619,628,608,611,613 );
my @id_results = (
    $tiling->identities('query', 'exact'),
    $tiling->identities('query', 'est'),
    $tiling->identities('query', 'max'),
    $tiling->identities('subject', 'exact'),
    $tiling->identities('subject', 'est'),
    $tiling->identities('subject', 'max')
    );
my @cn_results = (
    $tiling->conserved('query', 'exact'),
    $tiling->conserved('query', 'est'),
    $tiling->conserved('query', 'max'),
    $tiling->conserved('subject', 'exact'),
    $tiling->conserved('subject', 'est'),
    $tiling->conserved('subject', 'max')
    );
map { $_ = int($_) } @id_results, @cn_results;

is_deeply(\@id_results, \@orig_id_results, 'identities regression test');
is_deeply(\@cn_results, \@orig_cn_results, 'conserved regression test');

# tiling iterator regression tests

my ($qn, $sn)=(0,0);
while ($tiling->next_tiling('query')) {$qn++};
while ($tiling->next_tiling('subject')) {$sn++};
is ($qn, 8, 'tiling iterator regression test(1)');
is ($sn, 128, 'tiling iterator regression test(2)');
$tiling->rewind('subject');
while ($tiling->next_tiling('subject')) {$sn++};
is ($sn, 256, 'tiling iterator regression test(3, rewind)');

diag("Old blast.t tiling tests") if $VERBOSE;

ok($blio = Bio::SearchIO->new(
    '-format' => 'blast',
    '-file'   => test_input_file('ecolitst.wublastp')
   ), "ecolitst.wublastp");
$result = $blio->next_result;
$result->next_hit;
$hit = $result->next_hit;
$tiling = Bio::Search::Tiling::MapTiling->new($hit);
# Test HSP contig data returned by SearchUtils::tile_hsps()
# Second hit has two hsps that overlap.

# compare with the contig made by hand for these two contigs
# in t/data/contig-by-hand.wublastp
# (in this made-up file, the hsps from ecolitst.wublastp
#  were aligned and contiged, and Length, Identities, Positives 
#  were counted, by a human (maj) )
	
my $hand_hit = Bio::SearchIO->new(
    -format=>'blast', 
    -file=>test_input_file('contig-by-hand.wublastp')
    )->next_result->next_hit;
my $hand_hsp = $hand_hit->next_hsp;
my @hand_qrng = $hand_hsp->range('query');
my @hand_srng = $hand_hsp->range('hit');
my @hand_matches = $hand_hit->matches;

is(($tiling->range('query'))[0], $hand_qrng[0]);
is(($tiling->range('query'))[1], $hand_qrng[1]);
is(sprintf("%d",$tiling->identities('query')), $hand_matches[0]);
is(sprintf("%d",$tiling->conserved('query')), $hand_matches[1]);
is(($tiling->range('hit'))[0], $hand_srng[0]);
is(($tiling->range('hit'))[1], $hand_srng[1]);
is(sprintf("%d",$tiling->identities('hit')), $hand_matches[0]);
is(sprintf("%d",$tiling->conserved('hit')), $hand_matches[1]);

ok( $blio = Bio::SearchIO->new(
	'-format' => 'blast',
	'-file'   => test_input_file('dnaEbsub_ecoli.wublastx')
    ), "dnaEbsub_ecoli.wublastx");

$hit = $blio->next_result->next_hit;
$tiling = Bio::Search::Tiling::MapTiling->new($hit);
is(sprintf("%.3f",$tiling->frac_identical(-type=>'query',-denom=>'aligned',-context=>'p2')), '0.364');
is(sprintf("%.3f",$tiling->frac_identical(-type=>'hit',-denom=>'aligned',-context=>'all')), '0.366');
is(sprintf("%.3f",$tiling->frac_conserved(-type=>'query',-denom=>'aligned',-context=>'p2')), '0.537');
is(sprintf("%.3f",$tiling->frac_conserved(-type=>'hit',-denom=>'aligned',-context=>'all')), '0.540');
is(sprintf("%.2f",$tiling->frac_aligned_query(-context=>'p2')), '0.62');
is(sprintf("%.2f",$tiling->frac_aligned_hit(-context=>'all')), '0.71');

ok( $blio = Bio::SearchIO->new(
	'-format' => 'blast',
	'-file'   => test_input_file('tricky.wublast')
    ), "tricky.wublast");

$hit = $blio->next_result->next_hit;
$tiling = Bio::Search::Tiling::MapTiling->new($hit);
cmp_ok sprintf("%.3f",$tiling->frac_identical(-denom => 'aligned')), '>', 0.2, 'tricky.wublast(1)';
cmp_ok sprintf("%.3f",$tiling->frac_conserved(-denom => 'aligned')), '<=', 1, 'tricky.wublast(2)';
is(sprintf("%.2f",$tiling->frac_aligned_query), '0.92', 'tricky.wublast(3)');
is(sprintf("%.2f",$tiling->frac_aligned_hit), '0.91','tricky.wublast(4)');

diag("New tiling tests") if $VERBOSE;

# select test file set based on the environment variable
# BIOPERL_TILING_EXHAUSTIVE_TESTS

my $files = ($EXHAUSTIVE ? \%test_files : \%example_files);

foreach my $alg (@normal_formats, @xltd_formats) {
    diag("*******$alg files*******") if ($files->{$alg} && $VERBOSE);
    foreach my $tf (@{$files->{$alg}}) {
	ok( $blio = Bio::SearchIO->new( -format=>'blast', 
					-file=>test_input_file($tf)
	    ), "$tf" );
	$result = $blio->next_result;
	my $hit_count = 0;
	# compare the per-aligned-base identity avg over hsps
	# with frac_identical (bzw, conserved)
	
      HIT:
	while ( $hit = $result->next_hit ) {
	    ++$hit_count;
	    # quiet the "No HSPs" warning with -verbose => -1
	    ok( $tiling = Bio::Search::Tiling::MapTiling->new(-hit=>$hit,-verbose=>-1), "tile $tf hit $hit_count #hsps ".scalar $tiling->hsps );
	    my @hsps = $tiling->hsps;
	    
	    unless (@hsps) {
		diag( "--no hsps for $tf hit $hit_count") if $VERBOSE;
		next HIT;
	    }
	    my ($dpct, $est, $fast,$exact, $max);
	    my $tol = 0.10; # % difference accepted as approx. equal
	    
	    ## loop through contexts:
	    for my $type (qw( query hit )) {
		for my $context ($tiling->contexts($type)) {
		    diag(" --- $type $context ---") if $VERBOSE;
		    if (scalar($tiling->contexts($type, $context)) == 1) {
			# equality
			($dpct, $est, $fast) = $tiling->cmp_frac($type,'identical','aligned', 'est', 'fast', $context);
			is( $est,$fast, substr($type,0,1)." id: est ($est) = fast ($fast)");
			($dpct, $est, $fast) = $tiling->cmp_frac($type,'conserved','aligned', 'est', 'fast', $context);
			is( $est,$fast, substr($type,0,1)." cn: est ($est) = fast ($fast)");
		    }
		    else {
			# comparisons
			($dpct, $est, $fast) = $tiling->cmp_frac($type,'identical','aligned', 'est', 'fast', $context);
#			cmp_ok( $dpct, "<", $tol, 
#				substr($type,0,1)." id: est ($est) ~ fast ($fast)");
			($dpct, $exact, $max) = $tiling->cmp_frac($type,'identical','aligned', 'exact', 'max', $context);
			cmp_ok( abs($exact-$est)/$exact, "<" , $tol, 
				substr($type,0,1)." id: exact ($exact) ~ est ($est)");
			cmp_ok( $exact, "<=" , $max, 
				substr($type,0,1)." id: exact ($exact) <= max ($max)");
			
			($dpct, $est, $fast) = $tiling->cmp_frac($type,'conserved','aligned', 'est', 'fast', $context);
#			cmp_ok( $dpct, "<", $tol, 
#				substr($type,0,1)." cn: est ($est) ~ fast ($fast)");
			($dpct, $exact, $max) = $tiling->cmp_frac($type,'conserved','aligned', 'exact', 'max', $context);
			cmp_ok(  abs($exact-$est)/$exact, "<" , $tol, 
				 substr($type,0,1)." cn: exact ($exact) ~ est ($est)");
			cmp_ok( $exact, "<=" , $max, 
				substr($type,0,1)." cn: exact ($exact) <= max ($max)");
		    }
		}
	    }
	}
    }
}

# bug 2942

my %expected_ranges = ( 'm0' => [7, 11037], #query
			'm1' => [1770, 10865], #query 
			'm2' => [2462, 14599], #query
			'all' => [231, 3563] #subject
    );
$blio = Bio::SearchIO->new( -file=>test_input_file( $test_files{'bug2942'}->[0] ),
			    -format => 'blast' );
$hit = $blio->next_result->next_hit;
$tiling = Bio::Search::Tiling::MapTiling->new($hit);
for ( 'm0', 'm1', 'm2' ) {
    is_deeply( [$tiling->range('query',$_)], $expected_ranges{$_}, "bug2942: query $_: range correct");
}
is_deeply( [$tiling->range('subject', 'all')], $expected_ranges{'all'}, "bug2942: subject all : range correct" );

# test get_tiled_alns

$blio = Bio::SearchIO->new( -file=>test_input_file( 'dcr1_sp.WUBLASTP' ) );
$result = $blio->next_result;
while ($hit = $result->next_hit) {
    last if $hit->name =~ /ASPTN/;
}

$tiling = Bio::Search::Tiling::MapTiling->new($hit);

ok my @alns = $tiling->get_tiled_alns, "get_tiled_alns";
is scalar @alns, 6, "got all alns";

for my $aln ( @alns ) {
    my (@aint, @qint, @sint);
    my $qs = $aln->get_seq_by_id('query');
    my $ss = $aln->get_seq_by_id('subject');
    ok my @qfeats = $qs->get_SeqFeatures;
    foreach (@qfeats) {
	push @aint, [$_->start, $_->end];
	push @qint, [($_->get_tag_values('query_start'))[0],
		     ($_->get_tag_values('query_end'))[0] ];
    }
    is( eval(join('+', map {$$_[1]-$$_[0]+1} @aint)),
	eval(join('+', map {$$_[1]-$$_[0]+1} @qint)), 
	"aln and qfeat lengths correspond" );
    is( $qs->length - $qs->num_gaps('-'), eval(join('+', map {$$_[1]-$$_[0]+1} @qint)), "q length correct");
    ok my @hfeats = $ss->get_SeqFeatures;
    @aint = ();
    ok ( @qfeats == @hfeats, "features on q and s correspond");
    foreach (@hfeats) {
	push @aint, [$_->start, $_->end];
	push @sint, [($_->get_tag_values('subject_start'))[0],
		     ($_->get_tag_values('subject_end'))[0] ];
    }
    is( eval(join('+', map {$$_[1]-$$_[0]+1} @aint)),
	eval(join('+', map {$$_[1]-$$_[0]+1} @sint)), 
	"aln and hfeat lengths correspond" );
    is( $ss->length - $ss->num_gaps('-'), eval(join('+', map {$$_[1]-$$_[0]+1} @sint)), "s length correct");

}
1;



package Bio::Search::Tiling::MapTiling;

sub cmp_frac {
    my ($tiling, $type, $method, $denom, @actions) = @_;
    my ($a, $b);
    my $context = ($actions[2] ? $actions[2] : 'all');
    $a = $tiling->frac(-type=>$type, 
		       -method=>$method, 
		       -denom=>$denom,
		       -action=>$actions[0],
		       -context=>$context);
    $b = $tiling->frac(-type=>$type, 
		       -method=>$method, 
		       -denom=>$denom,
		       -action=>$actions[1],
		       -context=>$context);
    return ( abs($a-$b)/$a, f(5,$a), f(5,$b) );
}

sub f { my ($d,$val) = @_; sprintf("%.${d}f",$val) }       

    

    
1;

