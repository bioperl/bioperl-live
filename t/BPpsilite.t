# -*-Perl-*- mode

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
    plan tests => 11;
}

use Bio::Tools::BPlite::Iteration;
use Bio::Tools::BPpsilite;
use Bio::Root::IO;

my $report = Bio::Tools::BPpsilite->new(-file=>Bio::Root::IO->catfile("t","data", "psiblastreport.out"));
ok $report;
ok $report->query =~ /DICDI/;# " query not found";
ok $report->database =~ /swissprot/i;# " database name not found";

my $total_iterations = $report->number_of_iterations;
ok $total_iterations, 2, " wrong total iteration number";

my $last_iteration = $report->round($total_iterations);
my $oldhitarray_ref = $last_iteration->oldhits;

# Process initial newly identified hit only
my ($sbjct, $id, $new_hsp, @is_old);
 HIT: while($sbjct = $last_iteration->nextSbjct) {
	$id = $sbjct->name;
	@is_old =  grep  /\Q${id}\E/, @$oldhitarray_ref;
	next HIT if (@is_old);
 	$new_hsp = $sbjct->nextHSP;
	ok $new_hsp->score, 1097, " HSP score not found";
	last HIT;
 }
$report->close();
close FH;

# Verify parsing of PHI-PSI Blast reports
open FH, Bio::Root::IO->catfile("t","data","phipsi.out");
my $report2 = Bio::Tools::BPpsilite->new(-fh=>\*FH);

ok $report2;
ok $report2->pattern, "P-E-E-Q", " wrong phi pattern";
ok $report2->query_pattern_location->[1], 120, " wrong phi pattern location";

$total_iterations = $report2->number_of_iterations;
ok $total_iterations, 2, " wrong total iteration number in phiblast report";

my $last_iteration2 = $report2->round($total_iterations);
my $sbjct2 = $last_iteration2->nextSbjct;
ok $last_iteration2->newhits->[1] =~ /ARATH/;# " Hit not found in phiblast report";
my $hsp2 = $sbjct2->nextHSP;
ok $hsp2->hit->end, 343, " HSP start not found in phiblast report";
$report2->close();

close FH;





