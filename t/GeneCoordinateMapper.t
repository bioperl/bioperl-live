#-*-Perl-*-
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

    plan tests => 81;
}

use Bio::Location::Simple;
use Bio::Coordinate::Pair;
use Bio::Coordinate::Result;
use Bio::Coordinate::Result::Match;
use Bio::Coordinate::Result::Gap;
use Bio::Coordinate::ExtrapolatingPair;
use Bio::Coordinate::Collection;

use Data::Dumper;

use vars qw($DEBUG);
ok(1);



#
# Extrapolating pairs
#
#    No gaps returned, matches extrapolated
#     returns always a match or undef
#     -strict
#


# the  reverse strand pair
my $inr = new Bio::Location::Simple(-start=>2, -end=>5, -strand=>1);
my $outr = new Bio::Location::Simple(-start=>10, -end=>13, -strand=>-1);
ok my $pairr = Bio::Coordinate::ExtrapolatingPair->
    new(-in => $inr,
	-out => $outr
       );

my $posr = Bio::Location::Simple->new 
    (-start => 3, -end => 4, -strand=> 1 );
my $resr = $pairr->map($posr);
ok $resr->start, 11;
ok $resr->end, 12;
ok $resr->strand, -1;



# propepide
my $match1 = Bio::Location::Simple->new 
    (-seq_id => 'propeptide', -start => 21, -end => 40, -strand=>1 );
# peptide
my $match2 = Bio::Location::Simple->new
    (-seq_id => 'peptide', -start => 1, -end => 20, -strand=>1 );

ok my $pair = Bio::Coordinate::ExtrapolatingPair->
    new(-in => $match1,
	-out => $match2,
	-strict => 1
       );

ok $pair->test;
ok $pair->strand(), 1; #  = in->strand * out->strand
ok $pair->in->seq_id(), 'propeptide';
ok $pair->strict(), 1;

my ($count, $pos, $pos2, $res, $match, $res2);

# match within
$pos = Bio::Location::Simple->new 
    (-start => 25, -end => 25, -strand=> -1 );
$res = $pair->map($pos);

ok $res->isa('Bio::Location::Simple');
ok $res->start, 5;
ok $res->end, 5;
ok $res->strand, -1;
ok $res->seq_id, 'peptide';


# match outside = undef
$pos = Bio::Location::Simple->new (-start => 5, -end => 5 );
$res = $pair->map($pos);

ok $res, undef;

#
# partial match = match
#
$pos2 = Bio::Location::Simple->new
    (-start => 20, -end => 22, -strand=> -1 );

ok $res = $pair->map($pos2);

ok $res->start, 0;
ok $res->end, 2;
ok $res->seq_id, 'peptide';
ok $res->strand, -1;


#
# partial match2 =  match & gap
#
$pos2 = Bio::Location::Simple->new (-start => 40, -end => 41, -strand=> 1 );
ok $res = $pair->map($pos2);
ok $res->start, 20;
ok $res->end, 20;

#
#enveloping
#
$pos2 = Bio::Location::Simple->new (-start => 19, -end => 41, -strand=> 1 );
ok $res = $pair->map($pos2);
ok $res->start, 1;
ok $res->end, 20;

#
# testing the changing the strand
#

# chr
$match1 = Bio::Location::Simple->new 
    (-seq_id => 'chr', -start => 21, -end => 40, -strand=>1 );
# gene
$match2 = Bio::Location::Simple->new
    (-seq_id => 'gene', -start => 1, -end => 20, -strand=>-1 );

 $pair = Bio::Coordinate::ExtrapolatingPair->
#my $pair = Bio::Coordinate::Pair->
    new(-in => $match1,
	-out => $match2,
	-strict => 0
       );

$pos = Bio::Location::Simple->new 
    (-start => 38, -end => 40, -strand=> 1 );
$res = $pair->map($pos);
#print Dumper $res;
ok $res->start, 1;
ok $res->end, 3;
ok $res->strand, -1;

$pos = Bio::Location::Simple->new 
    (-start => 1, -end => 3, -strand=> 1 );
$res = $pair->map($pos);
#print Dumper $res;
ok $res->start, 38;
ok $res->end, 40;
ok $res->strand, -1;


#
#
# Gene Mapper
#
#

use Bio::Coordinate::GeneMapper;

ok my $m = new Bio::Coordinate::GeneMapper(-in => 'propeptide',
					   -out => 'peptide');
#$m->verbose(2);

ok $m->peptide_offset(5), 5;
#print Dumper $m;


# match within
$pos = Bio::Location::Simple->new 
    (-start => 25, -end => 25, -strand=> 1 );
$res = $m->map($pos);
#print Dumper $res;

ok $res->start, 20;
ok $res->end, 20;
ok $res->strand, 1;
ok $res->seq_id, 'peptide';


#
# nozero
#

# match within
$pos = Bio::Location::Simple->new 
    (-start => 4, -end => 5, -strand=> 1 );
$res = $m->map($pos);
ok $res->start, -1;
ok $res->end, 0;

ok $m->nozero('in&out'), 'in&out';
$res = $m->map($pos);
ok $res->start, -2;
ok $res->end, -1;
ok $m->nozero(0), 0;



ok $m->swap;
$pos = Bio::Location::Simple->new 
    (-start => 5, -end => 5, -strand=> 1 );
$res = $m->map($pos);
ok $res->start, 10;

# cds -> propeptide
ok $m->in('cds'), 'cds';
ok $m->out('propeptide'), 'propeptide';

$res = $m->map($pos);
ok $res->start, 2;
ok $res = $m->_translate($pos);
ok $res->start, 2;
ok $res = $m->_reverse_translate($pos);
ok $res->start, 13;
ok $res->end, 15;

$pos = Bio::Location::Simple->new 
    (-start => 26, -end => 26, -strand=> 1 );
$m->out('peptide');
$res = $m->map($pos);
ok $res->start, 4;


#
# frame
#

$pos = Bio::Location::Simple->new 
    (-start => 1, -end => 3, -strand=> 1 );
$res = $m->_frame($pos);
ok $res->start, 1;
ok $res->end, 3;


#$m->verbose(2);
#         5   9     10  14    15  19
#print "++++++++++++++++++++\n";
#
# Collection representing exons
#
#  cds    1   5     6   10    11  15
#  exon   1   5     1   5     1   5
#  gene   0   4    10   14   20   24
#         |---|     |---|     |---|
#-----|-----------------------|---|--
# chr 1   5   9    15   19   25   29
#         pair1     pair2     pair3

# gene
my $e1 = Bio::Location::Simple->new 
    (-seq_id => 'gene', -start => 5, -end => 9, -strand=>1 );
my $e2 = Bio::Location::Simple->new 
    (-seq_id => 'gene', -start => 15, -end => 19, -strand=>1 );
my $e3 = Bio::Location::Simple->new 
    (-seq_id => 'gene', -start => 25, -end => 29, -strand=>1 );
my @cexons = ($e1, $e2, $e3);

#$m->exons($genestruct);
$m= new Bio::Coordinate::GeneMapper;

$m->in('chr');
$m->out('gene');
my $off = $m->cds(5);
ok $off->start, 5; # start of the coding region
ok $m->exons(@cexons), 3;

$m->out('exon');
$pos = Bio::Location::Simple->new
    (-start => 6, -end => 7, -strand=> 1 );
$res = $m->map($pos);

ok $res->start, 2;
ok $res->end, 3;

$m->out('negative_intron');
#$m->out('exon');
$pos = Bio::Location::Simple->new 
    (-start => 12, -end => 14, -strand=> 1 );
$res = $m->map($pos);
#skip 'negative_intron', $res->start, -3;
#skip 'negative_intron', $res->end, -1;
#skip $res->seq_id, 'intron1';
#exit;

# cds
$m->out('cds');
$pos = Bio::Location::Simple->new
    (-start => 5, -end => 9, -strand=> 1 );
$res = $m->map($pos);
ok $res->start, 1;
ok $res->end, 5;

$pos = Bio::Location::Simple->new
    (-start => 15, -end => 19, -strand=> 1 );
$res = $m->map($pos);
ok $res->start, 6;
ok $res->end, 10;

$pos = Bio::Location::Simple->new
    (-start => 5, -end => 19, -strand=> 1 );
$res = $m->map($pos);
ok $res->start, 1;
ok $res->end, 10;
#$m->to_string;


ok $m->cds(3); # recalculating exons
#$m->to_string;

#
# Problem with negative numbers
#
# Collection representing exons
#
#  cds  -11  -7    -6  -2    -1   3  :27
#  cds   -6  -2    -1 1 3     4   8  :17
#  exon   1   5     1   5     1   5
#  gene -21  -17  -11  -7    -1 1 3  :27
#  gene -11  -7    -1 1 3     9   13 :17
#         |---|     |---|     |---|
#-----|-----------------------|---|--
# chr 1   5   9    15   19   25   29
#         pair1     pair2     pair3

$m= new Bio::Coordinate::GeneMapper;

$m->in('chr');
$m->out('gene');
$off = $m->cds(17);
ok $off->start, 17; # start of the coding region
ok $m->exons(@cexons), 3;
#$m->to_string;
#print Dumper $m;


# testing parameter handling in the constructor
ok $m = new Bio::Coordinate::GeneMapper(-in => 'gene',
					-out => 'peptide',
					-cds => 3,
					-exons => @cexons,
					-utr => 7,
					-peptide_offset => 5
				       );

#$m->to_string;

#
# Real life data
# Mapping SNPs into  human serum protein MSE55 and
# human galecting LGALS2 from Ensembl:
#

#Ensembl Gene ID	Exon Start (Chr bp)	Exon End (Chr bp)	Exon Coding Start (Chr bp)
#	Exon Coding End (Chr bp)	Strand

my @gene1_dump = split ( /\n/, qq {
ENSG00000128283	34571058	34571126			1
ENSG00000128283	34576610	34577350	34576888	34577350	1
ENSG00000128283	34578646	34579858	34578646	34579355	1
});


my @gene2_dump = split ( /\n/, qq {
ENSG00000100079	34590438	34590464			-1
ENSG00000100079	34582387	34582469	34582387	34582469	-1
ENSG00000100079	34581114	34581273	34581114	34581273	-1
ENSG00000100079	34580784	34580950	34580804	34580950	-1
}); # exon start should be less than end or is this intentional?

#Chromosome Name	Location (bp)	Strand	Reference ID
my @snp_dump = split ( /\n/, qq {
22	34572694	1	2235335
22	34572799	1	2235336
22	34572843	1	2235337
22	34574896	1	2076087
22	34575256	1	2076088
22	34578830	1	2281098
22	34579111	1	2281099
22	34580411	1	2235338
22	34580591	1	2281097
22	34580845	1	2235339
22	34581963	1	2281100
22	34583722	1	140057
22	34585003	1	140058
22	34587726	1	968725
22	34588207	1	2284055
22	34591507	1	1969639
22	34591949	1	140059
});
shift @snp_dump;

my ($cdsr, @exons) = read_gene_data(@gene1_dump);

ok my $g1 = new Bio::Coordinate::GeneMapper(-in=>'chr', -out=>'gene');
$g1->cds($cdsr);

#$pos = Bio::Location::Simple->new
#    (-start => 34576888, -end => 34579355, -strand=> 1 );
$res = $g1->map($cdsr);
ok $res->start, 1;
ok $res->end, 2468;

$g1->exons(@exons);
$g1->in('gene');
$g1->out('cds');
$res = $g1->map($res);
ok $res->start, 1;
ok $res->end, 1173;

#map_snps($g1, @snp_dump);


#gene 2 in reverse strand
($cdsr, @exons) = read_gene_data(@gene2_dump);
ok my $g2 = new Bio::Coordinate::GeneMapper(-in=>'chr', -out=>'gene');
$g2->cds($cdsr);

$pos = Bio::Location::Simple->new
    (-start => $cdsr->end-2, -end => $cdsr->end, -strand=> 1 );
$res = $g2->map($pos);
ok $res->start, 1;
ok $res->end, 3;
ok $res->strand, -1;
#print Dumper \@exons;

$g2->exons(@exons);

#map_snps($g2, @snp_dump);



#todo:
#  gene in opposite strand,
#  frame,
#  negative exon coordinates,
#  negative positions, no zero as in human genetics
#  strict mapping mode
#  test swapping
#  extrapolating pair code into Bio::Coordinate::Pair
#  split gene->inex mappings into start end... just athought




sub read_gene_data {
    my ($self,@gene_dump) = @_;
    my ($cds_start, $cds_end, $strand, @exons);

    #one line per exon
    my ($first, $first_line);
    for my $line ( @gene_dump ) {

	my ($geneid, $exon_start, $exon_end, $exon_cstart,
	    $exon_cend, $exon_strand) = split /\t/, $line;

	$strand = $exon_strand if $exon_strand;
	#print join (' ', $geneid, $exon_start, $exon_strand), "\n";

	# CDS location in chromosome coordinates
	$cds_start = $exon_cstart if !$cds_start and $exon_cstart;
	$cds_end = $exon_cend if $exon_cend;


	if ($exon_start > $exon_end) {
	    ($exon_start, $exon_end) = ($exon_end, $exon_start);
	}

	my $exon = Bio::Location::Simple->new
	    (-seq_id => 'gene', -start => $exon_start,
	     -end => $exon_end, -strand=>$strand, -verbose=>2);
	push @exons, $exon;
    }

    if ($cds_start > $cds_end) {
	($cds_start, $cds_end) = ($cds_end, $cds_start);
    }

    my $cdsr = Bio::Location::Simple->new (-start => $cds_start,
					   -end => $cds_end,
					   -strand=> $strand);

    return ($cdsr, @exons);
}


sub map_snps {
    my ($mapper, @snps) =@_;
    $mapper->in('chr');
    $mapper->out('cds');
    foreach my $line (@snps) {
	$mapper->out('cds');

	my ($chr, $start, $strand, $id) = split /\t/, $line;
	my $loc = Bio::Location::Simple->new
	    ( -start => $start,
	     -end => $start, -strand=>$strand );

	my $res = $mapper->map($loc);
	my $cds_start = 0;
	$cds_start = $res->start if defined $res;#defined $res->start;
	print $id, "\t", $cds_start, "\n";

	# coding
	if ($cds_start) {
	    $mapper->out('propeptide');
	    my $frame_obj = $mapper->_frame($res);
	    my $res = $mapper->map($loc);
	    my $cds_start = 0;
	    $cds_start = $res->start if defined $res;#defined $res->start;
	    print  "\t\t", $cds_start, " (", $frame_obj->start, ")\n";

	}

    }


}
