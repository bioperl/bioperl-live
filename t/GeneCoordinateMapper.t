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

    plan tests => 42;
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

my ($c, $value);


#
# Extrapolating pairs
#
#    No gaps returned, matches extrapolated
#     returns always a match or undef
#     -strict
#
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

ok $res->match->start, 5;
ok $res->match->end, 5;
ok $res->match->strand, -1;
ok $res->match->seq_id, 'peptide';


# match outside = undef
$pos = Bio::Location::Simple->new (-start => 5, -end => 5 );
$res = $pair->map($pos);

ok $res, undef;
#ok $res->match->strand, 1; # no guessing of strand
#ok $res->match->start, -16;
#ok $res->match->length, $pos->length;
#ok $res->match->seq_id, 'peptide';

#
# partial match = match
#
$pos2 = Bio::Location::Simple->new
    (-start => 20, -end => 22, -strand=> -1 );

ok $res = $pair->map($pos2);

ok $res->match->start, 1;
ok $res->match->end, 2;
ok $res->match->seq_id, 'peptide';
ok $res->match->strand, -1;


#
# partial match2 =  match & gap
#
$pos2 = Bio::Location::Simple->new (-start => 40, -end => 41, -strand=> 1 );
ok $res = $pair->map($pos2);
ok $res->match->start, 20;
ok $res->match->end, 20;

#
#enveloping
#
$pos2 = Bio::Location::Simple->new (-start => 19, -end => 41, -strand=> 1 );
ok $res = $pair->map($pos2);
ok $res->match->start, 1;
ok $res->match->end, 20;


#
#
# Gene Mapper
#
#

use Bio::Coordinate::GeneMapper;

ok my $m = new Bio::Coordinate::GeneMapper(-in => 'propeptide',
					   -out => 'peptide');

ok $m->peptide_offset(5), 5;


# match within
$pos = Bio::Location::Simple->new 
    (-start => 25, -end => 25, -strand=> 1 );
$res = $m->map($pos);
ok $res->start, 20;
ok $res->end, 20;
ok $res->strand, 1;
ok $res->seq_id, 'peptide';

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
ok $res = $m->_reverse_translate($pos);


$pos = Bio::Location::Simple->new 
    (-start => 26, -end => 26, -strand=> 1 );
$m->in('transcript');
ok $m->utr(7), 7;
$m->out('peptide');
$res = $m->map($pos);

#$m->verbose(2);

#print "++++++++++++++++++++\n";
#
# Collection representing exons
#
#         5   9     10  14    15  19
#         1   5     6   10    11  15
#         |---|     |---|     |---|
#-----|-----------------------|---|--
#     1   5   9     15  19    25  29
#         pair1     pair2     pair3

# gene
my $e1 = Bio::Location::Simple->new 
    (-seq_id => 'gene', -start => 5, -end => 9, -strand=>1 );
my $e2 = Bio::Location::Simple->new 
    (-seq_id => 'gene', -start => 15, -end => 19, -strand=>1 );
my $e3 = Bio::Location::Simple->new 
    (-seq_id => 'gene', -start => 25, -end => 29, -strand=>1 );
my $genestruct = [$e1, $e2, $e3];
$m->exons($genestruct);

$m->in('chr');
ok $m->gene_offset(7), 7;
ok $res = $m->map($pos);
ok $res->start, -3;

ok $m->gene_offset(3);
ok $res = $m->map($pos), undef;

#$m->to_string;

# testing parameter handling in the constructor
ok $m = new Bio::Coordinate::GeneMapper(-in => 'gene',
					-out => 'peptide',
					-gene_offset => 3,
					-exons => $genestruct,
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

#exon1 :-5829 -  -5761
#exon2 : -277 -   463

my @gene2_dump = split ( /\n/, qq {
ENSG00000100079	34590438	34590464			-1
ENSG00000100079	34582387	34582469	34582387	34582469	-1
ENSG00000100079	34581114	34581273	34581114	34581273	-1
ENSG00000100079	34580784	34580950	34580804	34580950	-1
});

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


ok my $g1 = new Bio::Coordinate::GeneMapper(-in=>'chr', -out=>'gene');

my ($gene1, $gene2);
my ($gene_offset, $cds_start, $cds_end, $strand, @exons);

#one line per exon
my ($first, $first_line);
for my $line ( @gene1_dump ) {
    unless ($first) {$first = 1; next;}

    my ($geneid, $exon_start, $exon_end, $exon_cstart,
	$exon_cend, $exon_strand) = split /\t/, $line;

#    print join (' ', $geneid, $exon_start, $exon_strand), "\n";

    # first line has gene_offset;
    unless ($first_line) {
	# gene coordinate start is defined by the start of the first exon;
	$gene_offset =
	$gene_offset = $exon_start - $exon_strand;
	$strand = $exon_strand;
	$first_line = 1;
    }

    # CDS location in chromosome coordinates
    $cds_start = $exon_cstart if !$cds_start and $exon_cstart;
    $cds_end = $exon_cend if !$cds_end and $exon_cend;

    my $exon = Bio::Location::Simple->new
	(-seq_id => 'gene', -start => $exon_start, -end => $exon_end, -strand=>$strand );
    push @exons, $exon;
}


$g1->gene_offset($cds_start - $strand);
my @gene_exons = map {$g1->map($_)} @exons;
$g1->exons(\@gene_exons);


# map CDS location into transcript coordinates
my $cds_in_chr = Bio::Location::Simple->new
    (-seq_id => 'gene', -start => $cds_start, -end => $cds_end, -strand=>$strand );
$g1->out('transcript');
my $cds_in_transctipt = $g1->map($cds_in_chr);

#$g1->utr($cds_in_transctipt->start - 1); #pass it aither a scalar or a Bio::Location!!

#print Dumper $g1, $cds_in_transctipt, \@exons, \@gene_exons;

#$g1->to_string;

# =strict(1);
# transcript($trloc = start of first exon, end of last exon)
# cds ($trloc = start of first exon, end of last exon);
