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

use Data::Dumper; #print Dumper $m;

# match within
$pos = Bio::Location::Simple->new 
    (-start => 25, -end => 25, -strand=> 1 );
$res = $m->map($pos);
#use Data::Dumper; print Dumper $res;
ok $res->start, 20;
ok $res->end, 20;
ok $res->strand, 1;
ok $res->seq_id, 'peptide';

ok $m->swap;
#use Data::Dumper; print Dumper $m;
$pos = Bio::Location::Simple->new 
    (-start => 5, -end => 5, -strand=> 1 );
$res = $m->map($pos);
ok $res->start, 10;

# cds -> propeptide
ok $m->in('cds'), 'cds';
ok $m->out('propeptide'), 'propeptide';

#print Dumper $m;
$res = $m->map($pos);
#print Dumper $m;
ok $res->start, 2;
#print Dumper $res;
ok $res = $m->_translate($pos);
#print Dumper $res;
ok $res = $m->_reverse_translate($pos);
#print Dumper $res;


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
#use Data::Dumper; print Dumper $genestruct;
#print ref ($genestruct), "\n";
$m->exons($genestruct);

$m->in('chr');
ok $m->gene_offset(7), 7;
ok $res = $m->map($pos);
ok $res->start, -3;
#use Data::Dumper; print Dumper $m;

ok $m->gene_offset(3);
#use Data::Dumper; print Dumper $m;
ok $res = $m->map($pos), undef;
#use Data::Dumper; print Dumper $m;

#$m->dump;

# testing parameter handling in the constructor
ok $m = new Bio::Coordinate::GeneMapper(-in => 'gene',
					-out => 'peptide',
					-gene_offset => 3,
					-exons => $genestruct,
					-utr => 7,
					-peptide_offset => 5
				       );

#$m->dump;
