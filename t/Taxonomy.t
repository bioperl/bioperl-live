# This is -*-Perl-*- code
# $Id$
use strict;
use vars qw($NUMTESTS $DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
my $error;

BEGIN { 
    eval { require Test; };
    if( $@ ) { 
	use lib 't';
    }
    use Test;
    $NUMTESTS = 8;
    $error = 0;
    plan tests => $NUMTESTS;
}
END {
    foreach ( $Test::ntest..$NUMTESTS) {
	skip('unable to run all of the tests because XML::Twig is not installed',1);
    }
}

eval { 
    require Bio::DB::Taxonomy;
    require XML::Twig; };
if( $@ ) {
    $error = 1;
}

my $actually_submit = $DEBUG > 0;
if( $error ==  1 ) {
    exit(0);
}

my $db = new Bio::DB::Taxonomy(-source => 'entrez');

ok($db);

if( $actually_submit == 0 ) { 
    print STDERR "skipping Taxonomy tests to avoid blocking\n" if( $DEBUG);
    foreach ( $Test::ntest..$NUMTESTS) {
	skip('Skip to avoid blocking',1);
    }
}  else { 
    my ($taxonid) = $db->get_taxonid('homo sapiens');
    ok($taxonid, 9606);

    my $n = $db->get_Taxonomy_Node($taxonid);
    ok($n);
    unless( $n ) {
	for ( 1..4 ) { skip(1,'no species object could be created') }
    } else {
	ok($n->species, 'sapiens');
	ok($n->common_name, 'human');
	ok($n->ncbi_taxid, 9606);
	ok($n->division, 'mammals');
    }
    sleep(3);
    my $yeastid = $db->get_taxonid('Saccharomyces cerevisiae');
    ok($yeastid, 4932);
}
