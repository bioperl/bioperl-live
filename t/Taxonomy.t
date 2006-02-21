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
	eval {
		require Bio::DB::Taxonomy;
		require XML::Twig;
	};
	if ( $@ ) {
		$error = 1;
		warn "Unable to run tests because XML::Twig is not installed\n";
	}
	$NUMTESTS = 20;
	$error = 0;
	plan tests => $NUMTESTS;
}

END {
	foreach ( $Test::ntest..$NUMTESTS) {
		skip('Unable to complete Taxonomy tests',1);
	}
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
		for(my $i=0;$i < 12;$i++ ) {skip(1,'no species object could be created')}
	} else {
		ok($n->species, 'sapiens', 'species name tested');
		ok($n->genus, 'Homo', 'genus tested');
		ok($n->common_name, 'human', 'common name tested');
		ok($n->ncbi_taxid, 9606, 'taxonomy id tested');
		ok($n->division, 'Primates', 'division tested');
		ok($n->genetic_code, 1, 'genetic code tested');
		ok($n->mitochondrial_genetic_code, 2, 'mitochondrial genetic code');
		ok(defined $n->pub_date);
		ok(defined $n->create_date);
		ok(defined $n->update_date);
		ok($n->rank, 'species', 'rank');	
		ok($n->parent_id, 9605,'parent id');
	}
	sleep(3);
	my $yeastid = $db->get_taxonid('Saccharomyces cerevisiae');
	ok($yeastid, 4932);

	# create a non-species node
	my $human_parent = $db->get_Taxonomy_Node('9605');
	ok($human_parent);
	unless( $human_parent ) {
		for( my $i = 0; $i < 3; $i++ ) {
			skip(1,'no species object could be created');
		}
	} else {
		ok($human_parent->parent_id, '207598');
		ok(($human_parent->classification)[0], 'Homo');
		ok($human_parent->genetic_code, 1);
	}
}
