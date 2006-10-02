# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## # $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'
use strict;
use vars qw($NUMTESTS $DEBUG);

BEGIN {
	$NUMTESTS = 7;
	$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
	
	eval {require Test::More;};
	if ($@) {
		use lib 't';
	}
	use Test::More;
	
	eval {
		require Graph;
	};
	if ($@) {
		plan skip_all => 'Graph not installed. This means that the module is not usable. Skipping tests';
	}
	else {
		plan tests => $NUMTESTS;
	}
	
	use_ok('Bio::Ontology::OntologyStore');
}

ok my $store = Bio::Ontology::OntologyStore->get_instance;

my $ontology;
eval {$ontology = $store->get_ontology(-name => 'Sequence Ontology');};
skip "Couldn't get sequence ontology, network problems? Skipping these tests", 5 if $@;

ok('got file okay');
ok(scalar($ontology->get_root_terms()) == 1);
my($txt) = $ontology->find_terms(-name => 'transcript');
is $txt->identifier, 'SO:0000673';
is $txt->name, 'transcript';
is $txt->definition, 'An RNA synthesized on a DNA or RNA template by an RNA polymerase.';
