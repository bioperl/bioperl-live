# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## # $Id$

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

    eval { require 'Graph.pm' };
    if( $@ ) {
	    print STDERR "\nGraph.pm doesn't seem to be installed on this system -- the GO Parser needs it...\n\n";
	    plan tests => 1;
	    ok( 1 );
	    exit( 0 );
    }

    plan tests => 6;
}


use Bio::Ontology::OntologyStore;

my $store = Bio::Ontology::OntologyStore->get_instance;
ok($store);

my $ontology;
eval {
  $ontology = $store->get_ontology(-name => 'Sequence Ontology');
};

if($@){
  skip("couldn't get sequence ontology, network down? $@",5);  
} else {
  ok('got file okay');
  ok(scalar($ontology->get_root_terms()) == 1);


  my($txt) = $ontology->find_terms(-name => 'transcript');
  ok($txt->identifier eq 'SO:0000673');
  ok($txt->name eq 'transcript');
  ok($txt->definition eq 'An RNA synthesized on a DNA or RNA template by an RNA polymerase.');
}
