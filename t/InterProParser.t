# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $GNF: projects/gi/symgene/src/perl/seqproc/t/InterProParser.t,v 1.7 2003/02/07 22:05:58 pdimitro Exp $

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
my $error;
use vars qw($NUMTESTS $DEBUG $HAVEGRAPHDIRECTED $errmsg);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

BEGIN {
	eval { require Test::More; };
	if( $@ ) {
	   use lib 't/lib';
	}
	use Test::More; 
	$NUMTESTS = 48;
	eval {
	   require XML::Parser::PerlSAX;
	   require XML::Parser;
	   require Graph::Directed;
	};
	if ( $@ ) {
	   plan skip_all => "XML::Parser, XML::Parser::PerlSAX, or Graph::Directed not installed. This means that InterPro Ontology Parsing module is not usable. Skipping tests";
	} else {
	   plan tests => $NUMTESTS;
	}
	use_ok('Bio::OntologyIO');
	use_ok('Bio::Root::IO');
}

my $io = Bio::Root::IO->new();

my $ipp = Bio::OntologyIO->new( -format => 'interpro',
										  -file =>
										  $io->catfile('t','data','interpro_short.xml'),
										  -ontology_engine => 'simple' );
ok ($ipp);

my $ip;
while(my $ont = $ipp->next_ontology()) {
    # there should be only one ontology
    is ($ip, undef);
    $ip = $ont;
}
#print $ip->to_string."\n";
my @rt = sort { $a->name cmp $b->name; } $ip->get_root_terms();

# there should be 4 root terms: InterPro Domain, InterPro Family,
# InterPro Repeat, and InterPro PTM (Post Translational Modification)
#
# I added 2 more terms, Active_site and Binding_site. -Juguang
is (scalar(@rt), 6);

# every InterPro term should have an ontology,
foreach ($ip->get_leaf_terms, @rt) {
	ok ($_->ontology);
	is ($_->ontology->name, "InterPro",
		 "term ".$_->name." in ontology InterPro");
}

# there are 6 fully instantiated InterPro terms in total, which should be returned as the leafs
is (scalar($ip->get_leaf_terms()), 8);
# roots and leafs together:
is (scalar($ip->get_all_terms()), 13);

# descendants and children (synonymous here because of depth 1)
# note that the sort should have placed Domain first and Family second
is (scalar($ip->get_descendant_terms($rt[2])), 4); # 4 InterPro Domains
is (scalar($ip->get_child_terms($rt[2])), 4);      # dto.
is (scalar($ip->get_descendant_terms($rt[3])), 3); # 3 Interpro Family
is (scalar($ip->get_child_terms($rt[3])), 3);      # dto.

# test for ancestors and parents (synonymous here because of depth 1)
foreach my $t ($ip->get_leaf_terms) {
	# every InterPro term has exactly one parent - namely either 
	# Domain, Family, Repeat, or PTM(Post Transl. Modification)
	if (!($t->identifier eq "Repeat" || $t->identifier eq "PTM" 
			|| $t->identifier eq 'Active_site' || $t->identifier eq 'Binding_site')) {
		is (scalar($ip->get_parent_terms($t)), 1);
		is (scalar($ip->get_ancestor_terms($t)), 1);
	}
}

# test for secondary accession map
is(scalar(keys %{$ipp->secondary_accessions_map}), 2);
