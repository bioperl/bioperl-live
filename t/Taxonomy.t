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
	$NUMTESTS = 63;
	$error = 0;
	plan tests => $NUMTESTS;
}

END {
    unlink("t/data/taxdump/nodes");
    unlink("t/data/taxdump/parents");
    unlink("t/data/taxdump/id2names");
    unlink("t/data/taxdump/names2id");
	foreach ( $Test::ntest..$NUMTESTS) {
		skip('Unable to complete Taxonomy tests',1);
	}
}

if( $error ==  1 ) {
    exit(0);
}

my $db_entrez = new Bio::DB::Taxonomy(-source => 'entrez');
ok $db_entrez;

my $dir = "t/data/taxdump";
my $db_flatfile = new Bio::DB::Taxonomy(-source => 'flatfile',
                               -directory => $dir,
                               -nodesfile => Bio::Root::IO->catfile('t','data','taxdump','nodes.dmp'),
                               -namesfile => Bio::Root::IO->catfile('t','data','taxdump','names.dmp'),
                               -force => 1);
ok $db_flatfile;

my $n;
foreach my $db ($db_entrez, $db_flatfile) {
    my $id;
    eval {
        $id = $db->get_taxonid('Homo sapiens');
    };
    if ($@) {
        for (1..32) {
            skip(1, 'Skip unable to connect to entrez database; no network or server busy?');
        }
        next;
    }
    ok $id, 9606;
    
    # easy test on human, try out the main node methods
    ok $n = $db->get_Taxonomy_Node(9606);
    ok $n->object_id, 9606;
    ok $n->ncbi_taxid, $n->object_id;
    ok $n->parent_id, 9605;
    ok $n->rank, 'species';
    
    ok $n->node_name, 'Homo sapiens';
    ok $n->scientific_name, $n->node_name;
    ok (($n->classification)[0], $n->node_name);
    ok ${$n->name('scientific')}[0], $n->node_name;
    ok $n->binomial, $n->node_name;
    
    my %common_names = map { $_ => 1 } $n->common_names;
    ok keys %common_names, 2;
    ok exists $common_names{human};
    ok exists $common_names{man};
    
    ok $n->division, 'Primates';
    ok $n->genetic_code, 1;
    ok $n->mitochondrial_genetic_code, 2;
    # these are entrez-only, data not available in dmp files
    if ($db eq $db_entrez) {
        ok defined $n->pub_date;
        ok defined $n->create_date;
        ok defined $n->update_date;
    }
    
    #*** deprecated?
    ok $n->genus, 'Homo';
    ok $n->species, 'Homo sapiens';
    ok ! $n->sub_species;
    
    # There needs to be more in-depth testing of Bio::Taxonomy::Node, but I may
    # be rewriting it soon, so holding off for now
    
    sleep(3) if $db eq $db_entrez;
    
    # do some trickier things...
    ok $n = $db->get_Taxonomy_Node('89593');
    ok $n->scientific_name, 'Craniata';
    
    sleep(3) if $db eq $db_entrez;
    
    ok $n = $db->get_Taxonomy_Node('1760');
    ok $n->scientific_name, 'Actinobacteria';
    
    sleep(3) if $db eq $db_entrez;
    
    # entrez isn't as good at searching as flatfile, so we have to special-case
    my @ids = $db->get_taxonids('Chloroflexi');
    $db eq $db_entrez ? (ok @ids, 1) : (ok @ids, 2);
    $id = $db->get_taxonids('Chloroflexi (class)');
    ok $id, 32061;
    
    @ids = $db->get_taxonids('Rhodotorula');
    ok @ids, 8;
    @ids = $db->get_taxonids('Rhodotorula <Microbotryomycetidae>');
    ok @ids, 1;
    ok $ids[0], 231509;
}