# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

# Note this uses Test::More; this should catch the few perl versions w/o
# this test suite

use strict;
use lib '..','.','./blib/lib';
use vars qw($NUMTESTS $DEBUG $error);

BEGIN { 
	$NUMTESTS = 90;
	$error = 0;
	$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
	# this seems to work for perl 5.6 and perl 5.8
	eval {require Test::More;};
	
	if ($@) {
		use lib 't';
	}
	
	use Test::More;
	
	eval {
		require XML::Simple;
		require LWP::UserAgent;
	};
	
	if (!$DEBUG) {
		plan skip_all => 'Set BIOPERLDEBUG=1 to run tests';
	} elsif ($@) {
		plan skip_all => 'Requires LWP::UserAgent and XML::Simple; skipping...';
	} else {
		plan tests => $NUMTESTS;
	}
	use_ok('Bio::DB::EUtilities');
}

require_ok('LWP::UserAgent');
require_ok('XML::Simple');

# NOTE : Bio::DB::EUtilities is just a specialized pipeline to get any 
# data available via NCBI's Entrez interface, with a few convenience methods
# to get UIDs and other additional information.  All data returned
# using EFetch is raw (not Bioperl objects) and is meant to be piped into
# other Bioperl modules at a later point for further processing

# protein acc
my @acc = qw(MUSIGHBA1 P18584 CH402638);
# protein GI
my @ids = qw(1621261 89318838 68536103 20807972 730439);
# test search term
my $term = 'dihydroorotase AND human';

my ($eutil, $response);

# Simple EFetch

SKIP: {
	$eutil = Bio::DB::EUtilities->new(
                                    -db         => 'protein',
                                    -id         => $ids[0],
									-rettype 	=> 'fasta'
                                      );
		  
	isa_ok($eutil, 'Bio::DB::GenericWebDBI');
	eval {$response = $eutil->get_response; };
	skip("EFetch HTTP error: $@", 2) if $@;
	isa_ok($response, 'HTTP::Response');
	my $content = $response->content;
	like($content, qr(PYRR \[Mycobacterium tuberculosis H37Rv\]),
		 'EFetch: Fasta format');
	
	# reuse the webagent
	$eutil->id($ids[1]);
	$eutil->rettype('gb');
	eval {$response = $eutil->get_response; };
	skip("EFetch HTTP error: $@", 2) if $@;
	isa_ok($response, 'HTTP::Response');
	$content = $response->content;
	like($content, qr(^LOCUS\s+EAS10332),'EFetch: GenBank format');
}

# EPost->EFetch with History (Cookie)

SKIP: {
	$eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'epost',
                                    -db         => 'protein',
                                    -id         => \@ids,
                                      );
		  
	isa_ok($eutil, 'Bio::DB::GenericWebDBI');
	eval {$response = $eutil->get_response; };
	skip("EPost HTTP error: $@", 12) if $@;
	isa_ok($response, 'HTTP::Response');
	my $cookie = $eutil->next_cookie;
	isa_ok($cookie, 'Bio::DB::EUtilities::Cookie');
	
	# set for epost, esearch, elink
	is($cookie->eutil, 'epost', '$epost->cookie->eutil()');
	is($cookie->database, 'protein', '$epost->cookie->database()');
	
	# these are not set using epost
	is($cookie->elink_dbfrom, undef, '$epost->cookie->elink_dbfrom()');
	is($cookie->esearch_total, undef, '$epost->cookie->esearch_total()');
	is($cookie->esearch_query, undef, '$epost->cookie->esearch_query()');
	is($cookie->elink_queryids, undef, '$epost->cookie->elink_queryids()');
	is($cookie->elink_linkname, undef, '$epost->cookie->elink_linkname()');
	
	# check the actual cookie
	my ($webenv, $key) = @{ $cookie->cookie };
	like($webenv, qr{^\S{50}}, '$epost->cookie->cookie() WebEnv');
	like($key, qr{^\d+}, '$epost->cookie->cookie() query key');
	
	# can we fetch the sequences using the cookie
	my $efetch = Bio::DB::EUtilities->new(
								-cookie		=> $cookie,
								-rettype  	=> 'fasta'
								  );
	# look for fasta headers
	my $total = grep(m{^>.*$}, split "\n", $efetch->get_response->content);
	is($total, 5, 'EPost to EFetch');
}

# ESearch, ESearch History

SKIP: {
	$eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'esearch',
                                    -db         => 'protein',
                                    -term       => $term,
									-retmax		=> 100
                                      );
		  
	isa_ok($eutil, 'Bio::DB::GenericWebDBI');
	eval {$response = $eutil->get_response; };
	skip("ESearch HTTP error:$@", 3) if $@;
	isa_ok($response, 'HTTP::Response');
	
	# can't really check for specific ID's but can check total ID's returned
	my @esearch_ids = $eutil->get_ids;
	is(scalar(@esearch_ids), 100, '$esearch->get_ids()');
	
	cmp_ok($eutil->esearch_count, '>', 117, '$esearch->esearch_count()');

	# usehistory (get a cookie)
	$eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'esearch',
                                    -db         => 'protein',
									-usehistory => 'y',
                                    -term       => $term,
                                      );
	
	eval {$response = $eutil->get_response; };
	skip("ESearch HTTP error:$@", 11) if $@;
	my $cookie = $eutil->next_cookie;
	isa_ok($cookie, 'Bio::DB::EUtilities::Cookie');
	is($cookie->eutil, 'esearch', '$esearch->cookie->eutil()');
	is($cookie->database, 'protein', '$esearch->cookie->database()');
	cmp_ok($cookie->esearch_total, '>', 117, '$esearch->cookie->esearch_total()');
	is($cookie->esearch_query, $term, '$esearch->cookie->esearch_query()');
	
	# these are not set using esearch
	is($cookie->elink_dbfrom, undef, '$esearch->cookie->elink_dbfrom()');
	is($cookie->elink_queryids, undef, '$esearch->cookie->elink_queryids()');
	is($cookie->elink_linkname, undef, '$esearch->cookie->elink_linkname()');
	
	# check the actual cookie
	my ($webenv, $key) = @{ $cookie->cookie };
	like($webenv, qr{^\S{50}}, '$esearch->cookie->cookie() WebEnv');
	like($key, qr{^\d+}, '$esearch->cookie->cookie() query key');
	
	# can we fetch the sequences using the cookie?
	my $efetch = Bio::DB::EUtilities->new(
								-cookie		=> $cookie,
								-rettype  	=> 'fasta',
								-retmax 	=> 5
								  );
	# look for the fasta headers
	my $total = grep(m{^>.*$}, split "\n", $efetch->get_response->content);
	is($total, 5, 'ESearch to EFetch'); 
}

# EInfo

SKIP: {
	$eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'einfo',
                                    -db  		=> 'protein',
                                      );
	isa_ok($eutil, 'Bio::DB::GenericWebDBI');
	eval {$response = $eutil->get_response; };
	skip("EInfo HTTP error:$@", 4) if $@;
	isa_ok($response, 'HTTP::Response');
	like($response->content, qr(<eInfoResult>), 'EInfo response');
	is($eutil->einfo_dbs->[0], 'protein', '$einfo->einfo_dbs()');
	like($eutil->einfo_db_lastupdate, qr(\d{4}\/\d{2}\/\d{2}\s\d{2}:\d{2}),
		 '$einfo->einfo_db_lastupdate()');
	cmp_ok($eutil->einfo_db_count, '>', 9200000, '$einfo->einfo_db_count()');
	is($eutil->einfo_db_desc, 'Protein sequence record', '$einfo->einfo_db_desc()');
	my @links = $eutil->einfo_dblink_info;
	my @fields = $eutil->einfo_dbfield_info;
	is(scalar(@links), 30, '$einfo->einfo_dblink_info()');
	is(scalar(@fields), 24, '$einfo->einfo_dbfield_info()');

	my %field = ('SingleToken' => 'Y',
				'Hierarchy' => 'N',
				'IsDate' => 'N',
				'TermCount' => '0',
				'Description' => 'Unique number assigned to each sequence',
				'Name' => 'UID',
				'IsNumerical' => 'Y');
	my %link = ('DbTo' => 'cdd',
				'Description' => 'Link to conserved domains within a protein',
				'Name' => 'protein_cdd',
				'Menu' => 'Conserved Domain Links');

	eq_hash($fields[1], \%field, '$einfo->einfo_dbfield_info()');
	eq_hash($links[1], \%link, '$einfo->einfo_dblink_info()');
	
	# einfo_dbfield_info, einfo_dblink_info 
	
	# all databases (list)
	$eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'einfo',
                                      );
	
	eval {$response = $eutil->get_response; };
	skip("EInfo HTTP error:$@", 1) if $@;
	
	my @db = sort qw(pubmed  protein  nucleotide  nuccore  nucgss  nucest  structure
	genome  books  cancerchromosomes  cdd  domains  gene  genomeprj  gensat
	geo  gds  homologene  journals  mesh  ncbisearch  nlmcatalog  omia  omim
	pmc  popset  probe  pcassay  pccompound  pcsubstance  snp  taxonomy
	unigene  unists);
	
	my @einfo_dbs = sort $eutil->einfo_dbs;
	is_deeply(\@einfo_dbs, \@db, 'All EInfo databases');
}

# ELink (normal; one db, one dbfrom) - ElinkData tests

SKIP: {
	my $eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'elink',
                                    -db  		=> 'taxonomy',
									-dbfrom		=> 'protein',
									-id			=> \@ids,
                                      );
		  
	isa_ok($eutil, 'Bio::DB::GenericWebDBI');
	my $response;
	eval {$response = $eutil->get_response; };
	skip("ELink HTTP error:$@", 4) if $@;
	isa_ok($response, 'HTTP::Response');
	like($response->content, qr(<eLinkResult>), 'ELink response');
	my @ids2 = qw(350054 306537 273068 83332 1394);
	is_deeply([sort $eutil->get_ids], [sort @ids2],'$elink->get_ids()');
	
	# Linkset tests
	is($eutil->get_linkset_count, 1, '$elink->get_linkset_count()');
	my $linkobj = $eutil->next_linkset;
	isa_ok($linkobj, 'Bio::DB::EUtilities::ElinkData');
	is($linkobj->elink_dbfrom, 'protein', '$linkdata->elink_dbfrom()');
	is_deeply([sort $linkobj->elink_queryids],
			  [sort @ids], '$linkdata->elink_queryids()');
	is($linkobj->elink_command, 'neighbor', '$linkdata->elink_command()');
	my $db = $linkobj->next_linkdb;
	is($db, 'taxonomy', '$linkdata->next_linkdb()');
	is_deeply([sort $linkobj->get_LinkIds_by_db($db)],
			  [sort @ids2], '$linkdata->get_LinkIds_by_db($db)');	
}

# To be added:

# ELink (normal, multiple db) 

# ELink (normal, multiple db, cookies)

# ELink (multi_id)

SKIP: {
	my $eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'elink',
                                    -db  		=> 'taxonomy',
									-dbfrom		=> 'protein',
									-multi_id 	=> 1,
									-id			=> \@ids,
                                      );
		  
	isa_ok($eutil, 'Bio::DB::GenericWebDBI');
	my $response;
	eval {$response = $eutil->get_response; };
	skip("ELink HTTP error:$@", 4) if $@;
	isa_ok($response, 'HTTP::Response');
	like($response->content, qr(<eLinkResult>), 'ELink response');
	my @ids2 = qw(350054 306537 273068 83332 1394);
	eval{$eutil->get_ids;};
	ok($@,'$elink->get_ids()');
	
	# Linkset tests
	is($eutil->get_linkset_count, 5, '$elink->get_linkset_count()');
	my $ct = 0;
	my @qids;
	my @retids;
	# ids may not be returned in same order as array, so need to grab and sort
	while (	my $linkobj = $eutil->next_linkset) {
		isa_ok($linkobj, 'Bio::DB::EUtilities::ElinkData');
		is($linkobj->elink_dbfrom, 'protein', '$linkdata->elink_dbfrom()');
		is($linkobj->elink_command, 'neighbor', '$linkdata->elink_command()');
		push @qids, $linkobj->elink_queryids;
		while (	my $db = $linkobj->next_linkdb) {
			is($db, 'taxonomy', '$linkdata->next_linkdb()');
			push @retids, $linkobj->get_LinkIds_by_db($db);
		}
	}
	is_deeply([sort @qids], [sort @ids], '$linkdata->elink_queryids()');
	is_deeply([sort @retids], [sort @ids2], '$linkdata->get_LinkIds_by_db($db)');
}

# ELink (multi_id, cookies)

# ELink (scores)

# Although the other EUtilities are available, no postprocessing is done on the
# returned XML yet

SKIP: {
	my $eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'esummary',
                                    -db         => 'protein',
                                    -id		    => \@ids,
                                      );
		  
	isa_ok($eutil, 'Bio::DB::GenericWebDBI');
	my $response;
	eval {$response = $eutil->get_response; };
	skip("ESummary HTTP error:$@", 11) if $@;
	isa_ok($response, 'HTTP::Response');
	like($response->content, qr(<eSummaryResult>), 'ESummary response');
	
	$eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'egquery',
                                    -term		=> $term,
                                      );
		  
	isa_ok($eutil, 'Bio::DB::GenericWebDBI');
	eval {$response = $eutil->get_response; };
	skip("EGQuery HTTP error:$@", 11) if $@;
	isa_ok($response, 'HTTP::Response');
	like($response->content, qr(<eGQueryResult>), 'EGQuery response');
}