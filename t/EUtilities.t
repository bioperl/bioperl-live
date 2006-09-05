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
	$NUMTESTS = 45;
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

# NOTE : Bio::DB::EUtilities is just a specialized pipeline to get any and all
# data available via NCBI's Entrez interface, with a few convenience methods
# to get UIDs and other additional information.  All data returned
# using EFetch is raw (not Bioperl objects) and is meant to be piped into
# other Bioperl modules at a later point for further processing

# protein acc
my @acc = qw(MUSIGHBA1 P18584 CH402638);
# protein GI
my @ids = (70733119,68347418,26991676,82736201);
# test search term
my $term = 'dihydroorotase AND human';

# Simple EFetch

SKIP: {
	my $efetch = Bio::DB::EUtilities->new(
                                    -db         => 'protein',
                                    -id         => $ids[0],
									-rettype 	=> 'fasta'
                                      );
		  
	isa_ok($efetch, 'Bio::DB::GenericWebDBI');
	my $response;
	eval {$response = $efetch->get_response; };
	skip("EFetch HTTP error: $@", 2) if $@;
	isa_ok($response, 'HTTP::Response');
	my $content = $response->content;
	like($content, qr(dihydroorotase \[Pseudomonas fluorescens Pf-5\]),'EFetch: Fasta format');
	
	# reuse the webagent
	$efetch->id($ids[1]);
	$efetch->rettype('gb');
	eval {$response = $efetch->get_response; };
	skip("EFetch HTTP error: $@", 2) if $@;
	isa_ok($response, 'HTTP::Response');
	$content = $response->content;
	like($content, qr(^LOCUS\s+AAY95024),'EFetch: GenBank format');
}

# EPost->EFetch with History (Cookie)

SKIP: {
	my $epost = Bio::DB::EUtilities->new(
                                    -eutil      => 'epost',
                                    -db         => 'protein',
                                    -id         => \@ids,
                                      );
		  
	isa_ok($epost, 'Bio::DB::GenericWebDBI');
	my $response;
	eval {$response = $epost->get_response; };
	skip("EPost HTTP error: $@", 12) if $@;
	isa_ok($response, 'HTTP::Response');
	my $cookie = $epost->next_cookie;
	isa_ok($cookie, 'Bio::DB::EUtilities::Cookie');
	
	# set for epost, esearch, elink
	is($cookie->eutil, 'epost', 'eutil()');
	is($cookie->database, 'protein', 'database()');
	
	# these are not set using epost
	is($cookie->elink_dbfrom, undef, 'elink_dbfrom()');
	is($cookie->esearch_total, undef, 'esearch_total()');
	is($cookie->esearch_query, undef, 'esearch_query()');
	is($cookie->elink_queryids, undef, 'elink_queryids()');
	is($cookie->elink_linkname, undef, 'elink_linkname()');
	
	# check the actual cookie
	my ($webenv, $key) = @{ $cookie->cookie };
	like($webenv, qr{^\S{50}}, 'cookie() WebEnv');
	like($key, qr{^\d+}, 'cookie() query key');
	
	# can we fetch the sequences using the cookie
	my $efetch = Bio::DB::EUtilities->new(
								-cookie		=> $cookie,
								-rettype  	=> 'fasta'
								  );
	# look for fasta headers
	my $total = grep(m{^>.*$}, split "\n", $efetch->get_response->content);
	is($total, 4, 'EPost->EFetch');
}

# ESearch, ESearch History

SKIP: {
	my $esearch = Bio::DB::EUtilities->new(
                                    -eutil      => 'esearch',
                                    -db         => 'protein',
                                    -term       => $term,
									-retmax		=> 100
                                      );
		  
	isa_ok($esearch, 'Bio::DB::GenericWebDBI');
	my $response;
	eval {$response = $esearch->get_response; };
	skip("ESearch HTTP error:$@", 3) if $@;
	isa_ok($response, 'HTTP::Response');
	
	my @ids = $esearch->get_ids;
	is(scalar(@ids), 100, 'ESearch retmax');
	
	cmp_ok($esearch->esearch_count, '>', 117, 'ESearch count');

	# usehistory (get a cookie)
	$esearch = Bio::DB::EUtilities->new(
                                    -eutil      => 'esearch',
                                    -db         => 'protein',
									-usehistory => 'y',
                                    -term       => $term,
                                      );
	
	eval {$response = $esearch->get_response; };
	skip("ESearch HTTP error:$@", 11) if $@;
	my $cookie = $esearch->next_cookie;
	isa_ok($cookie, 'Bio::DB::EUtilities::Cookie');
	is($cookie->eutil, 'esearch', 'eutil()');
	is($cookie->database, 'protein', 'database()');
	cmp_ok($cookie->esearch_total, '>', 117, 'esearch_total()');
	is($cookie->esearch_query, $term, 'esearch_query()');
	
	## these are not set using esearch
	is($cookie->elink_dbfrom, undef, 'elink_dbfrom()');
	is($cookie->elink_queryids, undef, 'elink_queryids()');
	is($cookie->elink_linkname, undef, 'elink_linkname()');
	
	## check the actual cookie
	my ($webenv, $key) = @{ $cookie->cookie };
	like($webenv, qr{^\S{50}}, 'cookie() WebEnv');
	like($key, qr{^\d+}, 'cookie() query key');
	
	# can we fetch the sequences using the cookie?
	my $efetch = Bio::DB::EUtilities->new(
								-cookie		=> $cookie,
								-rettype  	=> 'fasta',
								-retmax 	=> 5
								  );
	# look for the fasta headers
	my $total = grep(m{^>.*$}, split "\n", $efetch->get_response->content);
	is($total, 5, 'ESearch->EFetch'); 
}

# EInfo

SKIP: {
	my $eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'einfo',
                                    -db  		=> 'protein',
                                      );
		  
	isa_ok($eutil, 'Bio::DB::GenericWebDBI');
	my $response;
	eval {$response = $eutil->get_response; };
	skip("EInfo HTTP error:$@", 4) if $@;
	isa_ok($response, 'HTTP::Response');
	like($response->content, qr(<eInfoResult>), 'EInfo response');
	is($eutil->einfo_dbs, 'protein', 'einfo_dbs()');
	like($eutil->einfo_db_lastupdate, qr(\d{4}\/\d{2}\/\d{2}\s\d{2}:\d{2}),
		 'einfo_db_lastupdate()');
	
	# need tests for einfo_db_desc, einfo_db_count,
	# einfo_dbfield_info, einfo_dblink_info 
	
	# all databases
#	$eutil = Bio::DB::EUtilities->new(
#                                    -eutil      => 'einfo',
#                                    -db  		=> 'protein',
#                                      );
#	
#	eval {$response = $eutil->get_response; };
#	skip("EInfo HTTP error:$@", 5) if $@;
#	
#	my @db = qw(pubmed  protein  nucleotide  nuccore  nucgss  nucest  structure
#	genome  books  cancerchromosomes  cdd  domains  gene  genomeprj  gensat
#	geo  gds  homologene  journals  mesh  ncbisearch  nlmcatalog  omia  omim
#	pmc  popset  probe  pcassay  pccompound  pcsubstance  snp  taxonomy
#	unigene  unists);
#	
#	my @einfo_dbs = $eutil->einfo_dbs;
#	
#	for my $db (@db) {
#		is(grep(m{$db eq $_}, @einfo_dbs), 1, "DB: $db");
#	}
}

# To be added:

# ELink (normal) 

# ELink (normal, cookies) 

# ELink (multi_id)

# ELink (multi_id, cookies)

# ELink (scores)

# Although the other EUtilities are available, no postprocessing is done on the
# XML yet

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