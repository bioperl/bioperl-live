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
		  
	isa_ok($efetch, 'Bio::DB::EUtilities::efetch');
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
	like($content, qr(^LOCUS\s+AAY95024),'EFetch: Fasta format');
}

# EPost->EFetch with History (Cookie)

SKIP: {
	my $epost = Bio::DB::EUtilities->new(
                                    -eutil      => 'epost',
                                    -db         => 'protein',
                                    -id         => \@ids,
                                      );
		  
	isa_ok($epost, 'Bio::DB::EUtilities::epost');
	my $response;
	eval {$response = $epost->get_response; };
	skip("EPost HTTP error: $@", 13) if $@;
	isa_ok($response, 'HTTP::Response');
	my $cookie = $epost->next_cookie;
	isa_ok($cookie, 'Bio::DB::EUtilities::Cookie');
	
	# set for epost, esearch, elink
	is($cookie->eutil, 'epost', 'eutil');
	is($cookie->database, 'protein', 'database');
	
	# these are not set using epost
	is($cookie->elink_dbfrom, undef, 'dbfrom');
	is($cookie->esearch_total, undef, 'total');
	is($cookie->esearch_query, undef, 'query');
	is($cookie->elink_queryids, undef, 'queryids');
	is($cookie->elink_linkname, undef, 'linkname');
	
	# check the actual cookie
	my ($webenv, $key) = @{ $cookie->cookie };
	like($webenv, qr{^\S{50}}, 'cookie WebEnv');
	like($key, qr{^\d+}, 'cookie query key');
	
	# can we fetch the sequences using the cookie
	my $efetch = Bio::DB::EUtilities->new(
								-cookie		=> $cookie,
								-rettype  	=> 'fasta'
								  );
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
		  
	isa_ok($esearch, 'Bio::DB::EUtilities::esearch');
	my $response;
	eval {$response = $esearch->get_response; };
	skip("ESearch HTTP error:$@", 11) if $@;
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
	is($cookie->eutil, 'esearch', 'eutil');
	is($cookie->database, 'protein', 'database');
	cmp_ok($cookie->esearch_total, '>', 117, 'total');
	is($cookie->esearch_query, $term, 'query');
	
	## these are not set using esearch
	is($cookie->elink_dbfrom, undef, 'dbfrom');
	is($cookie->elink_queryids, undef, 'queryids');
	is($cookie->elink_linkname, undef, 'linkname');
	
	## check the actual cookie
	my ($webenv, $key) = @{ $cookie->cookie };
	like($webenv, qr{^\S{50}}, 'cookie WebEnv');
	like($key, qr{^\d+}, 'cookie query key');
	
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

# To be added:

# EInfo

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
		  
	isa_ok($eutil, 'Bio::DB::EUtilities::esummary');
	my $response;
	eval {$response = $eutil->get_response; };
	skip("ESummary HTTP error:$@", 11) if $@;
	isa_ok($response, 'HTTP::Response');
	like($response->content, qr(<eSummaryResult>), 'ESummary response');
	
	$eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'egquery',
                                    -term		=> $term,
                                      );
		  
	isa_ok($eutil, 'Bio::DB::EUtilities::egquery');
	eval {$response = $eutil->get_response; };
	skip("EGQuery HTTP error:$@", 11) if $@;
	isa_ok($response, 'HTTP::Response');
	like($response->content, qr(<eGQueryResult>), 'EGQuery response');
	
	$eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'einfo',
                                    -db  		=> 'protein',
                                      );
		  
	isa_ok($eutil, 'Bio::DB::EUtilities::einfo');
	eval {$response = $eutil->get_response; };
	skip("EInfo HTTP error:$@", 11) if $@;
	isa_ok($response, 'HTTP::Response');
	like($response->content, qr(<eInfoResult>), 'EInfo response');
}