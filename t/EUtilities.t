# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

# Note this uses Test::More; this should catch the few perl versions w/o
# this test suite

use strict;
use lib '..','.','./lib','./blib/lib';
use vars qw($NUMTESTS $DEBUG $error);

BEGIN { 
    $NUMTESTS = 756;
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
my @ids = sort qw(1621261 89318838 68536103 20807972 730439);

# test search term
my $term = 'dihydroorotase AND human';

my ($eutil, $response);

my %dbs = (taxonomy => 1,
           nucleotide =>1,
           pubmed => 1);
my %links = (protein_taxonomy => 1,
             protein_nucleotide => 1,
             protein_nucleotide_wgs => 1,
             protein_pubmed => 1,
             protein_pubmed_refseq => 1
             );
my %scores = (   1621261 =>   2147483647,
                20807972 =>          423,
                68536103 =>          554,
                  730439 =>          411,
                89318838 =>          725,);

# Simple EFetch

SKIP: {
    $eutil = Bio::DB::EUtilities->new(
                                    -db         => 'protein',
                                    -id         => $ids[0],
                                    -rettype    => 'fasta'
                                      );
          
    isa_ok($eutil, 'Bio::DB::GenericWebDBI');
    eval {$response = $eutil->get_response; };
    skip("EFetch HTTP error: $@", 2) if $@;
    isa_ok($response, 'HTTP::Response');
    my $content = $response->content;
    like($content, qr(PYRR \[Mycobacterium tuberculosis H37Rv\]),
         'EFetch: Fasta format');
    
    # reuse the EUtilities webagent
    $eutil->id($ids[1]);
    $eutil->rettype('gb');
    eval {$response = $eutil->get_response; };
    skip("EFetch HTTP error: $@", 2) if $@;
    isa_ok($response, 'HTTP::Response');
    $content = $response->content;
    like($content, qr(^LOCUS\s+NP_623143),'EFetch: GenBank format');
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
                                -cookie     => $cookie,
                                -rettype    => 'fasta'
                                  );
    # look for fasta headers
    my $total = grep(m{^>.*$}, split "\n", $efetch->get_response->content);
    is($total, 5, 'EPost to EFetch');
}

# ESummary

my %docsum = (1621261=> [ ['Caption','String','CAB02640'],
['Title','String','PROBABLE PYRIMIDINE OPERON REGULATORY PROTEIN PYRR '.
 '[Mycobacterium tuberculosis H37Rv]'],
['Extra','String','gi|1621261|emb|CAB02640.1|[1621261]'],
['Gi','Integer','1621261'],
['CreateDate','String','2003/11/21'],
['UpdateDate','String','2005/04/17'],
['Flags','Integer',''],
['TaxId','Integer','83332'],
['Status','String','live'],
['ReplacedBy','String',''],
['Comment','String',''], ],
20807972 => [ ['Caption','String','NP_623143'],
['Title','String','pyrimidine regulatory protein PyrR '.
 '[Thermoanaerobacter tengcongensis MB4]'],
['Extra','String','gi|20807972|ref|NP_623143.1|[20807972]'],
['Gi','Integer','20807972'],
['CreateDate','String','2002/05/09'],
['UpdateDate','String','2005/12/03'],
['Flags','Integer','512'],
['TaxId','Integer','273068'],
['Status','String','live'],
['ReplacedBy','String',''],
['Comment','String',''], ],
68536103 => [ ['Caption','String','YP_250808'],
['Title','String','putative pyrimidine operon regulatory protein '.
 '[Corynebacterium jeikeium K411]'],
['Extra','String','gi|68536103|ref|YP_250808.1|[68536103]'],
['Gi','Integer','68536103'],
['CreateDate','String','2005/07/04'],
['UpdateDate','String','2006/03/30'],
['Flags','Integer','512'],
['TaxId','Integer','306537'],
['Status','String','live'],
['ReplacedBy','String',''],
['Comment','String',''], ],
730439 => [ ['Caption','String','P41007'],
['Title','String','PyrR bifunctional protein '.
 '[Includes: Pyrimidine operon regulatory protein; '.
 'Uracil phosphoribosyltransferase (UPRTase)]'],
['Extra','String','gi|730439|sp|P41007|PYRR_BACCL[730439]'],
['Gi','Integer','730439'],
['CreateDate','String','1995/02/01'],
['UpdateDate','String','2006/07/25'],
['Flags','Integer',''],
['TaxId','Integer','1394'],
['Status','String','live'],
['ReplacedBy','String',''],
['Comment','String',''] ],
89318838 => [ ['Caption','String','EAS10332'],
['Title','String','Phosphoribosyltransferase '.
 '[Mycobacterium flavescens PYR-GCK]'],
['Extra','String','gi|89318838|gb|EAS10332.1|[89318838]'],
['Gi','Integer','89318838'],
['CreateDate','String','2006/03/09'],
['UpdateDate','String','2006/03/09'],
['Flags','Integer',''],
['TaxId','Integer','350054'],
['Status','String','live'],
['ReplacedBy','String',''],
['Comment','String',''] ]);

SKIP: {
	$eutil = Bio::DB::EUtilities->new(
									 -eutil      => 'esummary',
									 -db         => 'protein',
									 -id            => \@ids,
									   );
	isa_ok($eutil, 'Bio::DB::GenericWebDBI');
	
	eval {$response = $eutil->get_response; };
	skip("ESummary HTTP error:$@", 20) if $@;
	isa_ok($response, 'HTTP::Response');
	
	my @docs = $eutil->get_all_docsums();
	is(scalar(@docs), 5, '$esum->get_all_docsums()');
	while (my $ds = $eutil->next_docsum) {
		isa_ok($ds, 'Bio::DB::EUtilities::DocSum');
		
		my $id = $ds->esummary_id();
		ok(exists($docsum{$id}), '$docsum->esummary_id()');
		
		my @items = @{ $docsum{$id} };
		
		# iterate using item names
		
		for my $name ($ds->get_all_names()) {
			my $item = shift @items;
			my %data = $ds->get_item_by_name($name);
			is($data{Name}, $item->[0], 'get_all_names(),DocSum Name');
			is($data{Type}, $item->[1], 'get_all_names(),DocSum Type');
			is($data{Content}, $item->[2], 'get_all_names(),DocSum Content');
			is($ds->get_Type_by_name($name), $item->[1],
			   'get_all_names(),get_Type_by_name()');
			is($ds->get_Content_by_name($name), $item->[2],
			   'get_all_names(),get_Content_by_name()');
		}

		@items = @{ $docsum{$id} };
		
		# iterating through each item (only first 3)
		my $ct = 0;
		while (my %data = $ds->next_docsum_item()) {
			my $item = shift @items;
			is($data{Name}, $item->[0], 'next_docsum_item(),DocSum Name');
		    is($data{Type}, $item->[1], 'next_docsum_item(),DocSum Type');
		    is($data{Content}, $item->[2], 'next_docsum_item(),DocSum Content');
			last if $ct++ == 2;
		}
		
		#test rewind
		$ds->rewind_docsum_items;
		
		@items = @{ $docsum{$id} };
		# just check the first one...
		while (my %data = $ds->next_docsum_item) {
			my $item = shift @items;
			is($data{Name}, $item->[0], 'rewind_docsum_items(),DocSum Name');
		    is($data{Type}, $item->[1], 'rewind_docsum_items(),DocSum Type');
		    is($data{Content}, $item->[2], 'rewind_docsum_items(),DocSum Content');
			last;
		}
	}
}

# ESearch, ESearch History

SKIP: {
    $eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'esearch',
                                    -db         => 'protein',
                                    -term       => $term,
                                    -retmax     => 100
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
                                -cookie     => $cookie,
                                -rettype    => 'fasta',
                                -retmax     => 5
                                  );
    # look for the fasta headers
    my $total = grep(m{^>.*$}, split "\n", $efetch->get_response->content);
    is($total, 5, 'ESearch to EFetch'); 
}

# EInfo

SKIP: {
    $eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'einfo',
                                    -db         => 'protein',
                                      );
    isa_ok($eutil, 'Bio::DB::GenericWebDBI');
    eval {$response = $eutil->get_response; };
    skip("EInfo HTTP error:$@", 10) if $@;
    isa_ok($response, 'HTTP::Response');
    like($response->content, qr(<eInfoResult>), 'EInfo response');
    is($eutil->einfo_dbs->[0], 'protein', '$einfo->einfo_dbs()');
    like($eutil->einfo_db_lastupdate, qr(\d{4}\/\d{2}\/\d{2}\s\d{2}:\d{2}),
         '$einfo->einfo_db_lastupdate()');
    cmp_ok($eutil->einfo_db_count, '>', 9200000, '$einfo->einfo_db_count()');
    is($eutil->einfo_db_desc, 'Protein sequence record', '$einfo->einfo_db_desc()');
    my @links = $eutil->einfo_dblink_info;
    my @fields = $eutil->einfo_dbfield_info;
    cmp_ok(scalar(@links), '>',30, '$einfo->einfo_dblink_info()');
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
    
    # all databases (list)
    $eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'einfo',
                                      );
    
    eval {$response = $eutil->get_response; };
    skip("EInfo HTTP error:$@", 1) if $@;
    
    my @db = sort qw(pubmed  protein  nucleotide  nuccore  nucgss  nucest  structure
    genome  books  cancerchromosomes  cdd  domains  gene  genomeprj  gensat
    geo  gds  homologene  journals  mesh  ncbisearch  nlmcatalog  omia  omim
    pmc  popset  probe  pcassay  pccompound  pcsubstance  snp  taxonomy toolkit
    unigene  unists);
    
    my @einfo_dbs = sort $eutil->einfo_dbs;
    cmp_ok(scalar(@einfo_dbs), '>=', scalar(@db), 'All EInfo databases');
}

# ELink - normal (single ID array) - single db - ElinkData tests

SKIP: {
    $eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'elink',
                                    -db         => 'taxonomy',
                                    -dbfrom     => 'protein',
                                    -id         => \@ids,
                                      );
          
    isa_ok($eutil, 'Bio::DB::GenericWebDBI');
    eval {$response = $eutil->get_response; };
    skip("ELink HTTP error:$@", 10) if $@;
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

# ELink - normal (single ID array), multiple dbs 

SKIP: {
    # can use 'all' for db, but takes a long time; use named dbs instead
    $eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'elink',
                                    -db         => 'taxonomy,nucleotide,pubmed',
                                    -dbfrom     => 'protein',
                                    -id         => \@ids,
                                      );
    
    isa_ok($eutil, 'Bio::DB::GenericWebDBI');
    eval {$response = $eutil->get_response; };
    skip("ELink HTTP error:$@", 14) if $@;
    isa_ok($response, 'HTTP::Response');
    like($response->content, qr(<eLinkResult>), 'ELink response');;
    
    # This is designed to fail; grabbing IDs w/o knowing which DB
    # they belong to in a multiple DB search is fatal
    my @ids2;
    eval {@ids2 = $eutil->get_ids;};
    ok($@,'$elink->get_ids()');
    
    # Must grab the linkset first...
    is($eutil->get_linkset_count, 1, '$elink->get_linkset_count()');
    my $linkobj = $eutil->next_linkset;
    isa_ok($linkobj, 'Bio::DB::EUtilities::ElinkData');
    
    # then iterate through each database, grabbing the IDs for each database
    my %ids = (
            'taxonomy' => [sort qw(350054 306537 273068 83332 1394)],
            'nucleotide' => [sort qw(89318678 68535062 38490250 20806542)],
            'pubmed' => [sort qw(15968079 12368430 11997336 9634230 8206848)],
           );
    
    while (my $db = $linkobj->next_linkdb) {
        ok(exists $ids{$db}, "ElinkData database: $db");
        @ids2 = sort $linkobj->get_LinkIds_by_db($db);
        is_deeply($ids{$db}, \@ids2, "ElinkData database IDs: $db")
    }
    # other ElinkData methods
    is($linkobj->elink_dbfrom, 'protein', '$linkdata->elink_dbfrom()');
    is_deeply([sort $linkobj->elink_queryids],
              [sort @ids], '$linkdata->elink_queryids()');
    is($linkobj->elink_command, 'neighbor', '$linkdata->elink_command()');
}

# ELink - normal (single ID array), multiple dbs, cookies)

SKIP: {
    # can use 'all' for db, but takes a long time; use named dbs instead
    # this retrieves cookies instead (no ElinkData objects are stored)
    $eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'elink',
                                    -db         => 'taxonomy,nucleotide,pubmed',
                                    -dbfrom     => 'protein',
                                    -id         => \@ids,
                                    -cmd        => 'neighbor_history'
                                      );
    
    isa_ok($eutil, 'Bio::DB::GenericWebDBI');
    eval {$response = $eutil->get_response; };
    skip("ELink HTTP error:$@", 26) if $@;
    isa_ok($response, 'HTTP::Response');
    like($response->content, qr(<eLinkResult>), 'ELink response');;
    
    # This is designed to fail; grabbing IDs w/o knowing which DB
    # they belong to in a multiple DB search is fatal
    my @ids2;
    eval {@ids2 = $eutil->get_ids;};
    ok($@,'$elink->get_ids()');
    
    # No ElinkData objs
    is($eutil->get_linkset_count, 0, '$elink->get_linkset_count()');
    
    # There are ELink cookies instead
    is($eutil->get_cookie_count, 5, '$elink->get_cookie_count()');
    my $ct = 0;
    while (my $cookie = $eutil->next_cookie) {
        isa_ok($cookie, 'Bio::DB::EUtilities::Cookie');
        is($cookie->eutil, 'elink', '$elink->cookie->eutil()');
        ok(exists $dbs{$cookie->database},  '$elink->cookie->database()');
        is($cookie->elink_dbfrom, 'protein', '$elink->cookie->elink_dbfrom()');
        @ids2 = sort $cookie->elink_queryids;
        is_deeply(\@ids2, \@ids, '$elink->cookie->elink_queryids()');
        ok(exists $links{$cookie->elink_linkname}, '$elink->cookie->elink_linkname()');
        
        # these are not set using elink
        is($cookie->esearch_query, undef, '$elink->cookie->esearch_query()');
        is($cookie->esearch_total, undef, '$elink->cookie->esearch_total()');
        
        # check the actual cookie data
        my ($webenv, $key) = @{ $cookie->cookie };
        like($webenv, qr{^\S{50}}, '$elink->cookie->cookie() WebEnv');
        like($key, qr{^\d+}, '$elink->cookie->cookie() query key');
        
        # can we retrieve the data via efetch?  Test one...
        # Note the cookie has all the information contained to
        # retrieve data; no additional parameters needed
        if($cookie->database eq 'taxonomy') {
            my $efetch = Bio::DB::EUtilities->new(-cookie => $cookie);
            my $content = $efetch->get_response->content;
            like($content, qr(<TaxaSet>), 'ELink to EFetch : taxonomy');
        }
        last if $ct == 2;
        $ct++;
    }
}

# ELink (multi_id), single db
# this is a flag set to get one-to-one correspondence for ELink data

SKIP: {
    $eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'elink',
                                    -db         => 'taxonomy',
                                    -dbfrom     => 'protein',
                                    -multi_id   => 1,
                                    -id         => \@ids,
                                      );
          
    isa_ok($eutil, 'Bio::DB::GenericWebDBI');
    eval {$response = $eutil->get_response; };
    
    skip("ELink HTTP error:$@", 26) if $@;
    isa_ok($response, 'HTTP::Response');
    like($response->content, qr(<eLinkResult>), 'ELink response');
    my @ids2 = qw(350054 306537 273068 83332 1394);
    
    # This is designed to fail; IDs present in individual ElinkData objects
    # for one-to-one correspondence with ID groups
    eval{$eutil->get_ids;};
    ok($@,'$elink->get_ids()');
    
    # Linkset tests
    is($eutil->get_linkset_count, 5, '$elink->get_linkset_count()');
    my @qids;
    my @retids;
    my $ct = 0;
    # ids may not be returned in same order as array, so need to grab and sort
    while ( my $linkobj = $eutil->next_linkset) {
        isa_ok($linkobj, 'Bio::DB::EUtilities::ElinkData');
        is($linkobj->elink_dbfrom, 'protein', '$linkdata->elink_dbfrom()');
        is($linkobj->elink_command, 'neighbor', '$linkdata->elink_command()');
        push @qids, $linkobj->elink_queryids;
        while ( my $db = $linkobj->next_linkdb) {
            is($db, 'taxonomy', '$linkdata->next_linkdb()');
            push @retids, $linkobj->get_LinkIds_by_db($db);
        }
        last if $ct == 4;
        $ct++
    }
    is_deeply([sort @qids], [sort @ids], '$linkdata->elink_queryids()');
    is_deeply([sort @retids], [sort @ids2], '$linkdata->get_LinkIds_by_db($db)');
}

# ELink (multi_id, cookies)

# these need to be cleaned up

SKIP: {
    $eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'elink',
                                    -db         => 'taxonomy',
                                    -dbfrom     => 'protein',
                                    -multi_id   => 1,
                                    -id         => \@ids,
                                    -cmd        => 'neighbor_history'
                                      );
          
    isa_ok($eutil, 'Bio::DB::GenericWebDBI');
    eval {$response = $eutil->get_response; };
    
    # check this number, likely wrong
    skip("ELink HTTP error:$@", 27) if $@;
    isa_ok($response, 'HTTP::Response');
    like($response->content, qr(<eLinkResult>), 'ELink response');
    my @ids2 = qw(350054 306537 273068 83332 1394);
    
    # This is designed to fail; IDs present in individual ElinkData objects
    # for one-to-one correspondence with ID groups
    eval{$eutil->get_ids;};
    ok($@,'$elink->get_ids()');
    
    # Linkset tests (there aren't any)
    is($eutil->get_linkset_count, 0, '$elink->get_linkset_count()');
    is($eutil->get_cookie_count, 5, '$elink->get_cookie_count()');
    
    my $efetch = Bio::DB::EUtilities->new();
    my $ct = 0;
    while (my $cookie = $eutil->next_cookie) {
        isa_ok($cookie, 'Bio::DB::EUtilities::Cookie');
        is($cookie->eutil, 'elink', '$elink->cookie->eutil()');
        is($cookie->database, 'taxonomy',  '$elink->cookie->database()');
        is($cookie->elink_dbfrom, 'protein', '$elink->cookie->elink_dbfrom()');
        @ids2 = $cookie->elink_queryids;
        
        # should be single IDs, one per ElinkData obj
        is(scalar(@ids2), 1, '$elink->cookie->elink_queryids()');
        is($cookie->elink_linkname, 'protein_taxonomy',
           '$elink->cookie->elink_linkname()');
        # these are not set using elink
        is($cookie->esearch_query, undef, '$elink->cookie->esearch_query()');
        is($cookie->esearch_total, undef, '$elink->cookie->esearch_total()');
        
        # check the actual cookie data
        my ($webenv, $key) = @{ $cookie->cookie };
        like($webenv, qr{^\S{50}}, '$elink->cookie->cookie() WebEnv');
        like($key, qr{^\d+}, '$elink->cookie->cookie() query key');
        
        # can we retrieve the data via efetch?  Test one...
        # Note the cookie has all the information contained to
        # retrieve data; no additional parameters needed
        
        if($cookie->database eq 'taxonomy') {
            $efetch->add_cookie($cookie);
            my $content = $efetch->get_response->content;
            like($content, qr(<TaxaSet>), 'ELink to EFetch : taxonomy');
        }
        last if $ct == 2;
        $ct++;
    }
}

# ELink (multi_id, multidbs)

SKIP: {
    $eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'elink',
                                    -db         => 'taxonomy,nucleotide,pubmed',
                                    -dbfrom     => 'protein',
                                    -multi_id   => 1,
                                    -id         => \@ids,
                                      );
          
    isa_ok($eutil, 'Bio::DB::GenericWebDBI');
    eval {$response = $eutil->get_response; };
    
    # check this number, likely wrong
    skip("ELink HTTP error:$@", 20) if $@;
    isa_ok($response, 'HTTP::Response');
    like($response->content, qr(<eLinkResult>), 'ELink response');
    my @ids2 = qw(350054 306537 273068 83332 1394);
    
    # This is designed to fail; IDs present in individual ElinkData objects
    # for one-to-one correspondence with ID groups
    eval{$eutil->get_ids;};
    ok($@,'$elink->get_ids()');
    
    # Linkset tests (there aren't any)
    is($eutil->get_linkset_count, 5, '$elink->get_linkset_count()');
    is($eutil->get_cookie_count, 0, '$elink->get_cookie_count()');
    
    # Linkset tests
    my $ct = 0;
    while ( my $linkobj = $eutil->next_linkset) {
        isa_ok($linkobj, 'Bio::DB::EUtilities::ElinkData');
        is($linkobj->elink_dbfrom, 'protein', '$linkdata->elink_dbfrom()');
        is($linkobj->elink_command, 'neighbor', '$linkdata->elink_command()');
        my @dbs = $linkobj->get_all_linkdbs;
        cmp_ok(scalar(@dbs), '>=' , 2, '$linkobj->get_all_linkdbs()');
        while ( my $db = $linkobj->next_linkdb) {
            is($dbs{$db}, 1, '$linkdata->next_linkdb()');
            my @ids2 = $linkobj->get_LinkIds_by_db($db);
            cmp_ok(scalar(@ids2), '>=', 1, '$linkdata->get_LinkIds_by_db($db)');
        }
        last if $ct == 4;
        $ct++;
    }
}

# ELink (multi_id, multidb, cookies)

SKIP: {
    $eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'elink',
                                    -db         => 'taxonomy,nucleotide,pubmed',
                                    -dbfrom     => 'protein',
                                    -multi_id   => 1,
                                    -id         => \@ids,
                                    -cmd        => 'neighbor_history'
                                      );
          
    isa_ok($eutil, 'Bio::DB::GenericWebDBI');
    eval {$response = $eutil->get_response; };
    
    # check this number, likely wrong
    skip("ELink HTTP error:$@", 20) if $@;
    isa_ok($response, 'HTTP::Response');
    like($response->content, qr(<eLinkResult>), 'ELink response');
    my @ids2 = qw(350054 306537 273068 83332 1394);
    
    # This is designed to fail; IDs present in individual ElinkData objects
    # for one-to-one correspondence with ID groups
    eval{$eutil->get_ids;};
    ok($@,'$elink->get_ids()');
    
    # Linkset tests (there aren't any)
    is($eutil->get_linkset_count, 0, '$elink->get_linkset_count()');
    cmp_ok($eutil->get_cookie_count, '>', 15, '$elink->get_cookie_count()');
    my $ct = 0;
    while (my $cookie = $eutil->next_cookie) {
        isa_ok($cookie, 'Bio::DB::EUtilities::Cookie');
        is($cookie->eutil, 'elink', '$elink->cookie->eutil()');
        is($dbs{$cookie->database}, 1,  '$elink->cookie->database()');
        is($cookie->elink_dbfrom, 'protein', '$elink->cookie->elink_dbfrom()');
        @ids2 = $cookie->elink_queryids;
        
        # should be single IDs, one per ElinkData obj
        is(scalar(@ids2), 1, '$elink->cookie->elink_queryids()');
        is($links{$cookie->elink_linkname}, 1,
           '$elink->cookie->elink_linkname()');
        # these are not set using elink
        is($cookie->esearch_query, undef, '$elink->cookie->esearch_query()');
        is($cookie->esearch_total, undef, '$elink->cookie->esearch_total()');
        
        # check the actual cookie data
        my ($webenv, $key) = @{ $cookie->cookie };
        like($webenv, qr{^\S{50}}, '$elink->cookie->cookie() WebEnv');
        like($key, qr{^\d+}, '$elink->cookie->cookie() query key');
        last if $ct == 14;
        $ct++;
    }
}

# ELink (scores)

SKIP: {
    # an elink back to the same db (db eq dbfrom) returns similarity scores
    $eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'elink',
                                    -db         => 'protein',
                                    -dbfrom     => 'protein',
                                    -id         => $ids[0],
                                      );
          
    isa_ok($eutil, 'Bio::DB::GenericWebDBI');
    eval {$response = $eutil->get_response; };
    
    # check this number, likely wrong
    skip("ELink HTTP error:$@", 20) if $@;
    isa_ok($response, 'HTTP::Response');
    like($response->content, qr(<eLinkResult>), 'ELink response');
    
    # only one linkset, so this actually works (not recommended)
    my @ids2 = $eutil->get_ids;
    cmp_ok(scalar(@ids2), '>' ,765 ,'$elink->get_ids()');
    
    # Linkset tests (there aren't any)
    is($eutil->get_linkset_count, 1, '$elink->get_linkset_count()');
    is($eutil->get_cookie_count, 0, '$elink->get_cookie_count()');
    
    while ( my $linkobj = $eutil->next_linkset) {
        isa_ok($linkobj, 'Bio::DB::EUtilities::ElinkData');
        is($linkobj->elink_dbfrom, 'protein', '$linkdata->elink_dbfrom()');
        is($linkobj->elink_command, 'neighbor', '$linkdata->elink_command()');
        
        # get db with scores
        while ( my $db = $linkobj->next_scoredb) {
            is($db,'protein', '$linkdata->next_scoredb()');
            my @ids2 = $linkobj->get_LinkIds_by_db($db);
            cmp_ok(scalar(@ids2), '>', 765, '$linkdata->get_LinkIds_by_db($db)');
            for my $id (@ids) {
                is($linkobj->get_score($id), $scores{$id}, '$linkdata->get_score()');
            }
        }
    }
}

# Although the other EUtilities are available, no postprocessing is done on the
# returned XML yet

SKIP: {
    $eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'egquery',
                                    -term       => $term,
                                      );
          
    isa_ok($eutil, 'Bio::DB::GenericWebDBI');
    eval {$response = $eutil->get_response; };
    skip("EGQuery HTTP error:$@", 2) if $@;
    isa_ok($response, 'HTTP::Response');
    like($response->content, qr(<eGQueryResult>), 'EGQuery response');
}

