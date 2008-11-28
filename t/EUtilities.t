# -*-Perl-*- Test Harness script for Bioperl
# $Id$
#
# Note this uses Test::More; this should catch the few perl versions w/o
# this test suite

use strict;
our $NUMTESTS;
our $DEBUG;
our %EUTILS;

BEGIN {
    $NUMTESTS = 3;

    # I have set up eutils tests to run in sections for easier test maintenance
    # and keeping track of problematic tests. The below hash is the list of
    # tests, with test number and coderef.

    %EUTILS = (
        #'efetch'        => {'tests' => 5,
        #                    'sub'   => \&efetch},
        #'epost'         => {'tests' => 13,
        #                    'sub'   => \&epost},
        #'esummary'      => {'tests' => 254,
        #                    'sub'   => \&esummary},
        #'esearch'       => {'tests' => 15,
        #                    'sub'   => \&esearch},
        #'einfo'         => {'tests' => 10,
        #                    'sub'   => \&einfo},
        #'elink1'        => {'tests' => 9,
        #                    'sub'   => \&elink1},
        #'egquery'       => {'tests' => 3,
        #                    'sub'   => \&egquery},
        
        # The following tests either fail sporadically due to unknown client- or
        # server-side issues, contain volatile data, or are still being worked
        # on; uncomment to test
        
        #'elink2'        => {'tests' => 18,
        #                    'sub'   => \&elink2},
        #'elink3'        => {'tests' => 28,
        #                    'sub'   => \&elink3},
        #'elink4'        => {'tests' => 28,
        #                    'sub'   => \&elink4},
        #'multilink1'    => {'tests' => 40,
        #                    'sub'   => \&multilink1},
        #'multilink2'    => {'tests' => 49,
        #                    'sub'   => \&multilink2},
        #'multilink3'    => {'tests' => 0,
        #                    'sub'   => \&multilink3},
        #'scores'        => {'tests' => 0,
        #                    'sub'   => \&scores},
        );
    $NUMTESTS += $EUTILS{$_}->{'tests'} for (keys %EUTILS);
    $DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
    # this seems to work for perl 5.6 and perl 5.8
    eval {require Test::More;};
    
    if ($@) {
        use lib 't/lib';
    }
    
    use Test::More;
	use BioperlTest;
	
	test_begin(-tests               => $NUMTESTS,
			   -requires_modules    => [qw(XML::Simple LWP::UserAgent)],
			   -requires_networking => 1,
			  );
    
    use_ok('Bio::DB::EUtilities');
    use_ok('LWP::UserAgent');
    use_ok('XML::Simple');
}

# NOTE : Bio::DB::EUtilities is just a specialized pipeline to get any 
# data available via NCBI's Entrez interface, with a few convenience methods
# to get UIDs and other additional information.  All data returned
# using EFetch is raw (not Bioperl objects) and is meant to be piped into
# other Bioperl modules at a later point for further processing

#   protein acc
my @acc = qw(MUSIGHBA1 P18584 CH402638);

# protein GI
my @ids = sort qw(1621261 89318838 68536103 20807972 730439);

# test search term
my $term = 'dihydroorotase AND human';

my ($eutil, $response);

my %dbs = (taxonomy => 1,
           nucleotide => 1,
           pubmed => 1);
my %links = (protein_taxonomy => 1,
             protein_nucleotide => 1,
             protein_nucleotide_wgs => 1,
             protein_pubmed => 1,
             protein_pubmed_refseq => 1
             );

# this loops through the required tests, only running what is in %EUTILS
for my $test (keys %EUTILS) {
    $EUTILS{$test}->{'sub'}->();
}

# Simple EFetch

sub efetch {
    SKIP: {
        $eutil = Bio::DB::EUtilities->new(
                                        -db         => 'protein',
                                        -id         => [$ids[0]],
                                        -rettype    => 'fasta'
                                          );
              
        isa_ok($eutil, 'Bio::DB::GenericWebDBI');
        eval {$response = $eutil->get_response; };
        skip("EFetch HTTP error: $@", 4) if $@;
        isa_ok($response, 'HTTP::Response');
        my $content = $response->content;
        like($content, qr(PYRR \[Mycobacterium tuberculosis H37Rv\]),
             'EFetch: Fasta format');
        
        # reuse the EUtilities webagent
        $eutil->parameter_base->id([$ids[1]]);
        $eutil->parameter_base->rettype('gb');
        eval {$response = $eutil->get_response; };
        skip("EFetch HTTP error: $@", 2) if $@;
        isa_ok($response, 'HTTP::Response');
        $content = $response->content;
        like($content, qr(^LOCUS\s+NP_623143),'EFetch: GenBank format');
    }
}

# EPost->EFetch with History (Cookie)

sub epost {
    SKIP: {
        $eutil = Bio::DB::EUtilities->new(
                                        -eutil      => 'epost',
                                        -db         => 'protein',
                                        -id         => \@ids,
                                          );
              
        isa_ok($eutil, 'Bio::DB::GenericWebDBI');
        eval {$response = $eutil->get_response; };
        skip("EPost HTTP error", 12) if $@;
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
        my $total;
        eval{ $total = grep(m{^>.*$}, split "\n", $efetch->get_response->content);};
        skip("EPost HTTP error", 1) if $@;
        is($total, 5, 'EPost to EFetch');
    }
}

# ESummary

sub esummary {
    my %docsum = (1621261=> { 'Caption' => ['String','CAB02640'],
    'Title' => ['String','PROBABLE PYRIMIDINE OPERON REGULATORY PROTEIN PYRR '.
     '[Mycobacterium tuberculosis H37Rv]'],
    'Extra' => ['String','gi|1621261|emb|CAB02640.1|[1621261]'],
    'Gi' => ['Integer','1621261'],
    'CreateDate' => ['String','2003/11/21'],
    'UpdateDate' => ['String','2005/04/17'],
    'Flags' => ['Integer',''],
    'TaxId' => ['Integer','83332'],
    'Length' => ['Integer','193'],
    'Status' => ['String','live'],
    'ReplacedBy' => ['String',''],
    'Comment' => ['String',''], },
    20807972 => {'Caption' => ['String','NP_623143'],
    'Title' => ['String','pyrimidine regulatory protein PyrR '.
     '[Thermoanaerobacter tengcongensis MB4]'],
    'Extra' => ['String','gi|20807972|ref|NP_623143.1|[20807972]'],
    'Gi' => ['Integer','20807972'],
    'CreateDate' => ['String','2002/05/09'],
    'UpdateDate' => ['String','2005/12/03'],
    'Flags' => ['Integer','512'],
    'TaxId' => ['Integer','273068'],
    'Length' => ['Integer','178'],
    'Status' => ['String','live'],
    'ReplacedBy' => ['String',''],
    'Comment' => ['String',''], },
    68536103 => {'Caption' => ['String','YP_250808'],
    'Title' => ['String','putative pyrimidine operon regulatory protein '.
     '[Corynebacterium jeikeium K411]'],
    'Extra' => ['String','gi|68536103|ref|YP_250808.1|[68536103]'],
    'Gi' => ['Integer','68536103'],
    'CreateDate' => ['String','2005/07/04'],
    'UpdateDate' => ['String','2006/03/30'],
    'Flags' => ['Integer','512'],
    'TaxId' => ['Integer','306537'],
    'Length' => ['Integer','195'],
    'Status' => ['String','live'],
    'ReplacedBy' => ['String',''],
    'Comment' => ['String',''], },
    730439 => {'Caption' => ['String','P41007'],
    'Title' => ['String','PyrR bifunctional protein '.
     '[Includes: Pyrimidine operon regulatory protein; '.
     'Uracil phosphoribosyltransferase (UPRTase)]'],
    'Extra' => ['String','gi|730439|sp|P41007|PYRR_BACCL[730439]'],
    'Gi' => ['Integer','730439'],
    'CreateDate' => ['String','1995/02/01'],
    'UpdateDate' => ['String','2006/07/25'],
    'Flags' => ['Integer',''],
    'TaxId' => ['Integer','1394'],
    'Length' => ['Integer','179'],
    'Status' => ['String','live'],
    'ReplacedBy' => ['String',''],
    'Comment' => ['String',''] },
    89318838 => { 'Caption' => ['String','EAS10332'],
    'Title' => ['String','Phosphoribosyltransferase '.
     '[Mycobacterium gilvum PYR-GCK]'],
    'Extra' => ['String','gi|89318838|gb|EAS10332.1|[89318838]'],
    'Gi' => ['Integer','89318838'],
    'CreateDate' => ['String','2006/03/09'],
    'UpdateDate' => ['String','2006/03/09'],
    'Flags' => ['Integer',''],
    'TaxId' => ['Integer','350054'],
    'Length' => ['Integer','193'],
    'Status' => ['String','live'],
    'ReplacedBy' => ['String',''],
    'Comment' => ['String',''] } );
    SKIP: {
        $eutil = Bio::DB::EUtilities->new(
                                         -eutil      => 'esummary',
                                         -db         => 'protein',
                                         -id            => \@ids,
                                           );
        isa_ok($eutil, 'Bio::DB::GenericWebDBI');
        
        eval {$response = $eutil->get_response; };
        skip("ESummary HTTP error:$@", 253) if $@;
        isa_ok($response, 'HTTP::Response');
        
        my @docs = $eutil->get_all_docsums();
        is(scalar(@docs), 5, '$esum->get_all_docsums()');
        
        my $ct = 0;
        while (my $ds = $eutil->next_docsum) {
            isa_ok($ds, 'Bio::DB::EUtilities::DocSum');
            
            my $id = $ds->esummary_id();
            ok(exists($docsum{$id}), '$docsum->esummary_id()');
            
            my %items = %{ $docsum{$id} };
            
            # iterate using item names
            
            for my $name ($ds->get_all_names()) {
                $ct++;
                my %data = $ds->get_item_by_name($name);
                ok(exists $items{$name},'DocSum Name exists');
                is($data{Name}, $name, 'get_item_by_name(),DocSum Name');
                is($ds->get_Type_by_name($name), $items{$name}->[0],
                   'get_Type_by_name()');
                is($data{Type}, $items{$name}->[0], 'get_item_by_name(),DocSum Type');
            }
        }
        is($ct, 60);
    }
}

# ESearch, ESearch History

sub esearch {
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
}

# EInfo

sub einfo {
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
        cmp_ok(scalar(@fields), '>',24, '$einfo->einfo_dbfield_info()');
    
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
}


# ELink - normal (single ID array) - single db - ElinkData tests

sub elink1 {
    SKIP: {
        $eutil = Bio::DB::EUtilities->new(
                                        -eutil      => 'elink',
                                        -db         => 'taxonomy',
                                        -dbfrom     => 'protein',
                                        -id         => \@ids,
                                          );
              
        isa_ok($eutil, 'Bio::DB::GenericWebDBI');
        eval {$response = $eutil->get_response; };
        skip("ELink HTTP error:$@", 8) if $@;
        isa_ok($response, 'HTTP::Response');
        like($response->content, qr(<eLinkResult>), 'ELink response');
        # Data is too volatile to test; commenting for now...
        #my @ids2 = qw(350054 306537 273068 83332 1394);
        cmp_ok(@{$eutil->get_ids}, '>=', 4);
        #is_deeply([sort $eutil->get_ids], [sort @ids2],'$elink->get_ids()');
        
        # Linkset tests
        is($eutil->get_linkset_count, 1, '$elink->get_linkset_count()');
        my $linkobj = $eutil->next_linkset;
        isa_ok($linkobj, 'Bio::DB::EUtilities::ElinkData');
        is($linkobj->elink_dbfrom, 'protein', '$linkdata->elink_dbfrom()');
        #is_deeply([sort $linkobj->elink_queryids],
        #          [sort @ids], '$linkdata->elink_queryids()');
        is($linkobj->elink_command, 'neighbor', '$linkdata->elink_command()');
        my $db = $linkobj->next_linkdb;
        is($db, 'taxonomy', '$linkdata->next_linkdb()');
        #is_deeply([sort $linkobj->get_LinkIds_by_db($db)],
        #          [sort @ids2], '$linkdata->get_LinkIds_by_db($db)');   
    }
}

# ELink - normal (single ID array), multiple dbs 

sub elink2 {
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
        skip("ELink HTTP error:$@", 17) if $@;
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
        
        my $ct = 4;
        # tests per iteration
        my $ti = 2;
        while (my $db = $linkobj->next_linkdb) {
            ok(exists $ids{$db}, "ElinkData database: $db");
            @ids2 = sort $linkobj->get_LinkIds_by_db($db);
            is_deeply($ids{$db}, \@ids2, "ElinkData database IDs: $db");
        }
        is($ct, 0);
        # other ElinkData methods
        is($linkobj->elink_dbfrom, 'protein', '$linkdata->elink_dbfrom()');
        is_deeply([sort $linkobj->elink_queryids],
                  [sort @ids], '$linkdata->elink_queryids()');
        is($linkobj->elink_command, 'neighbor', '$linkdata->elink_command()');
        skip('No elink data: possible server problem',$ct*$ti);
    }
}

# ELink - normal (single ID array), multiple dbs, cookies)

sub elink3 {
    SKIP: {
        # can use 'all' for db, but takes a long time; use named dbs instead
        # this retrieves cookies instead (no ElinkData objects are stored)
        $eutil = Bio::DB::EUtilities->new(
                                        -eutil      => 'elink',
                                        -db         => 'taxonomy,nucleotide,pubmed',
                                        -dbfrom     => 'protein',
                                        -id         => \@ids,
                                        -cmd        => 'neighbor_history',
                                          );
        isa_ok($eutil, 'Bio::DB::GenericWebDBI');
        eval {$response = $eutil->get_response; };
        skip("ELink HTTP error", 27) if $@;
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
        my $ct = 2;
        my $content;
        
        # tests per iteration
        my $ti = 10;
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
                $content = $efetch->get_response->content;
            }
            $ct--;
            last if $ct == 0;
        }
        like($content, qr(<TaxaSet>), 'ELink to EFetch : taxonomy');
        is($ct,0,'Cookie count');
        skip('No cookies returned; possible server problem',$ct*$ti);
    }
}

# ELink (multi_id), single db
# this is a flag set to get one-to-one correspondence for ELink data

sub elink4 {
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
        
        skip("ELink HTTP error", 27) if $@;
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
        my $ct = 5;
        my $ti = 4;
        # ids may not be returned in same order as array, so need to grab and sort
        while ( my $linkobj = $eutil->next_linkset) {
            isa_ok($linkobj, 'Bio::DB::EUtilities::ElinkData');
            is($linkobj->elink_dbfrom, 'protein', '$linkdata->elink_dbfrom()');
            is($linkobj->elink_command, 'neighbor', '$linkdata->elink_command()');
            push @qids, $linkobj->elink_queryids;
            my $db = $linkobj->next_linkdb;
            is($db, 'taxonomy', '$linkdata->next_linkdb()');
            push @retids, $linkobj->get_LinkIds_by_db($db);
            $ct--;
            last if $ct == 0;
        }
        is($ct,0);
        is_deeply([sort @qids], [sort @ids], '$linkdata->elink_queryids()');
        is_deeply([sort @retids], [sort @ids2], '$linkdata->get_LinkIds_by_db($db)');
        skip('No Elink data returned; possible server problem',$ct*$ti);
    }
}

# ELink (multi_id, cookies)
# these need to be cleaned up

sub multilink1 {
    SKIP: {
        $eutil = Bio::DB::EUtilities->new(
                                        -eutil      => 'elink',
                                        -db         => 'taxonomy',
                                        -dbfrom     => 'protein',
                                        -multi_id   => 1,
                                        -id         => \@ids,
                                        -cmd        => 'neighbor_history',
                                        -verbose    => 2
                                          );
              
        isa_ok($eutil, 'Bio::DB::GenericWebDBI');
        eval {$response = $eutil->get_response; };
        
        # check this number, likely wrong
        skip("ELink HTTP error", 39) if $@;
        isa_ok($response, 'HTTP::Response');
        like($response->content, qr(<eLinkResult>), 'ELink response');
        my @ids2 = qw(350054 306537 273068 83332 1394);
        
        # This is designed to fail; IDs present in individual ElinkData objects
        # for one-to-one correspondence with ID groups
        eval{$eutil->get_ids;};
        ok($@,'$elink->get_ids()');
        
        # Linkset tests (there aren't any)
        is($eutil->get_linkset_count, 0, '$elink->get_linkset_count()');
        cmp_ok($eutil->get_cookie_count, '>=', 4, '$elink->get_cookie_count()');
        
        my $efetch = Bio::DB::EUtilities->new();
        my $ct = 2;
        my $content;
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
                $content = $efetch->get_response->content;
                like($content, qr(<TaxaSet>), 'ELink to EFetch : taxonomy');
            }
            last if $ct == 0;
            $ct--;
        }
        is($ct, 0);
    }
}

# ELink (multi_id, multidbs)

sub multilink2 {
    SKIP: {
        $eutil = Bio::DB::EUtilities->new(
                                        -eutil      => 'elink',
                                        -db         => 'taxonomy,nucleotide,pubmed',
                                        -dbfrom     => 'protein',
                                        -multi_id   => 1,
                                        -id         => \@ids,
                                        -verbose    => 2
                                          );
              
        isa_ok($eutil, 'Bio::DB::GenericWebDBI');
        eval {$response = $eutil->get_response; };
        
        # check this number, likely wrong
        skip("ELink HTTP error", 20) if $@;
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
        my $ct = 4;
        
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
            $ct--;
            last if $ct == 0;
        }
        is($ct, 0);
        skip();
    }
}

# ELink (multi_id, multidb, cookies)

sub multilink3 {
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
}

# ELink (scores)

sub scores {
    my %scores = (   1621261 =>   2147483647,
                    20807972 =>          423,
                    68536103 =>          554,
                      730439 =>          411,
                    89318838 =>          '',);

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
        skip("ELink HTTP error:$@", 20);# if $@;
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
                my $ct = 0;
                is($db,'protein', '$linkdata->next_scoredb()');
                my @ids2 = $linkobj->get_LinkIds_by_db($db);
                cmp_ok(scalar(@ids2), '>', 765, '$linkdata->get_LinkIds_by_db($db)');
                for my $id (@ids) {
                    last if $ct++ == 6;
                    is($linkobj->get_score($id), $scores{$id}, '$linkdata->get_score()');
                }
            }
        }
    }
}

# Although the other EUtilities are available, no postprocessing is done on the
# returned XML yet

sub egquery {
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

1;
