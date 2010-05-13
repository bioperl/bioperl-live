# -*-Perl-*- Test Harness script for Bioperl
# $Id: EUtilities.t 15112 2008-12-08 18:12:38Z sendu $
#

use strict;
our $NUMTESTS;
our $DEBUG;
our %EUTILS;

BEGIN {
    $NUMTESTS = 4; # base number of tests (those not in blocks)

    # I have set up eutils tests to run in sections for easier test maintenance
    # and keeping track of problematic tests. The below hash is the list of
    # tests, with test number and coderef.
    
    # these now run very simple tests for connectivity and data sampling
    # main tests now with the parser

    %EUTILS = (
        'efetch'        => {'tests' => 5,
                            'sub'   => \&efetch},
        'epost'         => {'tests' => 11,
                            'sub'   => \&epost},
        'esummary'      => {'tests' => 254,
                            'sub'   => \&esummary},
        'esearch'       => {'tests' => 13,
                            'sub'   => \&esearch},
        'einfo'         => {'tests' => 10,
                            'sub'   => \&einfo},
        'elink1'        => {'tests' => 8,
                            'sub'   => \&elink1},
        'egquery'       => {'tests' => 4,
                            'sub'   => \&egquery},
        );
    $NUMTESTS += $EUTILS{$_}->{'tests'} for (keys %EUTILS);
    $DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
    # this seems to work for perl 5.6 and perl 5.8

	use Bio::Root::Test;
	
	test_begin(-tests               => $NUMTESTS,
			   -requires_modules    => [qw(XML::Simple LWP::UserAgent)],
			   -requires_email      => 1,
			  );
    
    use_ok('Bio::DB::EUtilities');
    use_ok('LWP::UserAgent');
    use_ok('Bio::Tools::EUtilities');
    use_ok('Bio::Tools::EUtilities::EUtilParameters');
}

my $email = test_email();

diag("Using $email for tests") if $DEBUG;

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
                                        -rettype    => 'fasta',
                                        -email      => $email
                                          );
              
        isa_ok($eutil, 'Bio::DB::GenericWebAgent');
        eval {$response = $eutil->get_Response; };
        skip("EFetch HTTP error: $@", 4) if $@;
        isa_ok($response, 'HTTP::Response');
        my $content = $response->content;
        like($content, qr(PYRR \[Mycobacterium tuberculosis H37Rv\]),
             'EFetch: Fasta format');
        
        # reuse the EUtilities webagent
        $eutil->parameter_base->id([$ids[1]]);
        $eutil->parameter_base->rettype('gb');
        eval {$response = $eutil->get_Response; };
        skip("EFetch HTTP error: $@", 2) if $@;
        isa_ok($response, 'HTTP::Response');
        $content = $response->content;
        like($content, qr(^LOCUS\s+NP_623143),'EFetch: GenBank format');
    }
}

# EPost->EFetch with History

sub epost {
    SKIP: {
        $eutil = Bio::DB::EUtilities->new(
                                        -eutil      => 'epost',
                                        -db         => 'protein',
                                        -id         => \@ids,
                                        -email      => $email
                                          );
              
        isa_ok($eutil, 'Bio::DB::GenericWebAgent');
        eval {$response = $eutil->get_Response; };
        skip("EPost HTTP error: $@", 10) if $@;
        isa_ok($response, 'HTTP::Response');
        # Any parameters are passed in to the parser, so these should be set.
        # Databases and IDs always default back to the submitted ones unless
        # the data being retrieved are IDs or contain new IDs (esearch, elink)
        
        is($eutil->get_database, 'protein', '$epost->get_database()');
        is(join(',',$eutil->get_ids), '1621261,20807972,68536103,730439,89318838', '$epost->get_ids()');
        
        # these are not set using epost
        is($eutil->get_count, undef, '$epost->get_count()');
        is($eutil->get_term, undef, '$epost->get_term()');

        my $history = $eutil->next_History;
        is($history->eutil, 'epost', 'History->eutil()');
        isa_ok($history, 'Bio::Tools::EUtilities::HistoryI');
        
        # check the actual History
        my ($webenv, $key) = $history->history;
        like($webenv, qr{^\S{25}}, '$epost WebEnv');
        like($key, qr{^\d+}, '$epost query key');
        
        # can we fetch the sequences?
        $eutil->set_parameters(
            -eutil => 'efetch',
            -history     => $history,
            -rettype    => 'fasta'
        );
        # look for fasta headers
        my ($r, $t);
        eval{ $r = $eutil->get_Response->content;};
        skip("EPost HTTP error", 1) if $@;
        $t = grep m{^>.*$}, split("\n", $r);
        is($t, 5, 'EPost to EFetch');
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
                                         -email      => $email
                                           );
        isa_ok($eutil, 'Bio::DB::GenericWebAgent');
        
        eval {$response = $eutil->get_Response; };
        skip("ESummary HTTP error:$@", 253) if $@;
        isa_ok($response, 'HTTP::Response');
        
        my @docs = $eutil->get_DocSums();
        is(scalar(@docs), 5, '$esum->get_DocSums()');
        
        my $ct = 0;
        while (my $ds = $eutil->next_DocSum) {
            isa_ok($ds, 'Bio::Tools::EUtilities::Summary::DocSum');
            
            my $id = $ds->get_id();
            ok(exists($docsum{$id}), '$docsum->get_id()');
            
            my %items = %{ $docsum{$id} };
            
            # iterate using item names
            
            for my $name ($ds->get_all_names()) {
                $ct++;
                my ($it) = $ds->get_Items_by_name($name);
                ok(exists $items{$name},'DocSum Name exists');
                is($it->get_name, $name, 'get_name(),DocSum Name');
                is($ds->get_type_by_name($name), $items{$name}->[0],
                   'get_type_by_name() from DocSum');
                is($it->get_type, $items{$name}->[0], 'get_type() from Item');
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
                                        -retmax     => 100,
                                        -email      => $email
                                          );
              
        isa_ok($eutil, 'Bio::DB::GenericWebAgent');
        eval {$response = $eutil->get_Response; };
        skip("ESearch HTTP error:$@", 12) if $@;
        isa_ok($response, 'HTTP::Response');
        
        # can't really check for specific ID's but can check total ID's returned
        my @esearch_ids = $eutil->get_ids;
        is(scalar(@esearch_ids), 100, '$esearch->get_ids()');
        
        cmp_ok($eutil->get_count, '>', 117, '$esearch->get_count()');
    
        # usehistory
        $eutil = Bio::DB::EUtilities->new(
                                        -eutil      => 'esearch',
                                        -db         => 'protein',
                                        -usehistory => 'y',
                                        -term       => $term,
                                        -retmax     => 100,
                                        -email      => $email
                                          );
        
        eval {$response = $eutil->get_Response; };
        skip("ESearch HTTP error:$@", 9) if $@;
        is($eutil->eutil, 'esearch', 'eutil()');
        is($eutil->get_database, 'protein', 'get_database()');
        cmp_ok($eutil->get_count, '>', 117, 'get_count()');
        is($eutil->get_term, $term, 'get_term()');
        is($eutil->get_ids, 100, 'History->get_ids()');
        
        my $history = $eutil->next_History;
        isa_ok($history, 'Bio::Tools::EUtilities::HistoryI');
        
        # check the actual data
        my ($webenv, $key) = $history->history;
        like($webenv, qr{^\S{15}}, 'WebEnv');
        like($key, qr{^\d+}, 'query key');
        
        # can we fetch the sequences?
        $eutil->set_parameters(
            -eutil      => 'efetch',
            -history    => $history,
            -rettype    => 'fasta',
            -retmax     => 5
        );
        # look for fasta headers
        my ($r, $t);
        eval{ $r = $eutil->get_Response->content;};
        skip("EPost HTTP error", 1) if $@;
        $t = grep m{^>.*$}, split("\n", $r);
        is($t, 5, 'EPost to EFetch');
    }
}

# EInfo

sub einfo {
    SKIP: {
        $eutil = Bio::DB::EUtilities->new(
                                        -eutil      => 'einfo',
                                        -db         => 'protein',
                                        -email      => $email
                                          );
        isa_ok($eutil, 'Bio::DB::GenericWebAgent');
        eval {$response = $eutil->get_Response; };
        skip("EInfo HTTP error:$@", 10) if $@;
        isa_ok($response, 'HTTP::Response');
        like($response->content, qr(<eInfoResult>), 'EInfo response');
        is(($eutil->get_database)[0], 'protein', '$einfo->get_database()');
        like($eutil->get_last_update, qr(\d{4}\/\d{2}\/\d{2}\s\d{2}:\d{2}),
             '$einfo->get_last_update()');
        cmp_ok($eutil->get_record_count, '>', 9200000, '$einfo->get_record_count()');
        is($eutil->get_description, 'Protein sequence record', '$einfo->get_description()');
        my @links = $eutil->get_LinkInfo;
        my @fields = $eutil->get_FieldInfo;
        cmp_ok(scalar(@links), '>',30, '$einfo->get_LinkInfo()');
        cmp_ok(scalar(@fields), '>',24, '$einfo->get_FieldInfo()');
    
        # all databases (list)
        $eutil = Bio::DB::EUtilities->new(
                                        -eutil      => 'einfo',
                                        -email      => $email
                                          );
        
        eval {$response = $eutil->get_Response; };
        skip("EInfo HTTP error:$@", 1) if $@;
        
        my @db = sort qw(pubmed  protein  nucleotide  nuccore  nucgss  nucest  structure
        genome  books  cancerchromosomes  cdd  domains  gene  genomeprj  gensat
        geo  gds  homologene  journals  mesh  ncbisearch  nlmcatalog  omia  omim
        pmc  popset  probe  pcassay  pccompound  pcsubstance  snp  taxonomy toolkit
        unigene  unists);
        
        my @einfo_dbs = sort $eutil->get_databases;
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
                                        -email      => $email
                                          );
              
        isa_ok($eutil, 'Bio::DB::GenericWebAgent');
        eval {$response = $eutil->get_Response; };
        skip("ELink HTTP error:$@", 7) if $@;
        isa_ok($response, 'HTTP::Response');
        like($response->content, qr(<eLinkResult>), 'ELink response');
        # Data is too volatile to test; commenting for now...
        #my @ids2 = qw(350054 306537 273068 83332 1394);
        cmp_ok($eutil->get_ids, '>=', 4);
        #is_deeply([sort $eutil->get_ids], [sort @ids2],'$elink->get_ids()');
        
        # Linkset tests
        is($eutil->get_LinkSets, 1, '$elink->get_LinkSets()');
        my $linkobj = $eutil->next_LinkSet;
        isa_ok($linkobj, 'Bio::Tools::EUtilities::Link::LinkSet');
        is($linkobj->get_dbfrom, 'protein', '$linkdata->get_dbfrom()');
        #is_deeply([sort $linkobj->elink_queryids],
        #          [sort @ids], '$linkdata->elink_queryids()');
        my $db = $linkobj->get_dbto;
        is($db, 'taxonomy', '$linkdata->get_dbto()');
        #is_deeply([sort $linkobj->get_LinkIds_by_db($db)],
        #          [sort @ids2], '$linkdata->get_LinkIds_by_db($db)');   
    }
}

sub elink2 {
    my @genome_ids = qw(30807 33011 12997 16707 45843 31129 31141 31131 31133 32203 31135);
    SKIP: {
        $eutil = Bio::DB::EUtilities->new(
                                        -eutil      => 'elink',
                                        -db         => 'nuccore',
                                        -dbfrom     => 'genomeprj',
                                        -id         => @genome_ids,
                                        -email      => $email
                                          );
              
        eval {$response = $eutil->get_Response; };
        skip("ELink HTTP error:$@", 7) if $@;
        isa_ok($response, 'HTTP::Response');
        like($response->content, qr(<eLinkResult>), 'ELink response');
        # Data is too volatile to test; commenting for now...
        #my @ids2 = qw(350054 306537 273068 83332 1394);
        cmp_ok($eutil->get_ids, '>=', 4);
        #is_deeply([sort $eutil->get_ids], [sort @ids2],'$elink->get_ids()');
        
        # Linkset tests
        is($eutil->get_LinkSets, 1, '$elink->get_LinkSets()');
        my $linkobj = $eutil->next_LinkSet;
        isa_ok($linkobj, 'Bio::Tools::EUtilities::Link::LinkSet');
        is($linkobj->get_dbfrom, 'protein', '$linkdata->get_dbfrom()');
        #is_deeply([sort $linkobj->elink_queryids],
        #          [sort @ids], '$linkdata->elink_queryids()');
        my $db = $linkobj->get_dbto;
        is($db, 'taxonomy', '$linkdata->get_dbto()');
        #is_deeply([sort $linkobj->get_LinkIds_by_db($db)],
        #          [sort @ids2], '$linkdata->get_LinkIds_by_db($db)');   
    }
}

sub egquery {
    SKIP: {
    $eutil = Bio::DB::EUtilities->new(
                                    -eutil      => 'egquery',
                                    -term       => $term,
                                    -email      => $email
                                      );
          
    isa_ok($eutil, 'Bio::DB::GenericWebAgent');
    eval {$response = $eutil->get_Response; };
    skip("EGQuery HTTP error:$@", 3) if $@;
    isa_ok($response, 'HTTP::Response');
    like($response->content, qr(<eGQueryResult>), 'EGQuery response');
    my @gq = $eutil->get_GlobalQueries;
    cmp_ok(scalar(@gq), '>=', 30, 'get_GlobalQueries')
    }
}

1;
