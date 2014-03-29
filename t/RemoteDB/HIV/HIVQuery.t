# testing Bio::DB::Query::HIVQuery
# $Id: HIVQuery.t 232 2008-12-11 14:51:51Z maj $
use strict;
use warnings;

BEGIN {
    use Bio::Root::Test;
    test_begin(
        -tests => 41,
        -requires_modules => [qw( XML::Simple )]
    );
    use_ok('Bio::DB::Query::HIVQuery');
    use_ok('Bio::DB::HIV');
    use_ok('Bio::Annotation::Collection');
    use_ok('Bio::Annotation::Comment');
    use_ok('Bio::Annotation::Reference');
    use_ok('Bio::DB::HIV::HIVQueryHelper');

}

my $tobj= new Bio::DB::Query::HIVQuery(-RUN_OPTION=>0);

#object tests
isa_ok($tobj, 'Bio::DB::Query::HIVQuery');

#compliance tests
isa_ok($tobj, 'Bio::Root::Root');
can_ok($tobj, qw( count ids query ));

#methods
can_ok($tobj, qw( 
                  get_annotations_by_ids 
                  get_annotations_by_id 
                  add_annotations_for_id 
                  remove_annotations_for_ids
                  remove_annotations_for_id
                  remove_annotations
                  get_accessions
                  get_accessions_by_ids
                  get_accessions_by_id
                  )
    );
#internals
can_ok($tobj, qw(
                  _do_query
                  _reset
                  _session_id
                  _run_option
                  add_id
                  lanl_base
                  map_db
                  make_search_if
                  search_
                  _map_db_uri
                  _make_search_if_uri
                  _search_uri
                  _schema_file
                  _schema
                  _lanl_query
                  _lanl_response
                  _create_lanl_query
                  _do_lanl_request
                  _parse_lanl_response
                  _parse_query_string
                  _sorry 
                 )
    );

#defaults tests
ok($tobj->_map_db_uri, "_map_db_uri set in default object");
ok($tobj->_make_search_if_uri, "_make_search_if_uri set in default object");
ok($tobj->_search_uri, "_search_uri set in default object");
ok($tobj->_schema_file, "_schema_file set in default object");
ok(defined $tobj->_run_option, "_run_option set in default object");
ok($tobj->{_annotations}, "annotations container available");

# query syntax (no run)
$tobj->_reset;
$tobj->query(['coreceptor'=>['CCR5 CXCR4','CXCR4'], subtype=>['D', 'F'], country=>'Any']);
is($tobj->_do_query(0), 0, 'query syntax check 1');
$tobj->_reset;
$tobj->query({'coreceptor'=>['CCR5 CXCR4','CXCR4'], subtype=>['D', 'F'], country=>'Any'});
is($tobj->_do_query(0), 0, 'query syntax check 2');
$tobj->_reset;
$tobj->query({'query'=>{'coreceptor'=>['CCR5 CXCR4','CXCR4'], subtype=>['D', 'F']}, 'annot'=> ['country']});
is($tobj->_do_query(0),0, 'query syntax check 3');
$tobj->_reset;
$tobj->query("('CCR5 CXCR4', 'CXCR4')[coreceptor] (D F)[subtype] {[country]}");
is($tobj->_do_query(0),0, 'query parser check');

#multiquery parse
$tobj->_reset;
$tobj->query( "(SI[phenotype] ('CCR5 CXCR4')[coreceptor] C[subtype] OR NSI[phenotype] D[subtype]) AND ZA[country]");
is($tobj->_do_query(0),0, 'multiquery parse check');

SKIP: {
	test_skip(-requires_module => 'HTML::Parser', -tests => 3);
	no warnings; # HTML::Parser has an uninited value issue
	use_ok('HTML::Parser');
	use warnings;
	# help test; just tests that file can be written and that tags are in matching 
	# pairs, with reasonable placement of <html> and </html>
	my $hlpf = test_output_file();
	my $a = 0;
	my $html = 0;
	my $h = HTML::Parser->new(empty_element_tags => 1,
				  start_h => [sub {$a++; (shift eq 'html') && ($a==1) && $html++;}, "tagname"],
				  end_h   => [sub {$a--; (shift eq 'html') && ($a==0) && $html++;}, "tagname"] );
	ok($tobj->help($hlpf), "help html to file");
	{
		local $/;
		undef $/;
		open my $HP, '<', $hlpf or die "Could not read file '$hlpf': $!\n";
		$h->parse(<$HP>);
		is_deeply([$a, $html], [0, 2], "help html parsed");
		close $HP; # Always explicitly close filehandles
		1;
	}

}

#exceptions tests
#pre-run query exceptions
$tobj->verbose(2);
$tobj->_reset;
$tobj->query( "narb[scroob]" );
throws_ok {$tobj->_do_query} qr/BadParameter/, "bad field exception check";
$tobj->_reset;
$tobj->query( "narb[phenotype]" );
throws_ok {$tobj->_do_query} qr/BadParameter/, "bad match data exception check";
$tobj->_reset;
$tobj->query( [ 'cd4_count'  => "" ] );
throws_ok {$tobj->_do_query} qr/BadParameter/, "empty field not ok exception check";
$tobj->_reset;
$tobj->{_schema} = undef;
throws_ok {$tobj->_do_query} qr/SchemaNotInit/, "uninitialized schema exception check";
throws_ok {$tobj->count} qr/Query not yet run/, "query not run (level 1) warning check";
throws_ok {$tobj->get_annotations_by_id($tobj->ids)} qr/Requires query run/, "query not run (level 2) warning check";



# network tests
SKIP : {
    test_skip(-tests => 10,
          -requires_networking => 1);
    eval {$tobj  = Bio::DB::Query::HIVQuery->new(-QUERY=>"(SI[phenotype] ('CCR5 CXCR4')[coreceptor] C[subtype] OR NSI[phenotype] D[subtype]) AND ZA[country]",-RUN_OPTION=>2)};
    if ($@) {
        diag($@);
        skip("Network problems, skipping all", 10);
    }
    ok($tobj,"live query");
    cmp_ok( $tobj->count, ">=", 12, "enough sequences returned");
# test query object handling of Bio::DB::HIV   
    my ($tdb, $seqio, $seq);
    ok( $tdb = new Bio::DB::HIV, "create Bio::DB::HIV object");
        eval {$seqio = $tdb->get_Stream_by_query($tobj)};
    if ($@) {
        diag($@);
        skip("Network problems, skipping all", 7);
    }
    ok($seqio, "get SeqIO stream from query");
# test HIVAnnotProcessor indirectly
    ok($seq = $seqio->next_seq, "access sequence stream");
    ok($seq->annotation->get_value('Virus'), "'Virus' annotation present");
    like ($seq->annotation->get_value('Virus','phenotype'), qr/SI/, "'Virus' phenotype annotation present");
    like ($seq->annotation->get_value('Virus', 'subtype'), qr/[CD]/, "'Virus' subtype annotation present");
    ok($seq->accession_number, "GenBank accession available");

# exception checks
    $tobj->_reset;
    $tobj->verbose(2);
    $tobj->query( "2BR02B[accession]" );
    throws_ok {$tobj->_do_query(2)} qr/no sequences/i, "no sequences warning check";
}
