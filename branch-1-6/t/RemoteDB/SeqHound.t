# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 15,
			   -requires_modules => [qw(IO::String LWP::UserAgent)],
			   -requires_networking => 1);
	
	use_ok('Bio::DB::SeqHound');
}

END {
	unlink $Bio::DB::SeqHound::LOGFILENAME if -f $Bio::DB::SeqHound::LOGFILENAME;
}

my $verbose = test_debug();

my ($db,$seq,$seqio);
# get a single seq

$seq = $seqio = undef;

SKIP: {
    $db = Bio::DB::SeqHound->new(-verbose=>$verbose);
    eval {ok(defined($seq = $db->get_Seq_by_acc('J00522')));};
	skip('Could not connect to seqhound, skipping tests', 10) if $@;
    is( $seq->length, 408); 
    ok defined ($db->request_format('fasta'));
    eval {ok(defined($seq = $db->get_Seq_by_acc('NP_862707')));};
	skip('Could not connect to seqhound, skipping tests', 7) if $@;
    is( $seq->accession, 'NP_862707');
    is( $seq->length, 227); 
    ok( defined($db = Bio::DB::SeqHound->new(-verbose=>$verbose, 
					-retrievaltype => 'tempfile')));
    eval {ok(defined($seqio = $db->get_Stream_by_id(['BTACHRE'])));};
	skip('Could not connect to seqhound, skipping tests', 3) if $@;
    undef $db; # testing to see if we can remove db
    ok( defined($seq = $seqio->next_seq()));
    is( $seq->length, 1621);
}

$seq = $seqio = undef;

SKIP: {
    $db = Bio::DB::SeqHound->new(-verbose => $verbose,
			    -retrievaltype => 'tempfile',
			    -format => 'genbank'
			    ); 
    eval {ok( defined($seqio = $db->get_Stream_by_acc(['J00522', 'AF303112', 'J02231'])));};
	skip('Could not connect to seqhound for batch test, skipping tests', 4) if $@;
	is($seqio->next_seq->length, 408);
    is($seqio->next_seq->length, 1611);
    is($seqio->next_seq->length, 200);
}
