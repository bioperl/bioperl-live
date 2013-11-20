# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 16,
			   -requires_modules => [qw(IO::String HTTP::Request::Common)],
			   -requires_networking => 1);
	
	use_ok('Bio::DB::EMBL');
}

my $verbose = test_debug();

my ($db,$seq,$seqio);
# get a single seq

$seq = $seqio = undef;

SKIP: { 
    ok defined($db = Bio::DB::EMBL->new(-verbose=>$verbose)); 
    ok(defined($seq = $db->get_Seq_by_acc('J00522')));
    is( $seq->length, 408); 
    ok defined ($db->request_format('fasta'));
	
    eval {ok(defined($seq = $db->get_Seq_by_acc('J02231')))};
	skip('could not connect to embl',2) if $@;
    like( $seq->id, qr/J02231/);
    is( $seq->length, 200); 
    ok( defined($db = Bio::DB::EMBL->new(-verbose=>$verbose, 
					-retrievaltype => 'tempfile')));
    eval {ok(defined($seqio = $db->get_Stream_by_id(['AEE33958'])))};
	skip('could not connect to embl',2) if $@;
    undef $db; # testing to see if we can remove gb
    ok( defined($seq = $seqio->next_seq()));
    cmp_ok( $seq->length, '>=', 1);
}

$seq = $seqio = undef;

SKIP: {
    $db = Bio::DB::EMBL->new(-verbose => $verbose,
			    -retrievaltype => 'tempfile',
			    -format => 'fasta'
			    ); 
    eval{ok( defined($seqio = $db->get_Stream_by_acc(['J00522 AF303112 J02231'])))};
	skip('could not connect to embl',3) if $@;
    my %seqs;
    # don't assume anything about the order of the sequences
    while ( my $s = $seqio->next_seq ) {
		my ($type,$x,$name) = split(/\|/,$s->display_id);
		$seqs{$x} = $s->length;
    }
    is($seqs{'J00522'},408);
    is($seqs{'AF303112'},1611);
    is($seqs{'J02231'},200);
}
