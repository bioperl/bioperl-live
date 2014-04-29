# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests               => 36,
               -requires_modules    => [qw(IO::String
                                           LWP::UserAgent)],
               -requires_networking => 1);

    use_ok('Bio::DB::BioFetch');
}

my $verbose = test_debug();

my $dbwarn = "Warning: Couldn't connect to EMBL with Bio::DB::BioFetch!\n";

my ($db,$db2,$seq,$seqio);

SKIP :{
    # get a single seq
    ok defined($db = Bio::DB::BioFetch->new(-verbose => $verbose));
    # get a RefSeq entry
    ok $db->db('refseq');
    eval {
        $seq = $db->get_Seq_by_acc('NM_006732'); # RefSeq VERSION
    };
    skip($dbwarn, 4) if $@;
    isa_ok($seq, 'Bio::SeqI');
    is($seq->accession_number,'NM_006732');
    is($seq->accession_number,'NM_006732');
    is( $seq->length, 3776);
}

SKIP: {
    # EMBL
    $db->db('embl');
    eval {
        $seq = $db->get_Seq_by_acc('J02231');
    };
    skip($dbwarn, 3) if $@;
    isa_ok($seq, 'Bio::SeqI');
    is($seq->id, 'J02231');
    is($seq->length, 200);
}

SKIP: {
    eval {
        $seqio = $db->get_Stream_by_id(['AEE33958']);
    };
    skip($dbwarn, 3) if $@;
    undef $db; # testing to see if we can remove gb
    $seq = $seqio->next_seq();
    isa_ok($seqio, 'Bio::SeqIO');
    isa_ok($seq, 'Bio::SeqI');
    cmp_ok( $seq->length, '>=', 1);
}

SKIP: {
    #swissprot
    ok $db2 = Bio::DB::BioFetch->new(-db => 'swissprot');
    eval {
        $seq = $db2->get_Seq_by_id('YNB3_YEAST');
    };
    skip($dbwarn, 5) if $@;
    isa_ok($seq, 'Bio::SeqI');
    is($seq->length, 125);
    is($seq->division, 'YEAST');
    $db2->request_format('fasta');
    eval {
        $seq = $db2->get_Seq_by_acc('P43780');
    };
    skip($dbwarn, 2) if $@;
    isa_ok($seq, 'Bio::SeqI');
    is($seq->length,103);
}

$seq = $seqio = undef;

SKIP: {
    ok $db = Bio::DB::BioFetch->new(-retrievaltype => 'tempfile',
                                    -format        => 'fasta',
                                    -verbose       => $verbose
                                    );
    $db->db('embl');
    eval {
        $seqio = $db->get_Stream_by_id('J00522 AF303112 J02231');
    };
    skip($dbwarn, 7) if $@;
    my %seqs;
    # don't assume anything about the order of the sequences
    while ( my $s = $seqio->next_seq ) {
        isa_ok($s, 'Bio::SeqI');
        my ($type,$x,$name) = split(/\|/,$s->display_id);
        $seqs{$x} = $s->length;
    }
    isa_ok($seqio, 'Bio::SeqIO');
    is($seqs{'J00522'},408);
    is($seqs{'AF303112'},1611);
    is($seqs{'J02231'},200);
}

SKIP: {
    ok $db = Bio::DB::BioFetch->new(-db      => 'embl',
                                    -verbose => $verbose ? $verbose : -1);

    # check contig warning (WebDBSeqI)
    eval {
        $seq = $db->get_Seq_by_acc('NT_006732');
    };
    like($@, qr{contigs are whole chromosome files}, 'contig warning');
    eval {
        $seq = $db->get_Seq_by_acc('NM_006732');
    };
    skip($dbwarn, 2) if $@;
    isa_ok($seq, 'Bio::SeqI');
    is($seq->length,3776);
}

# unisave
SKIP: {
    ok $db = Bio::DB::BioFetch->new(-db      => 'unisave',
                                    -verbose => $verbose);
    eval {
        $seq = $db->get_Seq_by_acc('LAM1_MOUSE');
    };
    skip($dbwarn, 4) if $@;
    isa_ok($seq, 'Bio::SeqI');
    is($seq->display_id, 'LAM1_MOUSE');
    is($seq->accession, 'P14733');
    is($seq->length, 587);
}
