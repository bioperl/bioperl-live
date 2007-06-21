# This is -*-Perl-*- code
# $Id$

use strict;

BEGIN {
	use lib 't/lib';
	use BioperlTest;
	
	test_begin(-tests => 116,
			   -requires_modules => [qw(IO::String
									    LWP::UserAgent
										HTTP::Request::Common)],
			   -requires_networking => 1);
	
	use_ok('Bio::DB::GenBank');
	use_ok('Bio::DB::GenPept');
	use_ok('Bio::DB::SwissProt');
    use_ok('Bio::DB::GDB');
    use_ok('Bio::DB::MeSH');
}

my %expected_lengths = ('NDP_MOUSE' => 131,
                        'NDP_HUMAN' => 133,
                        'MUSIGHBA1' => 408,
                        'AF303112'  => 1611,
                        'J00522'    => 408,
                        'AF303112'  => 1611,
                        'AF303112.1' => 1611,
                        '2981014'   => 1156,
                        'AF041456'  => 1156,
                        'AY080910'  => 798,
                        'AY080909'  => 1042,
                        'AF155220'  => 1172,
                        '405830'    => 1743,
                        'CELRABGDI' => 1743,
                        '195055'    => 136,
                        'AAD15290'  => 136,
                        'AAC06201'  => 353,
                        'P43780'    => 103,
                        'BOLA_HAEIN'=> 103,
                        'YNB3_YEAST'=> 125,
                        'O39869'    => 56,
                        'P18584'    => 497,
                        'DEGP_CHLTR'=> 497,
                        'AF442768'  => 2547,
                        'P31383'    => 635,
                        'CH402638'  => 5041);

my ($gb, $seq, $seqio, $seqin, $query);

#
# Bio::DB::GenBank
#
ok $gb = Bio::DB::GenBank->new('-delay'=>0), 'Bio::DB::GenBank';

# get a single seq
SKIP: {
    eval {$seq = $gb->get_Seq_by_id('MUSIGHBA1');};
    skip "Couldn't connect to Genbank with Bio::DB::GenBank.pm. Do you have network access? Skipping GenBank tests", 4 if $@;
    is $seq->length, $expected_lengths{$seq->display_id};
    eval {$seq = $gb->get_Seq_by_acc('AF303112');};
    skip "Couldn't connect to Genbank with Bio::DB::GenBank.pm. Transient network problems? Skipping GenBank tests", 3 if $@;
    is $seq->length, $expected_lengths{$seq->display_id};
    eval {$seq = $gb->get_Seq_by_version('AF303112.1');};
    skip "Couldn't connect to Genbank with Bio::DB::GenBank.pm. Transient network problems? Skipping GenBank tests", 2 if $@;
    is $seq->length, $expected_lengths{$seq->display_id};
    eval {$seq = $gb->get_Seq_by_gi('405830');};
    skip "Couldn't connect to Genbank with Bio::DB::GenBank.pm. Transient network problems? Skipping GenBank tests", 1 if $@;
    is $seq->length, $expected_lengths{$seq->display_id};
}

$seq = $seqio = undef;

# batch mode
SKIP: {
    eval {$seqio = $gb->get_Stream_by_id([qw(J00522 AF303112 2981014)]);};
    skip "Batch access test failed for Genbank. Skipping those tests", 4 if $@;
    my $done = 0;
    while (my $s = $seqio->next_seq) {
        is $s->length, $expected_lengths{$s->display_id};
        $done++;
    }
    skip('No seqs returned', 4) if !$done;
    is $done, 3;
}

$seq = $seqio = undef;

# test the temporary file creation and fasta
ok $gb = Bio::DB::GenBank->new('-format' => 'fasta', '-retrievaltype' => 'tempfile', '-delay' => 0);
SKIP: {
    eval {$seq = $gb->get_Seq_by_id('MUSIGHBA1');};
    skip "Couldn't connect to complete GenBank tests with a tempfile with Bio::DB::GenBank.pm. Skipping those tests", 6 if $@;
    # last part of id holds the key
    is $seq->length, $expected_lengths{(split(/\|/,$seq->display_id))[-1]};
    eval {$seq = $gb->get_Seq_by_acc('AF303112');};
    skip "Couldn't connect to complete GenBank tests with a tempfile with Bio::DB::GenBank.pm. Skipping those tests", 5 if $@;
    # last part of id holds the key
    is $seq->length, $expected_lengths{(split(/\|/,$seq->display_id))[-1]};
    # batch mode requires genbank format
    $gb->request_format("gb");
    eval {$seqio = $gb->get_Stream_by_id([qw(J00522 AF303112 2981014)]);};
    skip "Couldn't connect to complete GenBank batch tests with a tempfile with Bio::DB::GenBank.pm. Skipping those tests", 4 if $@;
    my $done = 0;
    while (my $s = $seqio->next_seq) {
        is $s->length, $expected_lengths{$s->display_id};
        undef $gb; # test the case where the db is gone, 
        # but a temp file should remain until seqio goes away.
        $done++;
    }
    skip('No seqs returned', 4) if !$done;
    is $done, 3;
}

$seq = $seqio = undef;

# test pipeline creation
ok $gb = Bio::DB::GenBank->new('-retrievaltype' => 'pipeline', '-delay' => 0);
SKIP: {
    eval {$seq = $gb->get_Seq_by_id('MUSIGHBA1');};
    skip "Couldn't connect to complete GenBank tests with a pipeline with Bio::DB::GenBank.pm. Skipping those tests", 6 if $@;
    is $seq->length, $expected_lengths{$seq->display_id};
    eval {$seq = $gb->get_Seq_by_acc('AF303112');};
    skip "Couldn't connect to complete GenBank tests with a pipeline with Bio::DB::GenBank.pm. Skipping those tests", 5 if $@;
    is $seq->length, $expected_lengths{$seq->display_id};
    eval {$seqio = $gb->get_Stream_by_id([qw(J00522 AF303112 2981014)]);};
    skip "Couldn't connect to complete GenBank tests with a pipeline with Bio::DB::GenBank.pm. Skipping those tests", 4 if $@;
    my $done = 0;
    while (my $s = $seqio->next_seq) {
        is $s->length, $expected_lengths{$s->display_id};
        undef $gb; # test the case where the db is gone, 
        # but the pipeline should remain until seqio goes away
        $done++;
    }
    skip('No seqs returned', 4) if !$done;
    is $done, 3;
}

$seq = $seqio = undef;

# test query facility
ok $query = Bio::DB::Query::GenBank->new('-db'      => 'nucleotide',
                                         '-query'   => 'Onchocerca volvulus[Organism]',
                                         '-mindate' => '2002/1/1',
                                         '-maxdate' => '2002/12/31'), 'Bio::DB::Query::GenBank';
SKIP: {
    cmp_ok $query->count, '>', 0;
    my @ids = $query->ids;
    cmp_ok @ids, '>', 0;
    is @ids, $query->count;
    ok $gb = Bio::DB::GenBank->new('-delay' => 0);
    eval {$seqio = $gb->get_Stream_by_query($query);};
    skip "Couldn't connect to complete GenBank query tests. Skipping those tests", 5 if $@;
    my $done = 0;
    while (my $s = $seqio->next_seq) {
        is $s->length, $expected_lengths{$s->display_id};
        undef $gb; # test the case where the db is gone, 
        # but the pipeline should remain until seqio goes away
        $done++;
    }
    skip('No seqs returned', 5) if !$done;
    is $done, 4;
}

$seq = $seqio = undef;

# test query facility (again)
ok $query = Bio::DB::Query::GenBank->new('-db'  => 'nucleotide',
                                         '-ids' => [qw(J00522 AF303112 2981014)]);
SKIP: {
    cmp_ok $query->count, '>', 0;
    my @ids = $query->ids;
    cmp_ok @ids, '>', 0;
    is @ids, $query->count;
    $gb = Bio::DB::GenBank->new('-delay' => 0);
    eval {$seqio = $gb->get_Stream_by_query($query);};
    skip "Couldn't connect to complete GenBank query tests. Skipping those tests: $@", 4 if $@;
    my $done = 0;
    while (my $s = $seqio->next_seq) {
        is $s->length, $expected_lengths{$s->display_id};
        $done++;
    }
    skip('No seqs returned', 4) if !$done;
    is $done, 3;
    $seqio->close(); # the key to preventing errors during make test, no idea why
}

$seq = $seqio = undef;

# and yet again, for bug 2133
$query = Bio::DB::Query::GenBank->new('-query'  => 'AF303112',
                                      '-ids' => [qw(J00522 AF303112 2981014)]);
is $query->query, 'J00522[PACC]|AF303112[PACC]|2981014[UID]';

# test contig retrieval
ok $gb = Bio::DB::GenBank->new('-delay'  => 0, '-format' => 'gbwithparts');
SKIP: {
    eval {$seq = $gb->get_Seq_by_id('CH402638');};
    skip "Couldn't connect to GenBank with Bio::DB::GenBank.pm. Skipping those tests", 3 if $@;
    is $seq->length, $expected_lengths{$seq->display_id};
    # now to check that postprocess_data in NCBIHelper catches CONTIG...
    ok $gb = Bio::DB::GenBank->new('-delay' => 0, '-format' => 'gb');
    eval {$seq = $gb->get_Seq_by_id('CH402638');};
    skip "Couldn't connect to GenBank with Bio::DB::GenBank.pm. Skipping those tests", 1 if $@;
    is $seq->length, $expected_lengths{$seq->display_id};
}

$seq = $seqio = undef;

# bug 1405
my @result;
ok $gb = Bio::DB::GenBank->new(-format => 'Fasta', -seq_start  => 2, -seq_stop   => 7);
SKIP: {
    eval {$seq = $gb->get_Seq_by_acc("A11111");};
    skip "Couldn't connect to complete GenBank tests. Skipping those tests", 15 if $@;
    is $seq->length, 6;
    # complexity tests
    ok $gb = Bio::DB::GenBank->new(-format => 'Fasta', -complexity => 0);
    eval {$seqin = $gb->get_Stream_by_acc("5");};
    skip "Couldn't connect to complete GenBank tests. Skipping those tests", 13 if $@;
    @result = (1136, 'dna', 342, 'protein');
    while ($seq = $seqin->next_seq) {
        is $seq->length, shift(@result);
        is $seq->alphabet, shift(@result);
    }
    is @result, 0;
    # Real batch retrieval using epost/efetch 
    # these tests may change if integrated further into Bio::DB::Gen*
    # Currently only useful for retrieving GI's via get_seq_stream
    $gb = Bio::DB::GenBank->new();
    eval {$seqin = $gb->get_seq_stream(-uids => [4887706 ,431229, 147460], -mode => 'batch');};
    skip "Couldn't connect to complete GenBank batchmode epost/efetch tests. Skipping those tests", 8 if $@;
    my %result = ('M59757' => 12611 ,'X76083'=> 3140, 'J01670'=> 1593);
	my $ct = 0;
    while ($seq = $seqin->next_seq) {
		$ct++;
		my $acc = $seq->accession;
        ok exists $result{ $acc };
        is $seq->length, $result{ $acc };
		delete $result{$acc};
    }
    skip('No seqs returned', 8) if !$ct;
	is $ct, 3;
    is %result, 0;
}

$seq = $seqin = undef;

#
# Bio::DB::GenPept
#
ok $gb = Bio::DB::GenPept->new();
SKIP: {
    eval {$seqin = $gb->get_seq_stream(-uids => [2981015, 1621261, 195055], -mode => 'batch');};
    skip "Couldn't connect to complete GenPept tests. Skipping those tests", 8 if $@;
    my %result = ('AAC06201' => 353, 'CAB02640' => 193, 'AAD15290' => 136);
    my $ct = 0;
    while ($seq = $seqin->next_seq) {
		$ct++;
		my $acc = $seq->accession;
        ok exists $result{ $acc };
        is $seq->length, $result{ $acc };
		delete $result{$acc};
    }
    skip('No seqs returned', 8) if !$ct;
	is $ct, 3;
    is %result, 0;
}

$seq = $seqio = undef;

ok $gb = Bio::DB::GenPept->new('-delay' => 0);
SKIP: { 
    eval {$seq = $gb->get_Seq_by_id('195055');};
    skip "Couldn't connect to Genbank with Bio::DB::GenPept.pm. Skipping those tests", 10 if $@;
    is $seq->length, $expected_lengths{$seq->display_id};
    eval {$seq = $gb->get_Seq_by_acc('AAC06201');};
    skip "Couldn't connect to Genbank with Bio::DB::GenPept.pm. Skipping those tests", 9 if $@;
    is $seq->length, $expected_lengths{$seq->display_id};
    eval {$seqio = $gb->get_Stream_by_id([qw(AAC06201 195055)]);};
    skip "Couldn't connect to Genbank with Bio::DB::GenPept.pm. Skipping those tests", 8 if $@;
    my $done = 0;
    while( my $s = $seqio->next_seq ) {
        is $s->length, $expected_lengths{$s->display_id};
        $done++;
    }
    skip('No seqs returned', 8) if !$done;
    is $done, 2;
    # swissprot genpept parsing   
    eval {$seq = $gb->get_Seq_by_acc('2AAA_YEAST');};
    skip "Couldn't connect to Genbank with Bio::DB::GenPept.pm. Skipping those tests", 5 if $@;
    is $seq->length, $expected_lengths{$seq->display_id};
    
    # test dbsource stuff
    # small chance this might change but hopefully not
    my @annot = $seq->annotation->get_Annotations('dblink');
    cmp_ok(scalar(@annot), '>', 31);	
    is $annot[0]->database, 'swissprot';
    is $annot[0]->primary_id, '2AAA_YEAST';
    is (($seq->annotation->get_Annotations('swissprot_dates'))[0]->value, 'Jul 1, 1993');
}

$seq = $seqio = undef;

#
# Bio::DB::SwissProt
#
ok $gb = Bio::DB::SwissProt->new(-retrievaltype =>'pipeline', -delay => 0);
SKIP: {
    eval {$seq = $gb->get_Seq_by_id('YNB3_YEAST');};
    skip "Couldn't connect to SwissProt with Bio::DB::Swiss.pm. Skipping those tests", 14 if $@;
    is $seq->length, $expected_lengths{$seq->display_id};
    is $seq->division, 'YEAST';
    
    eval {$seq = $gb->get_Seq_by_acc('P43780');};
    skip "Couldn't connect to SwissProt with Bio::DB::Swiss.pm. Skipping those tests", 12 if $@;
    is $seq->length, $expected_lengths{$seq->display_id};
    eval {$seq = $gb->get_Seq_by_acc('O39869');};
    skip "Couldn't connect to SwissProt with Bio::DB::Swiss.pm. Skipping those tests", 11 if $@;
    is $seq->length, $expected_lengths{$seq->accession_number};
    is $seq->accession_number, 'O39869';
    is $seq->division, '9PICO';
    
    # test for bug #958
    eval {$seq = $gb->get_Seq_by_id('P18584');};
    skip "Couldn't connect to SwissProt with Bio::DB::Swiss.pm. Skipping those tests", 8 if $@;
    is $seq->length, $expected_lengths{$seq->display_id};
    is $seq->display_id, 'DEGP_CHLTR';
    is $seq->division, 'CHLTR';

    ok $gb = Bio::DB::SwissProt->new('-retrievaltype' => 'tempfile', '-delay' => 0);
    eval {$seqio = $gb->get_Stream_by_id(['NDP_MOUSE', 'NDP_HUMAN']);};
    skip "Couldn't connect to SwissProt with Bio::DB::Swiss.pm. Skipping those tests", 4 if $@;
    undef $gb; # testing to see if we can remove gb
    ok $seq = $seqio->next_seq();
    is $seq->length, $expected_lengths{$seq->display_id};
    ok $seq = $seqio->next_seq();
    is $seq->length, $expected_lengths{$seq->display_id};
}

$seq = $seqio = undef;

#
# Bio::DB::GDB
#
ok my $gdb = Bio::DB::GDB->new();
SKIP: {
    my $info; 
    eval {$info = $gdb->get_info(-type => 'marker', -id => 'D1S243');};
    skip "Couldn't connect to GDB with Bio::DB::GDB.pm. Skipping those tests", 1 if $@;
    is $info->{gdbid}, 'GDB:188393';
}

#
# Bio::DB::EntrezGene
#
SKIP: {
	test_skip(-tests => 8, -requires_modules => ['Bio::ASN1::EntrezGene']);
    use_ok('Bio::DB::EntrezGene');
    ok $gb = Bio::DB::EntrezGene->new(-retrievaltype => 'tempfile', -delay => 0);
    eval {$seqio = $gb->get_Stream_by_id([2,3064]);};
    skip "Couldn't connect to Entrez with Bio::DB::EntrezGene. Skipping those tests", 6 if $@;
    $seq = $seqio->next_seq;
    is $seq->display_id, "A2M";
    is $seq->accession_number, 2;
    $seq = $seqio->next_seq;
    is $seq->display_id, "HD";
    is $seq->accession_number, 3064;
    eval {$seq = $gb->get_Seq_by_id(6099);};
    skip "Couldn't connect to Entrez with Bio::DB::EntrezGene. Skipping those tests", 2 if $@;
    is $seq->display_id, "RP";
    is $seq->accession_number, 6099;
}

$seq = $seqio = undef;

#
# Bio::DB::MeSH
#
ok my $mesh = Bio::DB::MeSH->new();
SKIP: {
    my $t;
    eval {$t = $mesh->get_exact_term('Dietary Fats');};
    skip "Couldn't connect to MeSH with Bio::DB::MeSH. Skipping those tests", 3 if $@;
    is $t->each_twig(), 2;
    eval {$t = $mesh->get_exact_term("Sinus Thrombosis, Intracranial");};
    skip "Couldn't connect to MeSH with Bio::DB::MeSH. Skipping those tests", 2 if $@;
    is $t->description, "Thrombus formation in an intracranial venous sinus, including the superior sagittal, cavernous, lateral, and petrous sinuses. Etiologies include thrombosis due to infection,  DEHYDRATION, coagulation disorders (see  THROMBOPHILIA), and  CRANIOCEREBRAL TRAUMA.";
    is $t->id, "D012851";
}

$seq = $seqio = undef;
