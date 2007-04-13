# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $DEBUG);

use lib '.','./blib/lib';

my $error;

BEGIN { 
    eval { require Test::More; };
	$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;

    $NUMTESTS = 36;
	
	eval { require IO::String; require LWP::UserAgent; };
	if ($@) {
		plan skip_all => "No LWP::UserAgent or IO::String installed; skipping tests";
	} elsif (!$DEBUG) {
		plan skip_all => 'Must set BIOPERLDEBUG=1 for network tests';
	} else {
		plan tests => $NUMTESTS;
	}
}

require_ok('Bio::DB::BioFetch');

my $verbose = -1;

## End of black magic.

my $dbwarn = "Warning: Couldn't connect to EMBL with Bio::DB::EMBL.pm!\n";

my ($db,$db2,$seq,$seqio);

SKIP :{
# get a single seq
ok defined($db = new Bio::DB::BioFetch(-verbose => $verbose));
# get a RefSeq entry
ok $db->db('refseq');
eval {
	$seq = $db->get_Seq_by_acc('NM_006732'); # RefSeq VERSION
};
skip($dbwarn, 4) if $@;
isa_ok($seq, 'Bio::SeqI');
is($seq->accession_number,'NM_006732');
is($seq->accession_number,'NM_006732');
is( $seq->length, 3775);
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
	$seqio = $db->get_Stream_by_id(['BUM']);
};
skip($dbwarn, 3) if $@;
undef $db; # testing to see if we can remove gb
$seq = $seqio->next_seq();
isa_ok($seqio, 'Bio::SeqIO');
isa_ok($seq, 'Bio::SeqI');
is( $seq->length, 200);
}

SKIP: {
#swissprot
ok $db2 = new Bio::DB::BioFetch(-db => 'swissprot');
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
ok $db = new Bio::DB::BioFetch(-retrievaltype => 'tempfile',
				 -format => 'fasta',
				 -verbose => $verbose
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
$verbose = -1;
ok $db = new Bio::DB::BioFetch(-db => 'embl',
			   -verbose => $verbose);

# check contig warning (WebDBSeqI)
eval {
	$seq = $db->get_Seq_by_acc('NT_006732');
};
like($@, qr{contigs are whole chromosome files}, 'contig warning');
eval {
	$seq = $db->get_Seq_by_acc('NM_006732');
};
skip($dbwarn, 3) if $@;
isa_ok($seq, 'Bio::SeqI');
is($seq->length, 3775);   
}
    
# unisave
SKIP: {
ok $db = new Bio::DB::BioFetch(-db => 'unisave',
			   -verbose => $verbose);
eval {
	$seq = $db->get_Seq_by_acc('LAM1_MOUSE');
};
skip($dbwarn, 3) if $@;
isa_ok($seq, 'Bio::SeqI');
is($seq->display_id, 'LAM1_MOUSE');
is($seq->accession, 'P14733');   
is($seq->length, 587);   
}
