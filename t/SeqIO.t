# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..31\n";
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use lib '../';
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqIO::MultiFile;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run.


sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

$str = Bio::SeqIO->new(-file=> 't/test.fasta', '-format' => 'Fasta');
test 2, ( $str );

test 3, ($seq = $str->next_seq()), 'failed to read fasta seq from stream'; 

print "Sequence 1 of 2 from fasta stream:\n", $seq->seq, "\n";

test 4, ( $seq->id eq 'roa1_drome' );

test 5, ($seq->length == 358);

#####
## ChrisDag -- testing out Bio::SeqIO::Raw & SeqIO::GCG
##
## We open a file, test.raw which has 2 raw lines of
## sequence. No formatting at all. Raw sequences are delimited
## simply by a newline. This code tests to make sure we can
## create 2 sequential bioseq objects out of the raw file without
## breaking or getting confused.
##
## Not tested yet: ability to write a raw formatted stream
## Not tested yet: ability to write a GCG formatted stream

$str = Bio::SeqIO->new(-file=> 't/test.raw', '-format' => 'Raw');

test 6, ( $str ), 'unable to open stream from raw sequence DB';

test 7, ( $seq = $str->next_seq() ), 'failed to read 1st raw sequence from stream';
print "Sequence 1 of 2 from Raw stream:\n", $seq->seq, "\n\n";

test 8, ($seq = $str->next_seq()), 'failed to read 2nd raw sequence from stream'; 
print "Sequence 2 of 2 from Raw stream:\n", $seq->seq, $seq->seq, "\n";


## Now we test Bio::SeqIO::GCG

$str = Bio::SeqIO->new(-file=> 't/test.gcg', '-format' => 'GCG');

test 9, ( $str ), 'unable to open stream from GCG sequence file';

test 10,( $seq = $str->next_seq()),'failed to read GCG sequence from stream'; 
print "Sequence 1 of 1 from GCG stream:\n", $seq->seq, "\n";

## Now we test Bio::SeqIO::GCG output writing

$str = Bio::SeqIO->new(-file=> '>t/gcg.out', '-format' => 'GCG');

$str->write_seq($seq);

test 11, 1;

#####
## End of ChrisDag's SeqIO tests.
#####

## Now we test Bio::SeqIO::GenBank
$str = Bio::SeqIO->new( -file=> 't/test.genbank', '-format' => 'GenBank');

test 12, ( $str ), 'unable to open stream from GenBank sequence file';
$str->verbose(-1);    # Set to -1 for release version, so warnings aren't printed

test 13, ( $seq = $str->next_seq() ),'failed to read GenBank sequence from stream'; 
print "Sequence 1 of 1 from GenBank stream:\n", $seq->seq, "\n";


$str = Bio::SeqIO->new(-file=> '>t/genbank.out', '-format' => 'GenBank');

$str->write_seq($seq);
test 14, 1;

# please leave this as the last line:
$str = undef;

# EMBL format

$ast = Bio::SeqIO->new( '-format' => 'embl' , -file => 't/roa1.dat');

my $as = $ast->next_seq();
test 15, defined $as->seq;


$ast = Bio::SeqIO->new( '-format' => 'GenBank' , -file => 't/roa1.genbank');

$as = $ast->next_seq();
test 16, defined $as->seq;

$mf = Bio::SeqIO::MultiFile->new( '-format' => 'Fasta' , -files => ['t/multi_1.fa','t/multi_2.fa']);

test 17, 1;

# read completely to the end
eval { 
    while( $seq = $mf->next_seq() ) {
	$temp = $seq->display_id;
    }
};
test 18, ! $@;
$temp = undef;
$ast = Bio::SeqIO->new( '-format' => 'Swiss' , -file => 't/roa1.swiss');

$as = $ast->next_seq();
test 19,  defined $as->seq && $as->id eq 'ROA1_HUMAN', "id is ".$as->id;

# Keith James' tests for SeqIO reading EMBL features with:
#
#  * locations containing < and/or >
#  * zero width features denoted by ^
#  * qualifiers containing terminal " resulting from "" mistakenly
#    split across two lines (from AceDB, I gather)

($ent, $seq, $out) = undef;

$ent = Bio::SeqIO->new( -FILE => 't/test.embl', -FORMAT => 'embl');
$seq = $ent->next_seq();

# test reading file
test 20, defined $seq->seq(),'failure to read Embl with ^ location and badly split double quotes';

$out = Bio::SeqIO->new(-file=> '>t/embl.out', '-format' => 'embl');

# test writing the same
test 21, $out->write_seq($seq),'failure to write Embl format with ^ < and > locations';

# ACeDB flatfile (ace) sequence format tests
{
    my $t_file = 't/test.ace';
    my( $before );
    {
        local $/ = undef;
        local *BEFORE;
        open BEFORE, $t_file;
        $before = <BEFORE>;
        close BEFORE;
    }

    # Test reading
    my $a_in = Bio::SeqIO->new( -FILE => $t_file, -FORMAT => 'ace');
    my( @a_seq );
    while (my $a = $a_in->next_seq) {
        push(@a_seq, $a);
    }

    test 22, (@a_seq == 3), 'wrong number of sequence objects';

    my $esc_name = $a_seq[1]->display_id;
    test 23, ($esc_name eq 'Name; 4% strewn with \ various / escaped characters'), "bad unescaping of characters, $esc_name";
    
    test 24, ($a_seq[0]->moltype eq 'protein' and 
	      $a_seq[1]->moltype eq 'dna'), 'moltypes incorrectly detected';
    
    # Test writing
    my $o_file = 't/test.out.ace';
    my $a_out = Bio::SeqIO->new( -FILE => "> $o_file", -FORMAT => 'ace');
    my $a_out_ok = 1;
    foreach my $a (@a_seq) {
        $a_out->write_seq($a) or $a_out_ok = 0;
    }
    undef($a_out);  # Flush to disk
    test 25, ($a_out_ok),'error writing sequence';
    
    my( $after );
    {
        local $/ = undef;
        local *AFTER;
        open AFTER, $o_file;
        $after = <AFTER>;
        close AFTER;
    }
    unlink($o_file);
    
    # Test that input and output files are identical
    test 26, ($before and $after and ($before eq $after)), 
    'test output file differs from input';
}

#
# Tests for feature-rich GenBank-entries. Added by HL <hlapp@gmx.net> 05/07/00
#
my $stream = Bio::SeqIO->new('-file' => 't/test.genbank',
			     '-format' => 'GenBank');
$stream->verbose(-1);    # Set to -1 for release version, so warnings aren't printed
my $seqnum = 0;
my $species;
my @cl;
my $lasts;
while($seq = $stream->next_seq()) {
    $seqnum++;
    if($seqnum == 3) {
	test 27, ($seq->display_id() eq "HUMBDNF");
	# check for correct recognition of species
	$species = $seq->species();
	@cl = $species->classification();
	test 28, ($species->binomial() eq "Homo sapiens"), 'species parsing incorrect for genbank';
	test 29, ($cl[3] ne $species->genus()), 'genus duplicated in genbank parsing';
    }
    # features which used to screw up the genbank/feature table parser
    $lasts = $seq;
}
test 30, ($lasts->display_id() eq "HUMBETGLOA");
$stream->close();
#
# we add a test regarding duplication of genus for EMBL as well.
#
$ent = Bio::SeqIO->new( -FILE => 't/test.embl', -FORMAT => 'embl');
$ent->verbose(-1);
$seq = $ent->next_seq();
$species = $seq->species();
@cl = $species->classification();
test 31, ($cl[3] ne $species->genus()), 'genus duplicated in EMBL parsing';
$ent->close();
