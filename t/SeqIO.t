# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($DEBUG);
BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    plan tests => 35 }

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqIO::MultiFile;

ok(1);


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run.

my $verbosity = -1;   # Set to -1 for release version, so warnings aren't printed

my ($str, $seq,$ast,$temp,$mf,$ent,$out); # predeclare variables for strict
$str = Bio::SeqIO->new(-file=> 't/test.fasta', '-format' => 'Fasta');
ok $str;

ok ($seq = $str->next_seq());

print "Sequence 1 of 2 from fasta stream:\n", $seq->seq, "\n" if ( $DEBUG);

ok $seq->id, 'roa1_drome' ;

ok $seq->length, 358;

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

ok $str;

ok ($seq = $str->next_seq());
print "Sequence 1 of 2 from Raw stream:\n", $seq->seq, "\n\n" if( $DEBUG);

ok ($seq = $str->next_seq());
    
print "Sequence 2 of 2 from Raw stream:\n", $seq->seq, $seq->seq, "\n" 
    if( $DEBUG);


## Now we test Bio::SeqIO::GCG

$str = Bio::SeqIO->new(-file=> 't/test.gcg', '-format' => 'GCG');

ok $str;

ok ( $seq = $str->next_seq());
print "Sequence 1 of 1 from GCG stream:\n", $seq->seq, "\n" if( $DEBUG);

## Now we test Bio::SeqIO::GCG output writing

$str = Bio::SeqIO->new(-file=> '>t/gcg.out', '-format' => 'GCG');

$str->write_seq($seq);
ok(1);

#####
## End of ChrisDag's SeqIO tests.
#####

## Now we test Bio::SeqIO::GenBank
$str = Bio::SeqIO->new( -file=> 't/test.genbank', '-format' => 'GenBank');

ok $str;
$str->verbose($verbosity);

ok ( $seq = $str->next_seq() );
print "Sequence 1 of 1 from GenBank stream:\n", $seq->seq, "\n" if( $DEBUG);


my $strout = Bio::SeqIO->new(-file=> '>t/genbank.out', '-format' => 'GenBank');
while( $seq ) {
    $strout->write_seq($seq);
    $seq = $str->next_seq();
}
undef $strout;
ok(1);

# please leave this as the last line:
$str = undef;

# EMBL format

$ast = Bio::SeqIO->new( '-format' => 'embl' , -file => 't/roa1.dat');
$ast->verbose($verbosity);
my $as = $ast->next_seq();
ok defined $as->seq;


$ast = Bio::SeqIO->new( '-format' => 'GenBank' , 
			-file => 't/roa1.genbank');
$ast->verbose($verbosity);
$as = $ast->next_seq();
ok defined $as->seq;

$mf = Bio::SeqIO::MultiFile->new( '-format' => 'Fasta' , 
				  -files => ['t/multi_1.fa','t/multi_2.fa']);

ok defined $mf;

# read completely to the end
eval { 
    while( $seq = $mf->next_seq() ) {
	$temp = $seq->display_id;
    }
};
ok ! $@;
$temp = undef;
$ast = Bio::SeqIO->new( '-format' => 'Swiss' , -file => 't/roa1.swiss');
$ast->verbose($verbosity);
$as = $ast->next_seq();
ok defined $as->seq;
ok $as->id, 'ROA1_HUMAN', "id is ".$as->id;

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
ok defined $seq->seq(), 1, 'failure to read Embl with ^ location and badly split double quotes';

$out = Bio::SeqIO->new(-file=> '>t/embl.out', '-format' => 'embl');

# test writing the same
ok $out->write_seq($seq),1,'failure to write Embl format with ^ < and > locations';

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

    ok @a_seq, 3, 'wrong number of sequence objects';

    my $esc_name = $a_seq[1]->display_id;
    ok( $esc_name , 'Name; 4% strewn with \ various / escaped characters', 
	"bad unescaping of characters, $esc_name");
    
    ok $a_seq[0]->moltype, 'protein', 'moltypes incorrectly detected';
    ok $a_seq[1]->moltype, 'dna', 'moltypes incorrectly detected';
    
    # Test writing
    my $o_file = 't/test.out.ace';
    my $a_out = Bio::SeqIO->new( -FILE => "> $o_file", -FORMAT => 'ace');
    my $a_out_ok = 1;
    foreach my $a (@a_seq) {
        $a_out->write_seq($a) or $a_out_ok = 0;
    }
    undef($a_out);  # Flush to disk
    ok $a_out_ok,1,'error writing sequence';
    
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
    ok( ($before and $after and ($before eq $after)),1, 
	'test output file differs from input');
}

#
# Tests for feature-rich GenBank-entries. Added by HL <hlapp@gmx.net> 05/07/00
#
my $stream = Bio::SeqIO->new('-file' => 't/test.genbank',
			     '-format' => 'GenBank');
$stream->verbose($verbosity);
my $seqnum = 0;
my $species;
my @cl;
my $lasts;
while($seq = $stream->next_seq()) {
    $seqnum++;
    if($seqnum == 3) {
	ok $seq->display_id(), "HUMBDNF";
	# check for correct recognition of species
	$species = $seq->species();
	@cl = $species->classification();
	ok( $species->binomial(), "Homo sapiens", 
	    'species parsing incorrect for genbank');
	ok( $cl[3] ne $species->genus(), 1, 
	    'genus duplicated in genbank parsing');
    }
    # features which used to screw up the genbank/feature table parser
    $lasts = $seq;
}
ok $lasts->display_id(), "HUMBETGLOA";
$stream->close();
#
# we add a test regarding duplication of genus for EMBL as well.
#
$ent = Bio::SeqIO->new( -FILE => 't/test.embl', -FORMAT => 'embl');
$ent->verbose($verbosity);
$seq = $ent->next_seq();
$species = $seq->species();
@cl = $species->classification();
ok( $cl[3] ne $species->genus(), 1, 'genus duplicated in EMBL parsing');
$ent->close();


# let's test to see how well we handle FuzzyLocations
$seq = Bio::SeqIO->new( '-format' => 'GenBank' , 
			-file => 't/testfuzzy.genbank');
$seq->verbose($verbosity);
$as = $seq->next_seq();
ok defined $as->seq;

$seq = Bio::SeqIO->new( '-format' => 'GenBank' , 
			-file => '>t/genbank.fuzzyout');
$seq->verbose($verbosity);
ok($seq->write_seq($as));
