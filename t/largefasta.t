# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;    
    plan tests => 15;
}

use Bio::SeqIO;
use vars qw($tmpfile);
use Bio::Root::IO;
END { unlink $tmpfile; }

$tmpfile = 't/largefastatest.out';
my $seqio = new Bio::SeqIO('-format'=>'largefasta',
			   '-file'  =>Bio::Root::IO->catfile("t","genomic-seq.fasta"));
ok defined $seqio, 1, 'cannot instantiate Bio::SeqIO::largefasta';

my $pseq = $seqio->next_seq();
$pseq->moltype('dna');
$pseq->desc('this is my description');;
my $plength = $pseq->length();
my $last_3 = $pseq->subseq($plength-3,$plength);

ok defined $pseq, 1, 'could not call next_seq';
ok $plength > 0, 1, "could not call length, seq was empty";
ok length($pseq->subseq(100, 299)), 200, 'error in subseq'; 
ok $pseq->trunc(100,199)->length(), 100, 'error in trunc'; 
ok $pseq->moltype(), 'dna', 'moltype was ' . $pseq->moltype();
ok $pseq->display_id(), 'HSBA536C5',"no display id";
ok $pseq->accession_number(), 'unknown', "no accession";
ok $pseq->desc, 'this is my description', 'no description';

ok open(OUT, ">$tmpfile"), 1,'could not open output file';

my $seqout = new Bio::SeqIO('-format' => 'largefasta',
			    '-fh'     => \*OUT );
ok defined $seqout, 1,'could not open seq with outputstream';

ok $seqout->write_seq($pseq), 1,'could not write seq';
$seqout->close();
close(OUT);
my $seqin = new Bio::SeqIO('-format' => 'largefasta',
			'-file'   => $tmpfile);
my $pseq2 = $seqin->next_seq;
ok ($plength, $pseq2->length(), 
    "written file was not same length as expected");
ok ($pseq->display_id(), $pseq2->display_id(), 
    "display ids were not identical as expected");
ok ($pseq->desc(), $pseq2->desc() , 
    "description was not identical (" . $pseq->desc() . 
    "," . $pseq2->desc() . ")");
