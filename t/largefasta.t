# This is -*-Perl-*- code
# $Id$

use strict;
BEGIN { 
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;    
    plan tests => 17;
	use_ok('Bio::Root::IO');
	use_ok('Bio::SeqIO');
}

use vars qw($tmpfile);
END { unlink $tmpfile; }

$tmpfile = Bio::Root::IO->catfile("t","data","largefastatest.out");
my $seqio = new Bio::SeqIO('-format'=>'largefasta',
			   '-file'  =>Bio::Root::IO->catfile("t","data","genomic-seq.fasta"));
is defined $seqio, 1, 'Instantiate Bio::SeqIO::largefasta';

my $pseq = $seqio->next_seq();
$pseq->alphabet('dna');
$pseq->desc('this is my description');;
my $plength = $pseq->length();
my $last_3 = $pseq->subseq($plength-3,$plength);

is defined $pseq, 1;
is $plength > 0, 1;
is length($pseq->subseq(100, 299)), 200; 
is $pseq->trunc(100,199)->length(), 100; 
is $pseq->alphabet(), 'dna';
is $pseq->display_id(), 'HSBA536C5';
is $pseq->accession_number(), 'unknown';
is $pseq->desc, 'this is my description';

is open(OUT, ">$tmpfile"), 1;

my $seqout = new Bio::SeqIO('-format' => 'largefasta',
			    '-fh'     => \*OUT );
is defined $seqout, 1;

is $seqout->write_seq($pseq), 1;
$seqout->close();
close(OUT);
my $seqin = new Bio::SeqIO('-format' => 'largefasta',
			'-file'   => $tmpfile);
my $pseq2 = $seqin->next_seq;
is ($plength, $pseq2->length());
is ($pseq->display_id(), $pseq2->display_id());
is ($pseq->desc(), $pseq2->desc());
