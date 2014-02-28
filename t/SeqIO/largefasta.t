# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 16);

    use_ok('Bio::SeqIO::largefasta');
}

my $tmpfile = test_output_file();

my $seqio = Bio::SeqIO->new('-format' => 'largefasta',
                            '-file'   => test_input_file('genomic-seq.fasta'),
                            );
isa_ok($seqio, 'Bio::SeqIO');

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

is open(OUT, '>', $tmpfile), 1;

my $seqout = Bio::SeqIO->new('-format' => 'largefasta',
                             '-fh'     => \*OUT );
is defined $seqout, 1;

is $seqout->write_seq($pseq), 1;
$seqout->close();
close(OUT);
my $seqin = Bio::SeqIO->new('-format' => 'largefasta',
                            '-file'   => $tmpfile);
my $pseq2 = $seqin->next_seq;
is ($plength, $pseq2->length());
is ($pseq->display_id(), $pseq2->display_id());
is ($pseq->desc(), $pseq2->desc());
