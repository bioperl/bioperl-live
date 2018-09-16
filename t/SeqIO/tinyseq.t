# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 16,
               -requires_modules => [qw(XML::Parser::PerlSAX XML::Writer)]);

    use_ok('Bio::SeqIO::tinyseq');
}

my $file    = test_input_file('test.tseq');
my $outfile = test_output_file();

my $instream = Bio::SeqIO->new( -file    => $file,
                                -format  => 'tinyseq' );

my $outstream = Bio::SeqIO->new( -file   => ">$outfile",
                                 -format => 'tinyseq' );

my $seq = $instream->next_seq;
ok(defined $seq);
ok(defined $seq->seq);
is($seq->length, 5830);
is($seq->accession_number,'NM_002253');
ok($seq->species);
is($seq->species->binomial, 'Homo sapiens');
is($seq->species->ncbi_taxid, 9606);
$outstream->write_seq($seq);
undef $outstream;

ok(-s $outfile);

my $reread = Bio::SeqIO->new( -file   => $outfile,
                              -format => 'tinyseq' );

my $seq2 = $reread->next_seq;

ok($seq2);
ok($seq2->seq);
is($seq2->length, 5830);
is($seq2->accession_number, 'NM_002253');
ok($seq2->species);
is($seq2->species->binomial, 'Homo sapiens');
is($seq2->species->ncbi_taxid, 9606);
