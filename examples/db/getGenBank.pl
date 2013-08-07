#!/usr/local/bin/perl
#
# How to retrieve GenBank entries over the Web
#
# by Jason Stajich
#
use Bio::DB::GenBank;
use Bio::SeqIO;
my $gb = new Bio::DB::GenBank;

# the output stream for your seqs, this can be a file
# instead or STDOUT, see the Bio::SeqIO module for info

my $seqout = new Bio::SeqIO(-fh => \*STDOUT, -format => 'fasta');

# if you want a single seq
my $seq = $gb->get_Seq_by_id('J00522');
$seqout->write_seq($seq);
# or by accession
$seq = $gb->get_Seq_by_acc('AF303112');

$seqout->write_seq($seq);

# feel free to pull multiple sequences...
# if you want to get a bunch of sequences use the get_Stream_by_id/acc methods
my $seqio = $gb->get_Stream_by_id([ qw(J00522 AF303112 2981014)]);

while( defined ($seq = $seqio->next_seq )) {
        $seqout->write_seq($seq);
}
