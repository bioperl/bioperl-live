#!/usr/bin/perl

# PROGRAM  : rev_and_trans.pl
# PURPOSE  : Simple driver for Bio::Seq revcom and translate
# AUTHOR   : Ewan Birney birney@sanger.ac.uk 
# CREATED  : Tue Oct 27 1998
#
# INSTALLATION
#    If you have installed bioperl using the standard
#    makefile system everything should be fine and 
#    dandy.
#
#    if not edit the use lib "...." line to point the directory
#    containing your Bioperl modules.
#


use Bio::Seq;
use Bio::SeqIO;

# new sequence from raw memory...
# it is *very* important to get the type right so it
# is translated correctly.

$seq = Bio::Seq->new ( -id => "myseq",
		      -seq => "CGCCGAAGAAGCATCGTTAAAGTCTCTCTTCACCCTGCCGTCATGTCTAAGTCAGAGTCTCCT",
		      -type => 'Dna');

$seqout = Bio::SeqIO->new('-format' => 'fasta', -fh => \*STDOUT);

# make a reverse complement sequence

$rev = $seq->revcom();

# the actual sequence is here

$actual_bases = $rev->seq();

print "Reversed sequence as a string is [$actual_bases]\n";

# we could also write it as fasta formatted output

$seqout->write_seq($rev);

# make a translation

$trans = $seq->translate();

print "Translated sequence!\n";

$seqout->write_seq($trans);

