#!/usr/bin/perl

# PROGRAM  : rev_and_trans.pl
# PURPOSE  : Simple driver for Bio::Seq revcom and translate
# AUTHOR   : Ewan Birney birney@sanger.ac.uk 
# CREATED  : Tue Oct 27 1998
# REVISION : $Id$
#
# INSTALLATION
#    If you have installed bioperl using the standard
#    makefile system everything should be fine and 
#    dandy.
#
#    if not edit the use lib "...." line to point the directory
#    containing your Bioperl modules.
#

use lib "/nfs/disk21/birney/prog/bioperl/Bio";
use Bio::Seq;

# new sequence from raw memory...
# it is *very* important to get the type right so it
# is translated correctly.

$seq = Bio::Seq->new ( -id => "myseq",
		      -seq => "CGCCGAAGAAGCATCGTTAAAGTCTCTCTTCACCCTGCCGTCATGTCTAAGTCAGAGTCTCCT",
		      -type => 'Dna');

# make a reverse complement sequence

$rev = $seq->revcom();

# the actual sequence is here

$actual_bases = $rev->str();

print "Reversed sequence as a string is [$actual_bases]\n";

# we could also write it as fasta formatted output

print $rev->out_fasta();

# make a translation

$trans = $seq->translate();

print "Translated sequence!\n";

print $trans->out_fasta();

