#!/usr/bin/perl -w

# A minimal script using seqtools.pl
# See seqs2.pl, seqs3.pl, and seqs4.pl for some more advanced scripts.
# Author  : Steve A. Chervitz (sac@genome.stanford.edu)
# Revision: $Id$
# Usage   : seqs1.pl -h
# Modified: 
#  sac, 16 Jun 1998: Added installation comment, require statement comments.

# Using seqtools.pl in the examples/blast distribution directory:
require "seqtools.pl"; 
# Proper path to seqtools.pl after you install it in your system:
#require "/share/www-data/html/perlOOP/bioperl/lib/Bio/drivers/seq/seqtools.pl";

use vars qw($ID $VERSION);
$ID      = 'seqs1.pl';
$VERSION = 0.1;

&init_seq(\&_usage);
&load_ids();
&get_seq_objs();  # not passing any function ref means all seq objects will be saved.
&print_seqs();
&wrap_up_seq();


#-----------
sub _usage {
#-----------
   print STDERR <<"QQ_USAGE_QQ";

Usage: $ID seq1.fasta
       $ID -prot seq1.fasta >& err.fasta
       $ID seq.fasta.gz -eid 
       $ID *.fasta
       $ID < seq2.fasta > out.fasta

 This is a minimal script that uses the Bio::Tools::Fasta.pm module 
 via seqtools.pl. Doesn't do much beyond loading a set of sequences 
 from a Fasta sequence file.

QQ_USAGE_QQ

}
