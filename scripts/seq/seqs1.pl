#!/usr/bin/perl -w

#----------------------------------------------------------------
# seqs1.pl
#
# A minimal script that allows basic printing, reformatting, and
# filtering of sequence data. Illustrates use of seqtools.pl 
# (which uses Bio::SeqIO).
# See seqs2.pl, and seqs3.pl for some more advanced scripts.
#
# Author  : Steve Chervitz <sac@bioperl.org>
#
# Revision: $Id$
#
# Usage   : seqs1.pl -h
# Examples: seqs1.pl -eg
#
# Installation:
#    Change require "seqtools.pl" to point to the proper location
#    of the seqtools.pl file within your system.
#
# Modified: 
#  sac, 21 feb 2000: Updated for use with Bio::SeqIO. Little change.
#  sac, 16 Jun 1998: Added installation comment, require statement comments.
#
#----------------------------------------------------------------

# Using seqtools.pl in the examples/seq distribution directory.
# If you install this script on your system, modify this require
# statement to have the absolute path to seqtools.pl
require "seqtools.pl";


use vars qw($VERSION $DESC);
$VERSION = 0.1;
$DESC = "This is a minimal script that uses the Bio::SeqIO module\n". 
 "via seqtools.pl. Allows basic printing & reformatting.";

&init_seq();
&load_ids();
&load_seqs();  
&print_seqs();
&wrap_up_seq();


#-----------
sub examples {
#-----------
<<"QQ_EG_QQ";
(using the files in this directory)

 $0 seq1.fasta
 $0 *.fasta
 gzip -cd seq.fasta.gz | $0 
 $0 -prot seq1.fasta -outfmt genbank -out seq1.gb
 $0 < seq2.fasta -outfmt swiss > out.swiss

QQ_EG_QQ
}

