#!/usr/bin/perl -w

#---------------------------------------------------------------------------
# PROGRAM  : seqs2.pl
# PURPOSE  : To demonstrate basic sequence manipulation using seqtools.pl,
#            specifically:
#               -- how to screen a set of sequences given a list
#                  of sequence IDs.
#               -- the use of the get_seq_objs() function.
#               -- write each sequence to a separate file.
# AUTHOR   : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED  : 10 Apr 1998
# REVISION : $Id$
# USAGE    : seqs2.pl -h
# EXAMPLES : seqs2.pl -eg
#
# INSTALLATION: 
#    Set the require "../seqtools.pl" to point to the proper location
#    of the seqtools.pl file
#
# There are many other possible things you can do like
# sorting the sequences, collecting stats, further editing,
# comparisons, to name but a few.
#
# MODIFIED:
#  sac, 16 Jun 1998: Added installation comment, require statement comments.
#
# SEE ALSO : seqs1.pl, seqs3.pl, seqs4.pl
#---------------------------------------------------------------------------

# Using seqtools.pl in the examples/blast distribution directory:
require "seqtools.pl"; 
# Proper path to seqtools.pl after you install it in your system:
#require "/share/www-data/html/perlOOP/bioperl/lib/Bio/drivers/seq/seqtools.pl";

use vars qw($ID $VERSION $DESC $opt_write_files);

$ID      = 'seqs2.pl';
$VERSION = 0.1;
$DESC    = "Demonstrates the use of &get_seq_objs() from seqtools.pl";

&init_seq();
&load_ids();

# &print_seq or &write_file will get called on each seq object 
# as it gets parsed from the input sequence file.
# &print_seq and &write_files are defined in seqtools.pl

$opt_write_files ? &get_seq_objs(\&write_file)
                 : &get_seq_objs(\&print_seq);   


&wrap_up_seq();


#------------
sub examples {
#------------
<<"QQ_EG_QQ";
(using the files in this directory)

  # Basic read and write.
  ./$ID *.fasta > out.fasta

  # Screen a set of sequences with a set of sequence IDs
  ./$ID -incl id.list -noexact seq1.fasta 
 
  # Convert a set of Fasta sequences to a set of GCG formatted files.
  ./$ID -eid seq1.fasta -outfmt gcg -write_files tmp/

  # Verify protein sequence data:
  ./$ID -prot -strict -nomon < seq1.fasta > out.fasta

QQ_EG_QQ
}





