#!/usr/bin/perl -w

#---------------------------------------------------------------------------
# PROGRAM  : seqs2.pl
# PURPOSE  : To demonstrate basic sequence manipulation using seqtools.pl,
#            specifically:
#               -- how to screen a set of sequences given a list
#                  of sequence IDs.
#               -- write each sequence to a separate file.
# AUTHOR   : Steve Chervitz (sac@bioperl.org)
# CREATED  : 10 Apr 1998
# REVISION : $Id$
# USAGE    : seqs2.pl -h
# EXAMPLES : seqs2.pl -eg
#
# INSTALLATION: 
#    Set the require "seqtools.pl" to point to the proper location
#    of the seqtools.pl file
#
# There are many other possible things you can do like
# sorting the sequences, collecting stats, further editing,
# comparisons, to name but a few.
#
# MODIFIED:
#  sac, 22 Feb 2000: Converted to use Bio::SeqIO. Little change.
#  sac, 16 Jun 1998: Added installation comment, require statement comments.
#
# SEE ALSO : seqs1.pl, seqs3.pl
#---------------------------------------------------------------------------

# Using seqtools.pl in the examples/blast distribution directory:
require "seqtools.pl"; 
# Proper path to seqtools.pl after you install it in your system:
#require "/home/steve/perl/bioperl/examples/seqio/seqtools.pl";

use vars qw($VERSION $DESC $opt_write_files);

$VERSION = 0.1;
$DESC    = "Demonstrates sequence screening and writing using procedures in seqtools.pl";


&init_seq();
&load_ids();

# &print_seq or &write_file will get called on each seq object 
# as it gets parsed from the input sequence file.
# &print_seq and &write_files are defined in seqtools.pl

$opt_write_files ? &load_seqs(\&write_file)
                 : &load_seqs(\&print_seq);   


&wrap_up_seq();


#------------
sub examples {
#------------
<<"QQ_EG_QQ";
(using the files in this directory)

  # Screen a set of sequences with a set of sequence IDs to include
  $0 -incl id.list -noexact seq1.fasta 
 
  # Screen a set of sequences with a set of sequence IDs to exclude
  $0 -excl id.list -noexact seq1.fasta 
 
  # Convert a set of Fasta sequences to a set of GCG formatted files.
  $0 seq1.fasta -outfmt gcg -write_files tmp/

QQ_EG_QQ
}





