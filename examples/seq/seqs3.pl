#!/usr/local/bin/perl -w

#---------------------------------------------------------------------------
# PROGRAM  : seqs3.pl
# PURPOSE  : To demonstrate basic sequence manipulation using seqtools.pl,
#            specifically, showing how to use the get_seq_data() function
#            to work with sequence data directly rather than sequence objects.
# AUTHOR   : Steve A. Chervitz (sac@genome.stanford.edu).
# CREATED  : 10 Apr 1998
# REVISION : $Id$
# USAGE    : seqs3.pl -h
# EXAMPLES : seqs3.pl -eg
#
# INSTALLATION: 
#    Set the require "../seqtools.pl" to point to the proper location
#    of the seqtools.pl file
#
# MODIFIED:
#  sac, 16 Jun 1998: Added installation comment, require statement comments.
#
# SEE ALSO : seqs1.pl, seqs2.pl, seqs4.pl
#---------------------------------------------------------------------------

# Using seqtools.pl in the examples/blast distribution directory:
require "seqtools.pl"; 
# Proper path to seqtools.pl after you install it in your system:
#require "/share/www-data/html/perlOOP/bioperl/lib/Bio/drivers/seq/seqtools.pl";

use vars qw($ID $VERSION $DESC);

$ID      = 'seqs3.pl';
$VERSION = 0.1;
$DESC    = "Demonstrates the use of &get_seq_data() from seqtools.pl";

&init_seq();
&load_ids();
&get_seq_data(\&process_data);
&wrap_up_seq();


#------------
sub examples {
#------------
<<"QQ_EG_QQ";
(using the files in this directory)

  # Basic read and write.
  ./$ID *.fasta > out.fasta

QQ_EG_QQ
}


#----------------
sub process_data {
#----------------
# Process a single sequence dataset.
# Data is provided separately rather than in Seq objects.
# It's a bit faster and would be useful if you only need the seq data and
# don't particularly need to use sequence objects.

    my ($id, $desc, $seq) = @_;

    printf "ID   = $id\n";
    printf "DESC = $desc\n";
    printf "SEQ  = \n$seq\n\n";
}    



