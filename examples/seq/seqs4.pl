#!/usr/local/bin/perl -w

#---------------------------------------------------------------------------
# PROGRAM  : seqs4.pl
# PURPOSE  : To demonstrate basic sequence editing using Fasta.pm.
#            A set of sequences is input from a Fasta file, their descriptions
#            are edited, certain sequences are filtered out, and the new Fasta
#            sequences are output.
#            The basic point is to show how to exercise Bio::Tools::Fasta.pm.
# AUTHOR   : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED  : 29 Apr 1998
# REVISION : $Id$
# USAGE    : seqs4.pl -h
# EXAMPLES : seqs4.pl -eg
#
# INSTALLATION: 
#    Set the require "../seqtools.pl" to point to the proper location
#    of the seqtools.pl file
#
# Note the use of eval{} for error handling in process_seq()
#
# MODIFIED:
#  sac, 16 Jun 1998: Added installation comment, require statement comments.
#
# SEE ALSO : seqs1.pl, seqs2.pl, seqs3.pl
#---------------------------------------------------------------------------

# Using seqtools.pl in the examples/blast distribution directory:
require "seqtools.pl"; 
# Proper path to seqtools.pl after you install it in your system:
#require "/share/www-data/html/perlOOP/bioperl/lib/Bio/drivers/seq/seqtools.pl";

my $MIN_LENGTH = 250;
my @filtered   = ();

use vars qw($ID $VERSION $DESC);

$ID      = 'seqs4.pl';
$VERSION = 0.1;
$DESC    = "Demonstrates sequence editing & screening using procedures in seqtools.pl";

&init_seq();

#&load_ids();  Not loading any ids since we want to process all seqs

&get_seq_objs(\&process_seq);

if(@filtered) {
    printf STDERR "\n%d sequences were filtered:\n", scalar(@filtered);
    foreach(@filtered) { print STDERR "$_\n"; }
}

&wrap_up_seq();


#-------------
sub examples {
#-------------
    <<"QQ_EG_QQ";
(using the files in this directory)

 ./$ID seq1.fasta > out.fasta

 
QQ_EG_QQ
}


#----------------
sub process_seq {
#----------------
# Screens out short sequences and any seq with a description 
# containing the word 'nucleic'. Also edits the sequence description
# to include the length of the sequence.

    my $seq = shift;

    my $id = $seq->id;
    my $desc = $seq->desc;

    # print STDERR "\nProcessing seq $id: $desc\n";

    if($seq->seq_len >= $MIN_LENGTH and not $desc =~ /nucleic/i) {
	# Note the use of eval {} for safety. We could be doing something
	# here which could cause a runtime error. If any operation
	# fails, we don't want it to blow up the whole script but we do want
	# to know about it later.
	eval{
	    $desc .= "; LENGTH = ".$seq->seq_len;
	};
	if($@) { push @filtered, "$id  $desc\n  ERROR:\n$@"; }
	else {
	    $seq->desc($desc); # set the new description.
	    &print_seq($seq);  # print it (via seqtools.pl)
	}
    } else {
	push @filtered, "$id  $desc";
    }
}

