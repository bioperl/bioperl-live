#!/usr/bin/perl -w

#---------------------------------------------------------------------------
# PROGRAM  : run_blast_local.pl
# PURPOSE  : To run local Blast analyses using the Bio::Tools::Blast.pm module.
# AUTHOR   : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED  : 20 Apr 1999
# REVISION : $Id$
#
# COMMENTS:
# This is just a simple prototype for running local Blasts
# on a set of sequences in a fasta file.
# It doesn't actually use the Bio::Tools::Blast system yet.
#---------------------------------------------------------------------------

# Proper path to seqtools.pl after you install it in your system:
require "/home/steve/perl/bioperl/examples/seq/seqtools.pl";

use vars qw($ID $VERSION $DESC);

$ID      = 'run_blast_local.pl';
$VERSION = 0.1;
$DESC    = "Runs Blast locally using all sequences in a Fasta file.";
$count = 0;
$blast = "/home/steve/bin/blast/blastall"; # Blast executable.

&init_seq();

#&load_ids();  Not loading any ids since we want to process all seqs

&get_seq_objs(\&blast_seq);

&wrap_up_seq();

#----------------
sub blast_seq {
#----------------
# Executes blast using a sequence object as input.

    my $seq = shift;

    $count++;

    my $id = $seq->id;
    if( $id =~ s/\|/_/g) {
      my @ids = split ('_', $id);
      $id  = $ids[0].'_'.$ids[1];
    }
    my $tmpfile = "/tmp/$id.fasta";
    my $outfile = "./blastout/$id.blastp";
    open (TMPOUT, ">$tmpfile") || die "Can't write file $tmpfile: $!\n";
    print TMPOUT $seq->out_fasta;
    close TMPOUT;

    # print STDERR "\nProcessing seq $id\n";

# Here's a typical command
#  blastall -p blastp -d yeast.aa -i tpi1.fasta -o tpi.blast.pam 

    system($blast, '-p', 'blastp', '-d', 'yeast.aa', '-i', $tmpfile, '-o', $outfile) == 0 or warn "\nTrouble executing blastp: $!\n";

}



