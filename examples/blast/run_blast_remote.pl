#!/usr/bin/perl -w

# IMPORTANT NOTE:
#  After running this script, you will have a file containing
#  NCBI BLAST queue submission data. These can be read by 
#  retrieve_blast.pl to fetch the actual BLAST reports from NCBI.
#  See comments in retrieve_blast.pl for usage details.

#---------------------------------------------------------------------------
# PROGRAM  : run_blast_remote.pl
# PURPOSE  : To demonstrate how to run and parse Blast reports using the 
#            Bio::Tools::Blast.pm module.
# AUTHOR   : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED  : 3 Feb 1998
# REVISION : $Id$
# WEBSITE  : http://bio.perl.org/Projects/Blast/
# USAGE    : run_blast_remote.pl -h
# EXAMPLES : run_blast_remote.pl -eg
#
# INSTALLATION: 
#    Set the require ".../blast_config.pl" to point to the proper location
#    of the blast_config.pl file. See blast_config.pl for additional steps.
#
# This is basically identical to the parse.pl script 
# with additional options added for running a Blast report.
#
# This demo script does not exercise all of the functionality of the Blast object.
# See the POD for the Bio::Tools::Blast.pm, accessible from the above URL.
#
# MODIFIED:
#
#   For all further modifications, refer to the cvs log.
#
#  sac, 22 Mar 2000: Added usage note about using retrieve_blast.pl
#  sac, 16 Jun 1998: Added installation comment, require statement comments.
#                    Minor example change.
#---------------------------------------------------------------------------

# Using blast_config.pl in the examples/blast distribution directory:
require "blast_config.pl"; 
# Proper path to blast_config.pl after you install it in your system:
#require "/home/steve/bin/bioperl/examples/blast/blast_config.pl";

require Bio::Seq;
require Bio::SeqIO;
use Carp;

# Using vars from blast_config to prevent warning messages under -w.
use vars qw($ID $VERSION $DESC $MONITOR %runParam %blastParam 
	    $opt_seq $opt_seqfmt $opt_html $opt_parse 
	    $opt_dna $opt_table $opt_compress);

$ID      = 'run_blast_remote.pl';
$VERSION = 0.3;
$DESC    = "Demonstrates running a single Blast report using Bio::Tools::Blast.pm";

#--------------
sub run_usage {
#--------------
    print STDERR "$ID, v$VERSION\n$DESC.\n";
    print STDERR <<"QQ_USAGE_QQ";

Usage: run_blast_remote.pl [ parameters ] seq.file

run_blast_remote.pl obtains the query sequence to be Blasted from a 
file specified on the command line. 
For running multiple sequences, see blast_seqs.pl.

NOTE: After using run_blast_remote.pl, use the helper script 
      retrieve_blast.pl to get the actual Blast report from NCBI. 
      See documentation in that script for usage details.

QQ_USAGE_QQ

    print STDERR "Hit <RETURN> to view parameters."; <STDIN>;
    &blast_run_params;
    &blast_general_params;
}

#------------
sub examples {
#------------
<<"QQ_EG_QQ";
(Run these in the examples/blast/ directory of the distribution.)

  ./$ID seq/yel009c.fasta -prog blastp -db swissprot -signif 1e-15 -html > out.gap2
  ./$ID seq/ymr284w.fasta -prog blastp -nogap -db yeast > out.nogap2
  ./$ID seq/acc1.dna.fasta -prog blastx -db ecoli -vers 1 -signif 1e-20 -nomon > out1
  ./$ID -html -noparse -prog blastp -db nr seq/yel009c.fasta 

NOTE: After using run_blast_remote.pl, use the helper script 
      retrieve_blast.pl to get the actual Blast report from NCBI. 
      See documentation in that script for usage details.

QQ_EG_QQ
}

##### MAIN ######

&init_blast(\&run_usage );

if(not -s ($opt_seq ||= $ARGV[0]))  {
    die "\nPlease specify a file containing a sequence.\n",
    "    Use the -seq option or list the file as the last arg.\n\n";
}
    

&set_blast_params();

# Create a new sequence object (query sequence)
# Only using the first sequence from the input file.

my $in = Bio::SeqIO->new('-file' => $opt_seq,
                         '-format' => $opt_seqfmt,
                         # $opt_dna is no longer used.
                         # Automatic sequence type detection is still
                         # an open issue.
                         # -type => $opt_dna ? 'Dna' : 'Amino' 
                        );

my $qseq = $in->next_seq();

$qseq or croak "\nTrouble reading sequence file.\nCannot proceed.\n";

# For debugging:
# my $out = Bio::SeqIO->new('-fh' => \*STDOUT, 
#                           '-format' => 'embl');
# print "\nSEQ ", $qseq->id, " IN EMBL FORMAT:\n";
# print $out->write_seq($qseq); 


if($qseq->id =~ /^no_id/i) {
    $opt_seq =~ /\/([\w_.]+)$/;
    $qseq->id($1 || $opt_seq);
}

# Load the sequence object in to the run parameters
# and then load the run parameters into the Blast object
# parameters (used for constructing the Blast object).
$runParam{-seqs}  = [ $qseq ];
$blastParam{-run} = \%runParam;


# Create the Blast object with the specified parameters.
# This will launch the analysis and then parse it (if desired).
# Note this uses the global $blastObj from blast_config.pl.

eval { 
    $blastObj = &create_blast;  # provided by blast_config.pl

    # The Blast object is failing to be created because
    # the NCBI request no longer returns the report directly
    # but rather an HTML page instructing how to retrieve it
    # from the queueing system. 
    # The Blast object isn't yet smart about this, so we do 
    # a simple check to see if we got a good looking blast object.
    ref $blastObj or exit(0);

    $file = $blastObj->file;
    $MONITOR and print STDERR "Blast output saved to: $file\n\n";
    
    if($opt_html and not $opt_parse) {
	
        # Uncomment the following 4 lines to send the HTML version to STDOUT.
        #	 $file = $opt_parse ? "$file.html" : $file;
        #	 open( HTML, $file) or die "Can't open HTML file: $file: $!\n\n";
        #	 while(<HTML>) { print; }
        #	 close HTML;
        exit 0;
    }
};

$@ and die $@;

$MONITOR = 0;
# Run some functions provided by blast_config.pl:
$opt_table ? &print_table($blastObj) : &show_results($blastObj);

$opt_compress && (print STDERR "\nCompressing file.", $blastObj->compress_file);

&wrap_up_blast;




