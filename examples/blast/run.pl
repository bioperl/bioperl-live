#!/usr/local/bin/perl -w

#---------------------------------------------------------------------------
# PROGRAM  : run.pl
# PURPOSE  : To demonstrate how to run and parse Blast reports using the 
#            Bio::Tools::Blast.pm module.
# AUTHOR   : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED  : 3 Feb 1998
# REVISION : $Id$
# WEBSITE  : http://bio.perl.org/Projects/Blast/
# USAGE    : run.pl -h
# EXAMPLES : run.pl -eg
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
#  sac, 16 Jun 1998: Added installation comment, require statement comments.
#                    Minor example change.
#---------------------------------------------------------------------------

# Using blast_config.pl in the examples/blast distribution directory:
require "blast_config.pl"; 
# Proper path to blast_config.pl after you install it in your system:
#require "/share/www-data/html/perlOOP/bioperl/bin/blast/blast_config.pl";

require Bio::Seq;

# Using vars from blast_config to prevent warning messages under -w.
use vars qw($ID $VERSION $DESC $MONITOR %runParam %blastParam 
	    $opt_seq $opt_seqfmt $opt_html $opt_parse 
	    $opt_dna $opt_table $opt_compress);

$ID      = 'run.pl';
$VERSION = 0.1;
$DESC    = "Demonstrates running a single Blast report using Bio::Tools::Blast.pm";

#--------------
sub run_usage {
#--------------
    print STDERR "$ID, v$VERSION\n$DESC.\n";
    print STDERR <<"QQ_USAGE_QQ";

Usage: run.pl [ parameters ] seq.file

run.pl obtains the query sequence to be Blasted from a 
file specified on the command line. 
For running multiple sequences, see blast_seqs.pl.

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

  ./$ID seq/yel009c.fasta -db swissprot -signif 1e-15 -html > out.gap2
  ./$ID seq/ymr284w.fasta -nogap -db yeast > out.nogap2
  ./$ID seq/yel009c.fasta -vers 1 -signif 1e-20 -nomon > out1
  ./$ID -html -noparse  seq/yel009c.fasta 

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
$qseq = new Bio::Seq ( -file => $opt_seq,
			  -ffmt => $opt_seqfmt,
			  -type => $opt_dna ? 'Dna' : 'Amino');

#print "\nSEQ ${\$qseq->id} IN FASTA FORMAT:\n";
#print $qseq->layout('Fasta'); <STDIN>;

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
$opt_table ? &print_table : &show_results;  # provided by blast_config.pl

$opt_compress && (print STDERR "\nCompressing file.", $blastObj->compress_file);

&wrap_up_blast;




