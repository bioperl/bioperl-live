#!/usr/local/bin/perl  -w

# WARNING:
#
#  There is a memory leak in the stream parsing code of Bio::Tools::Blast.pm 
#  that can cause this script to run out of memory and crash when processing 
#  streams containing large numbers of reports (several thousand).
#  sac --- Tue Jul 21 15:35:56 1998.
#
#  The memory leak has been somewhat abated but is still a problem.
#  The severity of the problem depends on the nature of the reports
#  and the parsing parameters (e.g., saving all hits or only those 
#  below 1e-20).
#  sac --- Thu Dec  3 00:22:04 1998

#---------------------------------------------------------------------------
# PROGRAM : parse_stream.pl
# PURPOSE : To demonstrate parsing a stream of Blast reports using the 
#           Bioperl Bio::Tools::Blast.pm module.
# AUTHOR  : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED : 20 Apr 1998
# REVISION: $Id$
# WEBSITE : http://bio.perl.org/Projects/Blast/
# USAGE   : parse_stream.pl -h
# EXAMPLES: parse_stream.pl -eg
#
# INSTALLATION: 
#    Set the require ".../blast_config.pl" to point to the proper location
#    of the blast_config.pl file. See blast_config.pl for additional steps.
#
# For processing single Blast reports, see the the parse.pl script.
#
# This demo script does not exercise all of the functionality of the Blast object.
# See the POD for the Bio::Tools::Blast.pm, accessible from the above URL.
#
# Sample BLAST output files: drivers/blast/out
# Sample HTML-formatted BLAST output file: drivers/blast.html
#
# MODIFIED:
#  16 Jun 1998, sac: Added installation comment, require statement comments.
#  21 Jul 1998, sac:
#     * Added warning about memory usage problem.
#     * Added print_blast_params() call.
#---------------------------------------------------------------------------

# Using blast_config.pl in the examples/blast distribution directory:
require "blast_config.pl"; 
# Proper path to blast_config.pl after you install it in your system:
#require "/share/www-data/html/perlOOP/bioperl/bin/blast/blast_config.pl";

# Using vars from blast_config to prevent warning messages under -w.
use vars qw($ID $VERSION $DESC $MONITOR %blastParam @objects $opt_table);

$ID      = 'parse_stream.pl';
$VERSION = 0.1;
$DESC    = "Demonstrates parsing a stream of Blast report using Bio::Tools::Blast.pm";
$opt_table = 1;

#----------------------
sub parse_stream_usage {
#-----------------------
    print STDERR "$ID, v$VERSION\n$DESC.\n";
    print STDERR <<"QQ_USAGE_QQ";

Usage: gzcat blast*.gz | $ID [ parameters ] > outfile
       print_blasts.pl dir | $ID [ parameters ] > outfile

 print_blasts.pl is a convenience script for working with 
 very large numbers of reports. See comments in print_blasts.pl
 for details.

 Blast reports with no significant hits will be saved to a file named
 "no_hits.#####.out", where '#####' is the PID of this script.

 WARNING:
  THIS SCRIPT SHOULD NOT BE USED FOR PROCESSING STREAMS CONTAINING
  LARGE NUMBERS OF BLAST REPORTS (N > 1000 OR SO). INSTEAD, USE 
  THE SCRIPT parse_multi.pl 

QQ_USAGE_QQ

    print STDERR "Hit <RETURN> to view parameters."; <STDIN>;
    &blast_parse_params;
    &blast_general_params;
}


#------------
sub examples {
#------------
<<"QQ_EG_QQ";
(Run these in the examples/blast/ directory of the distribution.)

  gzcat ./out/blastp.2?.gz | ./$ID -signif 1e-10 -table 2 > blast.table2
  cat ./out/blastx* | ./$ID -table 1 > blast.table1
  print_blasts.pl ./out | ./$ID -best -noshare
    
  The '-noshare' argument is necessary because the out/ directory
  contains a mixed bag of Blast reports (version 1, 2, blastp, tblastn,
  gapped, ungapped, etc.). Most of the time, -noshare is unnecessary
  since all reports have the same program, version, gapping, ... data.

  The "print_blasts.pl dir" syntax is recommended when working with large
  numbers of Blast reports (thousands). The Blasts reports located in "dir"
  can be compressed or not. 

QQ_EG_QQ
}


##### MAIN #####

&init_blast(\&parse_stream_usage );

&set_blast_params();


# Define parameters necessary for stream parsing
$blastParam{-stream}    = 1;
# Uncomment one of the two following lines:
$blastParam{-exec_func} = \&print_table;
#$blastParam{-save_array}  = \@objects;     # if you need to save all the objects.

&print_blast_params() if $MONITOR;

# Start Parsing.
$MONITOR && print STDERR "\nParsing Blast stream\n";

eval { &parse_stream;  };

if($@) {
    die "\n*** TROUBLE:\n$@\n";

} elsif(@objects) {
    &display;
}

&wrap_up_blast;

exit 0;


