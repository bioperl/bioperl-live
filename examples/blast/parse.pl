#!/usr/bin/perl -w

#---------------------------------------------------------------------------
# PROGRAM : parse.pl
# PURPOSE : To demonstrate parsing features of the Bio::Tools::Blast.pm module.
# AUTHOR  : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED : 3 Feb 1998
# REVISION: $Id$
# WEBSITE : http://bio.perl.org/Projects/Blast/
# USAGE   : parse.pl -h
# EXAMPLES: parse.pl -eg
#
# INSTALLATION: 
#    Set the require ".../blast_config.pl" to point to the proper location
#    of the blast_config.pl file. See blast_config.pl for additional steps.
#
# COMMENTS:
#
# Sample BLAST output files can be found in examples/blast/out/ of the distribution.
# For processing a stream of Blast reports, see the the parse_stream.pl script.
#
# This demo script does not exercise all of the functionality of the Blast object.
# See parse2.pl and parse_positions.pl script for some other manipulations and 
# the POD for the Bio::Tools::Blast.pm, accessible from the above website.
#
# MODIFIED:
#  sac,  4 Sep 1998: Added example of using -filt_func option.
#  sac, 16 Jun 1998: Added installation comment, require statement comments.
#                    Minor alteration of seq_inds() calls.
#  sac, 15 Jul 1998: Segregated code into parse2.pl which was formerly in 
#                    parse.pl but commented out.
#---------------------------------------------------------------------------

# Using blast_config.pl in the examples/blast distribution directory:
require "blast_config.pl"; 
# Proper path to blast_config.pl after you install it in your system:
#require "/share/www-data/html/perlOOP/bioperl/bin/blast/blast_config.pl";

# Using vars from blast_config to prevent warning messages under -w.
use vars qw($ID $VERSION $DESC $MONITOR %blastParam 
	    $opt_in $opt_table $opt_compress $opt_filt_func);

$ID      = 'parse.pl';
$VERSION = 0.02;
$DESC    = "Demonstrates parsing Blast reports using Bio::Tools::Blast.pm";

@errs = ();
#$opt_filt_func =
#    sub { $hit=shift;
#	   $hit->frac_aligned_hit >= 0.8; };
#

#-----------------
sub parse_usage {
#-----------------
    &blast_usage;
    &blast_parse_params;
    &blast_general_params;
}

#------------
sub examples {
#------------
<<"QQ_EG_QQ";
(Run these in the examples/blast/ directory of the distribution.)

  ./$ID out/blastx.2 
  ./$ID out/blastp.2.gz -signif 1e-15 -table 1
  ./$ID out/blastp.2.gz -signif 1e-15 -table 1 -exponent -desc
  ./$ID out/blastp.2.gz -signif 1e-15 -table 2
  ./$ID out/blastp.2.wu -check_all -filt_func '\$hit->gaps == 0' -table 2
  ./$ID out/blastp.205.gz -signif 1e-1 -nostats
  ./$ID -signif 1e-5 -table 1 < out/tblastn.2 > parsed.out
  ./$ID out/blastx.2.email.gz -table 1 -signif 1e-4  
  ./$ID out/blastp.email.html.gz -signif 1e-10 
  ./$ID out/blastp.2* -table 1 -best -nostats > parsed.out2
  ./$ID out/tblastn.206.gz -table 2 -signif 0.1 
  ./$ID out/blastp.1.gz   # should issue some warnings.

QQ_EG_QQ
}

##### MAIN #####

&init_blast(\&parse_usage);

if(!@ARGV and $opt_in) {  push @ARGV, $opt_in; }

&set_blast_params();


my ($blast_obj);

if(@ARGV) {
    # Building object(s) from files specified on command line.
    # Note that we don't really need to capture the $blast_object 
    # created by create_blast() since we can always access it via
    # the global $blastObj defined in blast_config.pl.
    # However, doing so makes things more obvious.
    $MONITOR && print STDERR "\nParsing Blast report file(s).\n";
    $count = 0;
    while($_ = shift) {
	# Load the file into the Blast parameters.
	next unless -s;
	$blastParam{-file} = $_;

	eval { 
	    # Create the Blast object with the specified parameters.
	    # Using functions provided by blast_config.pl
	    # which also supplies $blastObj.
	    $blast_obj = &create_blast;  
	    $opt_table ? &print_table($blast_obj) : &show_results($blast_obj);

	    $opt_compress && $blast_obj->compress_file; 
	    $blast_obj->destroy();  # important when crunching lots of reports.
	    $count++;
	};
	if($@) {
	    my $er = "\nFILE: $blastParam{-file}\n$@\n";
	    push @errs, $er;
	}
	print STDERR ".", $count % 50 ? '' : "\n";
    }
} else {
    # Building object from STDIN. Expecting only one Blast report.
    # To parse a stream of Blast reports, use parse_stream.pl.
    print STDERR "\nParsing Blast stream from STDIN.\n";
    $opt_table 
      # Process each Blast as you go.
      ? $blastParam{-exec_func} = \&print_table
          # Alternatively, try this:
          #? $blastParam{-exec_func} = \&display_hit_info();
      # Save all the Blast objects.
      :$blastParam{-save_array} = \@objects;  
}

if(@errs) {
    printf STDERR "\n*** %d Blast reports produced fatal errors:\n", scalar(@errs);
    foreach(@errs) { print STDERR $_; }
}

&wrap_up_blast;



