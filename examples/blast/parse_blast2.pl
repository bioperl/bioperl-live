#!/usr/bin/perl -w

#---------------------------------------------------------------------------
# PROGRAM : parse_blast2.pl
# PURPOSE : To demonstrate additional parsing features of the Bio::Tools::Blast.pm module.
# AUTHOR  : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED : 14 Jul 1998
# REVISION: $Id$
# WEBSITE : http://bio.perl.org/Projects/Blast/
# USAGE   : parse_blast2.pl -h
# EXAMPLES: parse_blast2.pl -eg
#
# INSTALLATION: 
#    Set the require ".../blast_config.pl" to point to the proper location
#    of the blast_config.pl file. See blast_config.pl for additional steps.
#
# COMMENTS:
#
# This demo script does not exercise all of the functionality of the Blast object.
# See the parse_blast.pl and parse_positions.pl script for some other manipulations  
# and the POD for the Bio::Tools::Blast.pm, accessible from the above website.
#
# MODIFIED:
#   30 Jul 1998, sac:
#     * Added test for retrieval of strand information from nucleotide blasts.
#
#---------------------------------------------------------------------------

# Using blast_config.pl in the examples/blast distribution directory:
require "blast_config.pl"; 
# Proper path to blast_config.pl after you install it in your system:
#require "/share/www-data/html/perlOOP/bioperl/bin/blast/blast_config.pl";

# Using vars from blast_config to prevent warning messages under -w.
use vars qw($ID $VERSION $DESC $MONITOR %blastParam 
	    $opt_in $opt_table $opt_compress);

$ID      = 'parse_blast2.pl';
$VERSION = 0.03;
$DESC    = "Demonstrates additional parsing Blast reports using Bio::Tools::Blast.pm";

@errs = ();

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

  ./$ID out/blastp.2.gz -signif 1e-15 -table 1
  ./$ID out/blastn.2.gz -table 1

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
    $MONITOR && print STDERR "\nReading Blast report from file(s).\n";
    $count = 0;
    while($_ = shift) {
	# Load the file into the Blast parameters.
	next unless -s;
	$blastParam{-file} = $_;

	# This supplies a function for filtering hits based on
	# arbitrary criteria. The idea is that the hit has to pass through this filter 
	# to be considered significant.
	# See the Bio::Tools::Blast::Sbjct.pm documentation for working with hit objects.
	  $blastParam{-filt_func} = 
	       sub { $hit=shift; 
		     $hit->num_hsps == 1 and $hit->frac_conserved >= 0.4; 
		 };
	eval { 
	    # Create the Blast object with the specified parameters.
	    # Using functions provided by blast_config.pl
	    # which also supplies $blastObj.
	    $blast_obj = &create_blast;  
	    $blast_obj->display();

	    print STDERR "Hit <RETURN> to continue."; <STDIN>;

	    $opt_table ? &print_table($blast_obj) : &show_results($blast_obj);

	     my $hsp = $blast_obj->hit->hsp;  # gets the first HSP of the first hit.

	    if($blast_obj->program =~ /blast[nx]/i) {
		print "\nQuery strand = ", $hsp->strand('query') || 'UNKNOWN', "\n";
		print "Sbjct strand = ", $hsp->strand('sbjct') || 'UNKNOWN', "\n";
		@strands = $hsp->strand();
		print "Both strands = @strands\n";
	    }
	    
	    print "\nQuery HSP seq = ", $hsp->seq_str('query'), "\n";

	     print "\nQuery HSP identical indices = \n", 
		     join(', ', $hsp->seq_inds('query', 'identical')), "\n";
	     print "\nQuery HSP conserved indices = \n", 
		    join(', ', $hsp->seq_inds('query', 'cons')), "\n";

	     print "\n\nSbjct HSP seq = ", $hsp->seq_str('sbjct'), "\n";
	     print "\nSbjct HSP identical indices = \n", 
		   join(', ', $hsp->seq_inds('sbjct', 'identical', 1)), "\n";
	     print "\nSbjct HSP conserved indices = \n", 
		    join(', ', $hsp->seq_inds('sbjct', 'cons', 1)), "\n";
	     
	     print "\nMatch sequence = \n", 
		    join(', ', $hsp->seq_str('match')), "\n";

	     print STDERR "Hit <RETURN> to continue."; <STDIN>;

	     print "\nBio::Seq objects:\n";
	     print  "Query in Fasta format:\n", $hsp->seq('query')->layout('Fasta');
	     print  "\nSbjct in GCG format:\n", $hsp->seq('sbjct')->layout('GCG');

	     print STDERR "Hit <RETURN> to continue."; <STDIN>;

	     print "\nBio::UnivAln object:\n";
	     my $aln = $hsp->get_aln;
	     print " consensus:\n", $aln->consensus();
	     print "\n layout:\n", $hsp->get_aln->layout('fasta');
	     # ReadSeq is required for MSF format.
	     $ENV{READSEQ_DIR} = '/home/users/sac/bin/solaris';
	     $ENV{READSEQ} = 'readseq';
	     print "\n layout:\n", $hsp->get_aln->layout('msf');

	     print STDERR "Hit <RETURN> to continue."; <STDIN>;

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
    print STDERR "\nReading Blast report from STDIN.\n";
    $blast_obj = &create_blast;
    $opt_table ? &print_table($blast_obj) : &show_results($blast_obj);

    # Uncomment this line for an different way to display hit data.
    #$opt_table ? &print_table($blast_obj) : &display_hit_info($blast_obj);

}

if(@errs) {
    printf STDERR "\n*** %d Blast reports produced fatal errors:\n", scalar(@errs);
    foreach(@errs) { print STDERR $_; }
}

&wrap_up_blast;



