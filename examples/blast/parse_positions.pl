#!/usr/bin/perl -w

#---------------------------------------------------------------------------
# PROGRAM : parse_positions.pl
# PURPOSE : To extract identical + conserved positions for a hit.
# AUTHOR  : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED : 12 Jun 1998
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
#   Data for each HSP separately could be collected separately, if desired.
#
# Sample BLAST output files can be found in examples/out/.
#
# MODIFIED:
#  sac, 16 Jun 1998: Added installation comment, require statement comments.
#---------------------------------------------------------------------------

# Using blast_config.pl in the examples/blast distribution directory:
require "blast_config.pl"; 
# Proper path to blast_config.pl after you install it in your system:
#require "/share/www-data/html/perlOOP/bioperl/bin/blast/blast_config.pl";

# Using vars from blast_config to prevent warning messages under -w.
use vars qw($ID $VERSION $DESC $MONITOR %blastParam 
	    $opt_in $opt_table);

$ID      = 'parse_positions.pl';
$VERSION = 0.01;
$DESC    = "Parses a Blast reports to extract identical + conserved positions for a hit.";

$opt_hit      = 'all';
$opt_collapse = 1;

@errs = ();

#------------
sub _usage {
#------------
    print STDERR "$ID, v$VERSION\n$DESC.\n";
    print STDERR <<"QQ_USAGE_QQ";

Usage: $ID [ parameters ] blast.file
       $ID [ parameters ] < blast.file  

 blast.file  : Raw Blast report file. Can be compressed.
               (STDIN should be uncompressed).
QQ_USAGE_QQ
    print STDERR "\nHit <RETURN> for Blast parameters"; <STDIN>;

    &blast_parse_params;

print STDERR <<OPTIONS;
 -hit <name>  : Sequence identifier for a hit for which positions
                are to be extracted (default = all hits).
 -nocollapse  : Don't collapse consecutive hits into from-to strings.
              : Collapsing: 1,2,3,4,5 ---> 1-5  (default = collapse).

OPTIONS
#'
    &blast_general_params;
    print STDERR "\nThis script processes single Blast reports only.\n\n";
}

#------------
sub examples {
#------------
<<"QQ_EG_QQ";
(Run these in the examples/blast/ directory of the distribution.)

  ./$ID out/blastx.2 -signif 1e-4
  ./$ID out/blastp.1.gcg.gz -signif 1e-100 -hit YER103W

QQ_EG_QQ
}

##### MAIN #####

&init_blast(\&_usage, 'hit=s', 'sbjct!', 'query!', 'collapse!');
&set_blast_params();

my $file = $opt_in || $ARGV[0];
my (@qidentical, @sconserved, @qconserved);  
# my (@qidentical, @qidentical, @sconserved, @qconserved);   #ps 3/24/00

# Load the file into the Blast parameters.
$blastParam{-file} = $file || '';

eval { 
    my $blast_obj = &create_blast;  
#    $opt_table ? &print_table($blast_obj) : &show_results($blast_obj);
	
    print "IDENTICAL AND CONSERVED POSITIONS FROM BLAST REPORT\n";
    print "====================================================\n";
    print "FILE : ", $file || 'STDIN', "\n";
    print "QUERY: ", $blast_obj->query, "\n";
    print "HITS : ", $blast_obj->num_hits, "\n";

    $opt_hit ||= $blast_obj->hit->name;

    my ($hit, $hsp);
    foreach $hit ($blast_obj->hits) {
	next unless $opt_hit eq 'all' or $hit->name =~ /$opt_hit/i;
	@qidentical = @qconserved = @sidentical = @sconserved = ();
	# This will merge data for all HSPs together.
	push @qidentical, $hit->seq_inds('query', 'id',   $opt_collapse);
	push @qconserved, $hit->seq_inds('query', 'cons', $opt_collapse);
	push @sidentical, $hit->seq_inds('sbjct', 'id',   $opt_collapse);
	push @sconserved, $hit->seq_inds('sbjct', 'cons', $opt_collapse);
	&report_data($hit);
	last unless $opt_hit eq 'all';
    }
};
if($@) {
    my $er = "\nFILE: $blastParam{-file}\n$@\n";
    push @errs, $er;
}

#---------------
sub report_data {
#---------------
    my $hit = shift;
    print "---------------------------------------------------------------\n";
    print "SBJCT: ", $hit->name, "\n";
    print "HSPS : ", $hit->num_hsps, "\n";

    print "QUERY_IDENTICAL: ";
    print join(",", @qidentical);

    print "\n\nQUERY_CONSERVED: ";
    print join(',', @qconserved); 
    print "\n\n";
    
    print "SBJCT_IDENTICAL: ";
    print join(",", @sidentical);
    
    print "\n\nSBJCT_CONSERVED: ";
    print join(",", @sconserved);
    print "\n\n";
}


if(@errs) {
    print STDERR "\n*** Blast report produced fatal error:\n";
    foreach(@errs) { print STDERR $_; }
}


