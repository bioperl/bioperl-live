BEGIN {
        use vars qw($INSTALL_PATH);

        #### NOTE: $INSTALL_PATH needs to be customized
        ###         your system. It should point
        ###         to the directory right above Bio/ in order
        ###         for perl to be able to locate the .pm files. 

	$INSTALL_PATH = "/home/steve/perl/bioperl";

        ###
        ####
    }

#---------------------------------------------------------------------------
# LIBRARY      : blast_config.pl
# PURPOSE      : To provide a set of standard functions & command-line options 
#                processing for working with Bio::Tools::Blast.pm objects.
# AUTHOR       : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED      : 15 May 1998
# REVISION     : $Id$
# WEBSITE      : http://bio.perl.org/Projects/Blast/
# INSTALLATION : Edit $INSTALL_PATH to point to the directory containing
#                your Bioperl modules (the Bio/ directory).
#
# USAGE: This script provides a library of handy routines that can be used
#        by other scripts that need to run and/or parse Blast reports.
#        All of the scripts in the examples/blast directory of the Bioperl
#        central distribution do this. See them for examples.
#
#        To use, simply require blast_config.pl into your script.
#        This script uses Perl's Getopt::Long module for processing command-line
#        parameters.
#        The -h and -eg command-line options provide usage and examples, respectively.
#
#  Here is a minimal script that uses this package:
#
#  #!/usr/bin/perl -w
#
#  ## adjust path according to your installation:
#  require "/share/www-data/html/perlOOP/bioperl/bin/blast/blast_config.pl"; 
#            
#  &init_blast();
#  &set_blast_params();
#  $blastParam{-file} = $ARGV[0];
#  &create_blast();
#  &show_results();
#  &wrap_up_blast();
#
# PUBLIC FUNCTIONS DEFINED IN blast_config.pl:
#
#    &init_blast
#    &set_blast_params 
#    &blast_usage 
#    &blast_params
#    &blast_run_params
#    &blast_parse_params 
#    &create_blast 
#    &parse_stream
#    &display
#    &file_manip 
#    &wrap_up_blast
#
# EXAMPLES : See the files parse.pl, parse_multi.pl, parse_stream.pl, run.pl,
#            and html.pl in this dir for some working examples.
#
# SCREENING HITS BASED ON ARBITRARY CRITERIA:
#
# When parsing Blast reports, the -signif parameter can be used to screen based on 
# statistical significance (P or Expect value). But hits can also be screened using 
# arbitrary criteria. This is accomplished by supplying a function reference
# in the -filt_func %blastParam:
#    $blastParam{-filt_func} = 
#                  sub { $hit=shift;
#                        $hit->gaps == 0 and $hit->frac_conserved > 0.5; };
#
# The idea is that the hit has to pass through this filter to be considered
# significant. Refer to the documentation for Bio::Tools::Blast::Sbjct.pm 
# for other ways to work with hit objects.
#
# MODIFICATIONS:
#   2 Feb 1999, sac:
#      * Removed default values for the -prog and -db command-line parameters.
#        These are no longer optional when running Blasts. Removed the -dna option
#        which was just a convenience for setting a default $opt_prog.
#
#   4 Sep  4 1998, sac:
#      * Fixed support for the -filt_func option for supplying a custom
#        filtering function via the command-line that is eval-ed to screen 
#        Blast hits.
#
#   21 Jul 1998, sac: 
#      * Changed the -err option to -log (more apropos).
#        (the -err option is still available for backward compatibility).
#      * Added the -exponent option.
#      * Added the method print_blast_params().
#
#   16 Jun 1998, sac: Added alternate $INSTALL_PATH for testing purposes.
#
#---------------------------------------------------------------------------

use lib $INSTALL_PATH;
use lib '.','..';  # fail-safe incase you forget to edit $INSTALL_PATH.

use Bio::Tools::Blast qw(:obj);
use Bio::Root::Global qw(:devel);
use Getopt::Long;
use Carp;

select(STDOUT); $|=1;

# Global:
$ID         = 'blast_config.pl';
$VERSION    = $BLAST_CONFIG_VERSION   = 0.23;
$DESC       = "Provides standard functions and variables for working with Bio::Tools::Blast.pm";
#$blastObj   = undef;  # The object created by create_blast().  # NO LONGER USED
%blastParam = ();
%runParam   = ();
@objects    = ();     # for saving all Blast objects created.
@nohits     = ();     # for saving all Blasts without hits.
@argv       = @ARGV;  # copy ARGV before GetOptions() massacres it.
$countBlast = 0;

# Command-line options:
# Parsing blast opts:
$opt_signif    = undef;  # -signif is used for screening Blast reports for significant hits
$opt_filt_func = undef;  # function reference to be executed for filtering out hits
                         # based on arbitrary criteria.
$opt_p      = undef;  # based on the significance value reported in the hit description lines.
$opt_e      = undef;  # -p and -e are provided for convenience (same as using -signif).
$opt_table  = undef;
# Booleans:
$opt_check_all = 0;   # check all hits for significance 
                      # (default = false, stops parsing after first non-signif hit).
$opt_stats     = 1;
$opt_best      = 0;
$opt_log       = 0;
$opt_err       = 0;  # deprecated (same as opt_log)
$opt_html      = 0;
$opt_stream    = 0;
$opt_share     = 1;
$opt_desc      = 0;  # do not include sbjct descriptions in table output.
$opt_wait      = undef;  # Amount of time to wait during read() before timing out.
                         # (you should need to alter this, but it's here just in case.)

# Run blast opts:
$opt_prog   = '';    # Blast flavor: blastp, blastn, tblastn, blastx, or tblastx
$opt_loc    = 0;     # Boolean = run locally requires customization of 
                     #           Bio::Tools::Blast::LocalBlast.pm for your site.
$opt_rem    = 1;     # Boolean = run remote Blast using Webblast.pm (default).
$opt_seq    = undef;
$opt_db     = '';    # The database to be queried.
$opt_seqfmt = 'fasta';
$opt_expect = 10;
$opt_v      = 100;
$opt_b      = 100;
$opt_gap_c  = undef;
$opt_gap_e  = undef;
$opt_word   = undef;
$opt_vers   = 2;
$opt_gap    = 1;
$opt_mat    = '';
$opt_filt   = 'default';
$opt_email  = undef;

# General options.
$opt_parse  = 1;    # parse the Blast report.
$opt_compress = 0;
$opt_min_len  = 15; # query sequences shorter than this are not processed.
$opt_mon    = 1;
$opt_debug  = 0;
$opt_strict = 0;
$opt_exponent = 0;
$opt_h      = !@ARGV;
$opt_eg     = 0;
$opt_params = 0;

my $_printed_labels  = 0;
my $filt_func     = undef;  # function ref based on $opt_filt_func

#----------------
sub blast_usage {
#----------------

    print STDERR "$ID, v$VERSION (blast_config.pl: v$BLAST_CONFIG_VERSION)\n$DESC.\n";
    print STDERR <<"QQ_USAGE_QQ";

For PARSING BLAST reports:
 Usage: $ID [ parsing parameters ] blast.files.*
        $ID [ parsing parameters ] < blast.file  

 blast.file  : Raw Blast report file. Can be compressed.
               (STDIN should be uncompressed).

For RUNNING BLAST analyses:
 Usage: $ID -prog <type> -db <db> [optional parameters] 

 (blast_seq.pl offers an interface for Blasting Fasta sequences.)

QQ_USAGE_QQ
    print STDERR "\nHit <RETURN> for Blast parameters"; <STDIN>;
}

#----------------
sub blast_params {
#----------------
# Prints usage information for running and parsing.
    
    &blast_run_params;
    print STDERR "\nHit <RETURN> for Blast parsing parameters"; <STDIN>;
    &blast_parse_params;
}


#-----------------------
sub blast_general_params {
#-----------------------
# Prints usage information for general parameters.

    print STDERR <<"QQ_GENERAL_QQ";

 GENERAL PARAMETERS:
 --------------------
 -noparse       : Do not parse the Blast results (default = parse).
 -compress      : Compress the Blast report file (if not already compressed).
 -min_len <int> : Query sequences shorter than this length are not processed. 
 -nomon         : Supress screen progress monitor.
 -params        : Print parsing/running parameters.
 -debug         : Activate debug mode (provides many details about object
                  construction, method invocations, etc.)
 -log <file>    : Save script log and error output to file (append mode).
                  Also will save sequences without significant hits.
 -exponent      : Report only the exponent portion of P/Expect values (integer).
 -eg            : Print examples.
 -h             : Print this help/usage information.

QQ_GENERAL_QQ
}


#--------------------
sub blast_run_params {
#---------------------
# Prints usage information for running Blasts.

    print STDERR << "QQ_RUN_QQ";

 BLAST RUN PARAMETERS:
 -------------------------------
 REQUIRED:
 -prog <type> : Type of Blast to run (type = blastp|blastn|tblastn|blastx|tblastx)
 -db   <db>   : Database to Blast against. (See list below for available dbs.)

 OPTIONAL:
 -seqfmt <str>: Format of the sequence file. Can be 'raw', 'fasta', or 'gcg'
                default seqfmt = $opt_seqfmt.
 -vers <str>  : Version of the Blast program to use (1, 2)
                default version = $opt_vers. PSI blasts not yet supported. 
 -nogap       : Perform an ungapped alignment with Blast2 
                default = gapping on.
 -expect <float>: Expect value (default = $opt_expect)
 -v <int>     : Number of one-line descriptions to show (default = $opt_v)
 -b <int>     : Number of alignments to show (default = $opt_b)
 -gap_c <int> : Gap creation penalty (default = 11)
 -gap_e <int> : Gap extension penalty (default = 1)
 -word  <int> : Word length (default = 11 for blastn, 3 for all other programs)
 -mat   <str> : Substitution scoring matrix (BLOSUM62) 
 -filt  <str> : Sequence complexity filter ($opt_filt) 
 -email <str> : Email address to send results (no parsing).
 -html        : Save the HTML-formated version of the report.
 -rem         : Run Blast at a remote site (default).
 -loc         : Run Blast at a local site (requires site-specific module).
 -compress    : Compress the Blast report file.

QQ_RUN_QQ

print &_list_dbs();

}

#--------------
sub _list_dbs {
#--------------
    my $dbn = join(', ', $Blast->db_remote('n'));
    my $dbp = join(', ', $Blast->db_remote('p'));

    sprintf "%s\n%s\n%s\n\n%s\n%s\n",
    '  Available Remote Databases (Peptide)',
    '-----------------------------------',
    $dbp,
    
    'Available Remote Databases (Nucleotide)',
    '--------------------------------------',
    $dbn;
      
  }


#--------------------
sub blast_parse_params {
#---------------------
# Prints usage information for parsing Blasts.

    print STDERR << "QQ_PARSE_QQ";

 BLAST PARSING PARAMETERS:
 --------------------------
 -signif <number>
             : Significance value to be used as a cut-off value for screening 
               hits. Number = integer, float, or sci.notation ('1e-25').
               Any hit above this number is skipped (value is an Expect value 
               for Blast2; P-value for Blast1 & WashU-Blast).
 -p  <float> : (same as -signif).
 -e  <float> : (same as -signif).
               (NOTE: For screening hits based on arbitrary criteria, see
	        comments in the blast.config.pl file.)
 -check_all  : Check all hits for significance when using significance criteria
               (default = stop checking hits after the first non-significant hit).
 -filt_func  : String to be eval-ed for custom filtering of a Blast hit, 
               E.g., -filt_func '\$hit->frac_identical > 0.9'. ("\$hit" must be used).
 -strict     : Activate strict mode for the Blast object (report additional 
               warnings and be more senstitive to various potential problems).
 -nostats    : Do not extract statistical information for the Blast run.
               (e.g., matrix, filters, Karlin-Altschul parameters)
 -noshare    : Do not share stats across a set of reports when parsing a
               Blast stream (default = share stats).
 -table <int>: Prints data from Blast report in a tab-delimited line(s).
               Table type is specified by <int> = 1,2,3. Type 1 = list HSP
               separately, type 2 = sum over all HSPs (tiling).
               Type 3 demonstrates a custom table.
 -desc       : Include sbjct descriptions in table output (default: do not include).
 -best       : Only process the best hit for each Blast report.
 -wait <int> : Amount of seconds to wait before timing out when reading in Blast
               reports. Shouldn't have to worry about this but it's provided
               just in case (default = 3sec).


QQ_PARSE_QQ

}


#----------------
sub init_blast {
#----------------
    ($usage_fref, @opts) = @_;

    &GetOptions('prog=s', 'vers=s', 'signif=s', 'p=s', 'e=s', 'parse!', 
		'mon!', 'debug!', 'h!', 'strict!', 'stats!',  'table=s', 
		'best!', 'err=s', 'log=s', 'html!', 'seq=s', 'prog=s', 'db=s', 
		'seqfmt=s', 'v=s', 'b=s', 'gap_c=s', 'gap_e=s', 
		'gap!', 'word=s', 'mat=s', 'filt=s', 'email=s', 'expect=s',
		'share!', 'stream!', 'tile_hsps!', 'residues!', 'desc!',
		'compress!', 'eg!', 'rem!', 'loc!', 'check_all!', 'min_len=s',
		'filt_func=s', 'exponent!', 'params!', 'wait=s',
		@opts);
    

    $opt_h and do{
	(ref($usage_fref) =~ /CODE/) ? &$usage_fref : (&blast_usage, &blast_params, &blast_general_params);
	exit 1;
    };

    $opt_eg and do{
	eval{ print STDERR "\n$ID EXAMPLES:\n",'-'x30,"\n\n", &examples(); };
	$@ and print STDERR "Sorry. No examples available.\n";
	exit 1;
    };

    if($opt_log or $opt_err) {
	$opt_log ||= $opt_err;
	print STDERR "\nLog file: $opt_log\n" if $opt_log;
	system('touch', $opt_log);
	open (STDERR, ">>$opt_log") or croak "Can't open log file $opt_log: $!\n\n";
    }

    if($MONITOR) {
	print STDERR "$ID, v$VERSION (blast_config.pl: v$BLAST_CONFIG_VERSION)\n",'-'x50,"\n";
	print STDERR "Started: ", `date` if $opt_stream;
    }

    if($opt_filt_func) {
	if (ref $opt_filt_func) {
	    $filt_func = $opt_filt_func;
	} else {
	    $filt_func = sub { my $hit=shift;
			       eval $opt_filt_func; };
	}
    }

    monitor($opt_mon);
    debug($opt_debug);
    
    $opt_signif = ($opt_signif || $opt_p || $opt_e);

#    $opt_prog ||= $opt_dna ? 'blastn' : 'blastp';  # -prog must now be specified when running a Blast.
    
    # -loc overrides -rem.
    $opt_loc and $opt_rem = 0;
    
    if($opt_parse) {
	# If user wants to parse it AND save an HTML-formated version
	# we need to save both versions.
	$opt_html = 'both' if $opt_html;
	
    } elsif($opt_html) {
	print STDERR "\nHTML-only mode. Blast results will not be parsed.\n\n";
    }

    # don't really need to do this since it can be set during parsing.
    $opt_exponent and $Blast->signif_fmt('exp');  
}


#---------------------
sub set_blast_params {
#---------------------
# There are quite a lot of different parameters.
# See documentation in Bio::Tools::Blast.pm for run() and parse().

# (Parameter tags may be in upper or lowercase.
#  I find uppercase to be more readable and safer
#  but lowercase seems more standard).

# Define parameters for the Blast run.
%runParam = ( 
	      -method   => $opt_rem ? 'remote' : 'local',
	      -prog     => $opt_prog,
	      -version  => $opt_vers,
	      -database => $opt_db,
	      -html     => $opt_html,
	      -seqs     => [ ],
	      -descr    => $opt_v,
	      -align    => $opt_b,
	      -expect   => $opt_expect,
	      -gap      => $opt_gap ? 'on' : 'off',
	      -filter   => $opt_filt,
	      -matrix   => $opt_mat,
	      -email    => $opt_email,
	      -gap_c    => $opt_gap_c,
	      -gap_e    => $opt_gap_e,
	      -word     => $opt_word,
	      -min_len  => $opt_min_len,
	      );     

# Define parsing parameters for the Blast object.
%blastParam = ( 
		-run        => '',
		-file       => '',
		-parse      => $opt_parse,
		-signif     => $opt_signif, 
		-filt_func  => $filt_func, 
		-min_len    => $opt_min_len, 
		-check_all_hits => $opt_check_all,
		-strict     => $opt_strict,
		-stats      => $opt_stats,
		-best       => $opt_best,
#		-stream     => $opt_stream,   # No longer used.
		-share      => $opt_share,
		-signif_fmt => $opt_exponent,
		-exec_func  => '',
		-save_array => '',   # Not used if exec_func is defined
                                     # (experimental).
  	        '-wait'     => $opt_wait,  # Use -WAIT and save the quotes
	       );

    &print_blast_params if $opt_params;
}

#-----------------------
sub print_blast_params {
#-----------------------
    my ($runp);

    print STDERR "\nBlast Parse Parameters:\n", '-'x25,"\n";
    foreach(sort keys %blastParam) {
	if($_ eq '-run' and ref($blastParam{$_})) { 
	    %runp = %{$blastParam{$_}}; 
	    next; 
	}
	if($_ eq '-filt_func') {
	    printf STDERR "  %15s => %s\n", $_, defined($blastParam{$_}) ? $opt_filt_func : '';
	} else {
	    printf STDERR "  %15s => %s\n", $_, defined($blastParam{$_}) ? $blastParam{$_} : '';
	}
    } 

    if(%runp) {
	print STDERR "\nBlast Run Parameters:\n", '-'x25,"\n";
	foreach(sort keys %runp) {
	    printf STDERR "  %15s => %s\n", $_, defined($runp{$_}) ? $runp{$_} : '';
	} 
    }
    print STDERR "\n", '-'x50,"\n\n";
}



#----------------
sub create_blast {
#----------------
    my $blast_obj;

    if(scalar $blastParam{'-run'} and not ($opt_prog and $opt_db)) {
      print STDERR "\nBlast program and/or database not defined.\n";
      my $msg = '';
      if(not $opt_prog)  {
	$msg = sprintf "Please defined a -prog parameter of\n".
	   "  blastp | blastn | tblastn | blastx | tblastx\n\n";
      }
      if(not $opt_db) {
	$msg .= &_list_dbs();
      }
      die $msg;
    }

    eval { 
	$countBlast++;
	$blast_obj = new Bio::Tools::Blast (%blastParam);
    };
    if($@) {
	croak "\n*** Trouble creating Blast object:\n$@\n\n";
    }
    return $blast_obj;
}

#----------------
sub parse_stream {
#----------------
    # Using the static Blast object exported by Bio::Tools::Blast.pm
    $countBlast = $Blast->parse(%blastParam);
}


#----------------
sub print_table {
#----------------
# Purpose  : Provides an interface to the table() functions of the Blast object.
#            Print a tab-delimited set of data for each significant hit.
#            Objects without significant hits are recorded.
# Usage    : &print_table( [blast_obj], [table_type])
# Argument : blast_obj  = Blast.pm object reference
#            table_type = integer: 1, 2, or 3 (default = $opt_table or 1)
#                         Type 1: Print data for each HSP on separate lines
#                                 by calling Blast::table().
#                         Type 2: Print tiled data for each hit (sums across HSPs)
#                                 by calling Blast::table_tiled().
# Returns  : n/a, prints line to STDOUT.
# Comments : Uses the $opt_desc value as an argument to the table() calls on
#            the Blast object.

    my $bo   = shift || die "blast_config::print_table(): No Blast object was supplied at line ".__LINE__."\n";
    local $_ = shift || $opt_table || 1;

    if($bo->is_signif) {
	SWITCH: {
	    /1/ && do{ print $bo->table_labels($opt_desc) if not $_printed_labels; 
		       print $bo->table($opt_desc); last SWITCH;};
	    /2/ && do{ print $bo->table_labels_tiled($opt_desc) if not $_printed_labels; 
		       print $bo->table_tiled($opt_desc); last SWITCH;};
	    /3/ && do{ print $_table_custom($bo); last SWITCH;};
	}
    } else {
	push @nohits, $bo->name.", File: ".($bo->file || '<STDIN>');
#	print STDERR "No significant hits: ${\$bo->name}\n";
    }
    $_printed_labels = 1;

    $bo->destroy if $opt_stream;  # important when crunching lots of reports.
    
}


#------------------
sub _table_custom {
#------------------
# Shows how to generate a custom set of Blast data when the standard
# table output methods provided by the Blast object (table() and table_tiled(), q.v.)
# won't do.
# Returns : string containing tab-delimited set of data for each hit.

    my $bo  = shift;
    my $str = '';

    $str = &_table_custom_labels unless $_custom_labels;

    foreach $hit($bo->hits) {
	$str .= printf "%s\t%s\t%.1e\t%.2f\t%s\n", 
	               $bo->name, $hit->name, $hit->expect, 
	               $hit->frac_identical, $hit->desc;
    }
    $str;
}

#-------------------------
sub _table_custom_labels {
#-------------------------
# Return   : string containing  column labels. Used in conjunction with table_custom()
    my $str = sprintf "\n%s\t%s\t%s\t%s\t%s\n", 
                      'QUERY', 'SBJCT', 'EXPCT', 'FR_ID', 'DESC';
   $str .= sprintf "%s\t%s\t%s\t%s\t%s\n", 
                      '-----', '-----', '------', '-----', '-----';

    $_custom_labels = 1;    
    $str;
}


#----------------
sub show_results {
#----------------
    my $bo  = shift || die "blast_config::show_results(): No Blast object was supplied.\n";

    $opt_parse or return print $bo->read();

    if (not $bo->is_signif) {
	print STDERR "\nNo significant hits.\n";
	$bo->display();
	return;
    }
    
    printf "\nQuery: %s, Length = %d\n", $bo->name, $bo->length;
    printf "%d total hits with significance <= %s\n\n", $bo->num_hits, 
    $bo->signif;
    
    ## Display stats and standard information for each hit.
    ## You need a wide screen to see this properly.
    $bo->display(); 

    $MONITOR && do{ print STDERR "\nHit <RETURN> to see hits."; <STDIN>; };
    $bo->display(-show=>'hits');
    
    &display_hit_info($bo);
}


#--------------------
sub display_hit_info {
#--------------------
## Create a custom display of information about each hit.
    my $bo  = shift || die "blast_config::display_hit_info(): No Blast object was supplied at line ".__LINE__."\n";

    $opt_parse or return print $bo->read();

    print "\n\nCUSTOM LISTING.\n";
    printf( "\n%-5s %-20s %-7s %-9s %-3s %-5s %-7s %-7s %-7s %-4s %s\n", 'HIT#', '   ID', 'SCORE', 'EXPECT', 'N', 'GAPS', 'F_ID', 'F_CONS', 'F_ALQ', 'QLEN', 'SLEN');
    printf( "%-5s %-20s %-7s %-9s %-3s %-5s %-7s %-7s %-7s %-4s %s\n", '----', '-------', '----', '-------', '--', '----', '-----', '-----', '-----', '----', '----');
    
    my $hit_count = 0;
    foreach $hit ($bo->hits) {
	
	# Skipping self-hit based on sequence identifier only.
	next if $hit->name eq $bo->name;
	
	$hit_count++;
	
	# Print out some data for the hit: 
	printf( "%-5d %-20s %-7d %-9.1e %-3d %-5d %-7.2f %-7.2f %-7.2f %4d %4d \n",
		$hit_count,
		$hit->name(),
		$hit->score(),
		$hit->expect(),
		$hit->num_hsps(),
		$hit->gaps(),
		$hit->frac_identical(), 
		$hit->frac_conserved(), 
		$hit->frac_aligned_query(), 
		$bo->length(), 
		$hit->length(), 
		);
	
#   Print HSP data:
	my $hsp_count = 0;
	my (@hsp_datum);
	foreach $hsp ($hit->hsps) {
	    $hsp_count++;
	    @hsp_datum = ($hsp->range('query'), $hsp->range('sbjct'));
	    push @hspData, [ $hit_count, $hsp_count, @hsp_datum];
#		printf "\nrecording hsp datum (n= %d):\n   @hsp_datum", scalar @hsp_datum;<STDIN>;
	}
    }
    print "\n\nHSP DATA FOR ALL HITS:\n";
    printf "%6s %6s %12s %12s\n", "HIT#", "HSP#", "QUERY-RANGE", "SBJCT-RANGE";
    printf "%6s %6s %12s %12s\n", "----", "----", "-----------", "------------";
    
    foreach (@hspData) {
	printf "%4d %6d %7d - %-5d %5d - %-5d\n", @$_;
    }
}    

#------------
sub display {
#------------
# Works with a set of Blast objects, when @objects is used.
    print STDERR "\n$countBlast Blast reports parsed.\n";
    foreach(@objects) {
	next unless $_->is_signif;
	print "\n\n",'-'x70,"\n";
	$_->display(-show=>'hits');
    }
}


#--------------
sub file_manip {
#--------------
    my $bo = shift || die "blast_config::file_manip(): No Blast object was supplied.\n";

    print STDERR "\n\nBASIC FILE MANIPULATIONS:\n";
    print STDERR "----------------------------\n\n";
    
    my ($file);
    printf STDERR "Current Blast file: %s\n", $bo->file(); 

    if($file = $bo->compress_file()) {
	printf STDERR "\nCompressed Blast file: %s\n", $file;
	$file = $bo->uncompress_file();
    } elsif($file = $bo->uncompress_file()) {
	printf STDERR "\nUncompressed Blast file: %s\n", $file;
	$file = $bo->compress_file();
    }
    print STDERR "\nLeaving file in its original state: $file\n";
}


#------------------
sub wrap_up_blast {
#------------------
# Report/record which Blast reports had no significant hits.

    if(@nohits) {
	print "\nNo significant hits for the following sequences:\n";
	my $count = 0;
	foreach(@nohits) {
	    $count++;
	    printf " %-10s", $_;
	    $count%6==0 and print "\n";
	}
	if($opt_log) {
	    my $file = "no_hits.$$.out";
	    open( NOHITS, ">$file") || croak "\nTrouble saving Blast reports without significant hits.\n".
		"Can't open file $file: $!\n";
	    print NOHITS join("\n", @nohits);
	    close NOHITS;
	    $MONITOR && do{
		printf STDERR "\n%d reports had no significant hits.\n", scalar(@nohits);
		print STDERR "saved to $file\n";
	    };
	}
    }
    $MONITOR and do{
	print STDERR "\n$countBlast Blast report(s) processed.\n";
	print STDERR "Finished parsing at ", `date` if $countBlast > 50;
	print STDERR "\nDone. $ID @argv\n";
    };

    exit 0;
}



1;
