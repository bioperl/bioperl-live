#!/usr/bin/perl  -w

# WARNING:
#
#  There is a memory leak in the stream parsing code of Bio::Tools::Blast.pm 
#  that can cause this script to run out of memory and crash when processing 
#  streams containing large numbers of reports (several thousand).
#  sac --- Tue Jul 21 15:35:56 1998.
#
# The good news is that there is now a workaround!
# Memory usage seems to be a problem only when a -signif parameter is supplied
# separately. By placing the significance criterion within a -filt_func,
# memory usage is not a problem. The cause is still under investigation.
#
# For example, the command-line argument of:
#      -signif 1e-5 
# is equivalent to
#      -filt_func '$hit->signif <= 1e-5'
#
# For more information, see the documentation for the Blast.pm module.

#---------------------------------------------------------------------------
# PROGRAM : parse_multi.pl
# PURPOSE : To parse a set of Blast report files.
# AUTHOR  : Steve Chervitz
# CREATED : 21 Jul 1998
# REVISION: $Id$
# WEBSITE : http://bio.perl.org/Projects/Blast/
# USAGE   : parse_multi.pl -h
# EXAMPLES: parse_multi.pl -eg
#
# INSTALLATION: 
#    Set the require ".../blast_config.pl" to point to the proper location
#    of the blast_config.pl file. See blast_config.pl for additional steps.
#
# COMMENTS:
#
# This is an alternative to parse_stream.pl for parsing multiple Blast reports
# on a file-by-file basis. parse_stream.pl tends to run into excessive memory 
# usage problems when processing many (>1000) reports. Memory requirements
# of parse_multi.pl should be more stable and lower than those of 
# parse_stream.pl.
#
# This script is essentially a union of print_blasts.pl and parse.pl.
#
# MODIFIED: 
#     
# TODO: Possibly add a recurse option.
#
#---------------------------------------------------------------------------

# Using blast_config.pl in the examples/blast distribution directory:
require "blast_config.pl"; 
# Proper path to blast_config.pl after you install it in your system:
#require "/share/www-data/html/perlOOP/bioperl/bin/blast/blast_config.pl";

# Using vars from blast_config to prevent warning messages under -w.
use vars qw($ID $VERSION $DESC $MONITOR %blastParam 
	    $opt_in $opt_table $opt_compress );

$ID      = 'parse_multi.pl';
$VERSION = 0.01;
$DESC    = "Demonstrates parsing multiple Blast report files using Bio::Tools::Blast.pm";

@errs = ();
$count = 0;

select(STDOUT); $|=1;

$ID        = 'parse_multi.pl';
$opt_d     = undef;
$opt_incl  = undef;
$opt_excl  = undef;
$opt_done  = undef;
$opt_tar   = 0;
$opt_mon   = 1;
@archives = ();
$blast_obj = undef;

#----------------------
sub parse_multi_usage {
#----------------------
    print STDERR "$ID, v$VERSION\n$DESC.\n";
    print STDERR <<"QQ_USAGE_QQ";

Usage: $ID -d DIR [-h]

 -d <dir>       : Directory containing Blast reports to be parsed 
                  Reports in the dur may be compressed. 
                  (the '-d' flag is optional if dir is the last argument) 
 -incl <string> : Include only file names containing string.
                  string can be a 'regexp' (case-sensitive).
 -excl <string> : Exclude file names contain string.
                  string can be a 'regexp' (case-sensitive).
 -tar           : Access files contained in gzipped, tar archives.
 -done <dir>    : Directory for where to move Blast reports that 
                  have been processed. mkdir is called if necessary.

QQ_USAGE_QQ

    print STDERR "Hit <RETURN> to view Blast parsing parameters."; <STDIN>;
    &blast_parse_params;
    &blast_general_params;
}


&init_blast(\&parse_multi_usage, 'd=s', 'h!', 'mon!', 'incl=s', 'excl=s', 'tar!', 'done=s');
&set_blast_params;

$opt_d ||= $ARGV[0] || die &parse_multi_usage;
$opt_mon and print STDERR "\nParsing files in directory $opt_d\n";
if ($opt_done and not -d $opt_done) {
    system('mkdir', $opt_done) == 0 or
    die "Can't create directory $opt_done: $!\n";
}

if($opt_mon) {
    $opt_incl and print STDERR "Including only files with names containing: $opt_incl\n";
    $opt_excl and print STDERR "Excluding files with names passing: $opt_excl\n";
    $opt_done and print STDERR "Moving processed files to $opt_done\n";
    $opt_tar and print STDERR "Accessing any tar archives\n";
    print STDERR "\n";
}


#------------
sub examples {
#------------
<<"QQ_EG_QQ";
(Run these in the examples/blast/ directory of the distribution.)

  ./$ID -d out -incl blastp -p 1e-10 -table 2 > blast.out
  ./$ID -d out -incl blastx -done out/done -table 1 -log out/parse.log > blast.out2

QQ_EG_QQ
}

&print_blast_params() if $opt_mon;

# Process the supplied directory.
opendir (DIR, $opt_d) or die "Can't opendir $opt_d: $!\n\n";
&get_files($opt_d);
closedir DIR;

# Process any archives contained in the supplied directory.
if($opt_tar and @archives) {
    $opt_tar = 0;  # Not recursively exploring archives.
    $opt_mon and print STDERR "\nParsing files in archives:\n";
    $opt_d .= '/' unless $opt_d =~ /\/$/;
    my ($newdir, $arch);
    foreach $arch (@archives) {
	$opt_mon and print STDERR "  $arch\n";
	$newdir = $arch;
	if(not $newdir =~ s/\.tar\.gz//) {
	    print STDERR "Bad dir name: $newdir\n";
	    next;
	}
	system("(cd $opt_d; gzip -cd $arch | tar -xf -)")== 0 or die "Can't access $opt_d or $arch: $!\n";
	opendir (DIR, $newdir) or die "Can't opendir $newdir: $!\n\n";
	&get_files($newdir);
	closedir DIR;
	# Removed the opened archive.
	system("rm -rf $newdir &")== 0 or die "Can't remove opened archive dir $newdir: $!\n";
    }
}

exit 0;

#--------------
sub get_files {
#--------------
    my $d = shift;
    
    my ($f, $full);

    foreach $f (readdir DIR) {
	$full = "$d/$f";
	next if ($full =~ /\/\..*$/ or -d $full or not -s $full);
	# Extra test to process only certain files:
	next if $opt_incl and $f !~ /$opt_incl/o;
	next if $opt_excl and $f =~ /$opt_excl/o;
	
#	$opt_mon and print STDERR "$f\n";
	
	if( $full =~ /\.tar\.gz/) {
	    push @archives, $full if $opt_tar;
	    next;
	}
	&get_blast($full);

	# Check to see if the file is to be moved.
	next unless $opt_done;  
	system('mv', '-f', $full, $opt_done) == 0 or
	    die "Can't move file to $opt_done: $!\n";
    }
}



#--------------
sub get_blast {
#--------------
    my ($file) = shift;

    # Load the file into the Blast parameters.
    $blastParam{-file} = $file;

#    print "\nCREATING BLAST OBJECT FOR $file";<STDIN>;

    eval { 
	# Create the Blast object with the specified parameters.
	$blast_obj = &create_blast;  
#	$blast_obj->display();<STDIN>;
	&print_table($blast_obj);
	$opt_compress && $blast_obj->compress_file; 
	$blast_obj->destroy();  # important when crunching lots of reports.
	undef $blast_obj;
	$count++;
	};
    if($@) {
	my $er = "\nFILE: $blastParam{-file}\n$@\n";
	push @errs, $er;
	print STDERR "\n$er\n";
    }
    $opt_mon and print STDERR ".", $count % 50 ? '' : "\n";

}
