#!/usr/local/bin/perl  -w

#---------------------------------------------------------------------------
# PROGRAM : print_blasts.pl
# PURPOSE : To print a set of Blast reports, generating a STDOUT stream
#           to be used as input by parse_stream.pl
# AUTHOR  : Steve A. Chervitz
# CREATED : 20 Apr 1998
# REVISION: $Id$
# USAGE   : print_blasts.pl -h
#
#           For modest numbers of reports, use the command line instead of 
#           print_blasts.pl as in:
#               cat blast.out.* | parse_stream.pl [args] > file
#
#           This will work fine for low to moderate numbers of Blast reports but
#           may cause the shell to complain for very large numbers of reports
#           (>500 or so), depending on your shell.
# COMMENTS:
#
# SEE ALSO: parse_stream.pl for more about processing Blast streams:
#
# MODIFIED: 
#     14 May 1998, sac.
#     20 Jul 1998, sac:
#        * Added ability to work with gzipped tar archives.
#     
# TODO: Possibly add a recurse option.
#
#---------------------------------------------------------------------------

use Getopt::Long;

select(STDOUT); $|=1;

$ID        = 'print_blasts.pl';
$opt_d     = undef;
$opt_incl  = undef;
$opt_excl  = undef;
$opt_done  = undef;
$opt_tar   = 0;
$opt_mon   = 0;
$opt_h     = !@ARGV;
@archives = ();

&GetOptions('d=s','h!','mon!', 'incl=s', 'excl=s', 'tar!', 'done=s');

$opt_d ||= $ARGV[0];

if($opt_h or !$opt_d) {
    die <<"QQ_USAGE_QQ";

Usage: $ID -d DIR [-h]

 -d DIR         : Directory containing Blast reports to be parsed 
                  Reports in the DIR may be compressed. 
                  The contents of each file in DIR is sent to STDOUT.
                  (the '-d' flag is optional if DIR is the last argument) 
 -incl <string> : Include only file names containing string.
                  string can be a 'regexp' (case-sensitive).
 -excl <string> : Exclude file names contain string.
                  string can be a 'regexp' (case-sensitive).
 -tar           : Print files contained in gzipped .tar archives.
 -done DIR      : Directory for where to move Blast reports that 
                  have been processed. mkdir is called if necessary.
 -mon           : Print name of each file processed to STDERR.
 -h             : Print help.

QQ_USAGE_QQ
}

if($opt_mon) {
    print STDERR '='x60, "\n$ID\n\n";
    $opt_incl and print STDERR "Including only files with names passing: $opt_incl\n";
    $opt_excl and print STDERR "Excluding files with names passing: $opt_excl\n";
    $opt_done and print STDERR "Moving processed files to $opt_done\n";
    $opt_tar and print STDERR "Accessing any tar archives\n";
    print STDERR "\n";
}

# Process the supplied directory.
opendir (DIR, $opt_d) or die "Can't opendir $opt_d: $!\n\n";
&print_files($opt_d);
closedir DIR;

# Process any archives contained in the supplied directory.
if($opt_tar and @archives) {
    $opt_tar = 0;  # Not recursively exploring archives.
    $opt_mon and print STDERR "\nPrinting files in archives:\n";
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
	&print_files($newdir);
	closedir DIR;
	# Removed the opened archive.
	system("rm -rf $newdir &")== 0 or die "Can't remove opened archive dir $newdir: $!\n";
    }
}

exit 0;

#----------------
sub print_files {
#----------------
    my $d = shift;
    
    $opt_mon and print STDERR "\nPrinting files in directory $d\n";

    my ($f, $full);

    foreach $f (readdir DIR) {
	$full = "$d/$f";
	next if ($full =~ /\/\..*$/ or -d $full or not -s $full);
	# Extra test to process only certain files:
	next if $opt_incl and $f !~ /$opt_incl/o;
	next if $opt_excl and $f =~ /$opt_excl/o;
	
	$opt_mon and print STDERR "$f\n";
	
	if(-B $full) {
	    if($full =~ /\.tar/) {
		push @archives, $full if $opt_tar;
		next;
	    } else {
		print `gzcat $full`;
	    }
	} else {
	    print `cat $full`;
	}
	next unless $opt_done;
	system('mkdir', $opt_done) if not -d $opt_done;
	system('mv', '-f', $full, $opt_done) == 0 or
	    die "Can't move file to $opt_done: $!\n";
    }
}





