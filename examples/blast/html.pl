#!/usr/local/bin/perl -w

#---------------------------------------------------------------------------
# PROGRAM : html.pl
# PURPOSE : To produce HTML-formated version of report using the 
#           Bio::Tools::Blast.pm module.
# AUTHOR  : Steve A. Chervitz
# CREATED : 14 Apr 1998
# REVISION: $Id$
# WEBSITE : http://bio.perl.org/Projects/Blast/
# USAGE   : html.pl -h
# EXAMPLES: html.pl -eg
#
# INSTALLATION: 
#    Set the require ".../blast_config.pl" to point to the proper location
#    of the blast_config.pl file. See blast_config.pl for additional steps.
#
# For more information see the POD for Bio::Tools::Blast.pm,
# accessible from the above website.
#
# Sample BLAST output files can be found in examples/out/
# Sample HTML-formatted BLAST output file: examples/blast.html
#
# MODIFIED:
#  sac, 16 Jun 1998: Added installation comment, require statement comments.
#---------------------------------------------------------------------------

# Using blast_config.pl in the examples/blast distribution directory:
require "blast_config.pl"; 
# Proper path to blast_config.pl after you install it in your system:
#require "/share/www-data/html/perlOOP/bioperl/bin/blast/blast_config.pl";

# Using vars from blast_config to prevent warning messages under -w.
use vars qw($ID $VERSION $DESC $opt_compress);

$ID      = 'html.pl';
$VERSION = 0.11;
$DESC    = "Demonstrates HTML-formatting of Blast reports using Bio::Tools::Blast.pm";
$opt_method = 1;

#-----------------
sub html_usage {
#-----------------
    print STDERR "$ID, v$VERSION\n$DESC.\n";
    print STDERR <<"QQ_USAGE_QQ";

Usage: $ID blast.file [-method int] > out.html
       $ID < blast.file [-method int] > out.html

 blast.file:  Complete path name of Blast file to be processed.

 OPTIONS
 ----------
 -method <int> : The method to use (1, 2, 3, or 4):
               :   1 = Generate HTML without a special header (default).
               :   2 = Generate HTML adding a special header.
               :   3 = Trap HTML, don't print.
               :   4 = Build Blast object and then generate output.
 -h            : Print help/usage.
 -examples     : Print examples.

QQ_USAGE_QQ
#'
}

#------------
sub examples {
#------------

<<"QQ_EG_QQ";
(Run these in the examples/blast/ directory of the distribution.)

  ./$ID out/blastx.2 > out.html
  ./$ID -method 2 < out/tblastn.206.out.gz > out2.html
  ./$ID -method 3 out/blastx.2 
  ./$ID -method 4 out/blastp.2.wu > out4.html 
  ./$ID out/blastp.email.html.gz    # Should produce an error.

QQ_EG_QQ
}

##### MAIN #####

&init_blast(\&html_usage, 'method=s');

my $file = $ARGV[0] || '';

$file || print STDERR "\nReading Blast report from STDIN.\n";

if($opt_method == 1) {
    ## EXAMPLE 1:
    ##
    # Uses STDIN if no file is supplied.
    print STDERR "Not providing any special header.\n";
    # Results are sent directly to STDOUT unless an -out parameter is supplied.
    # This permits faster response time. Blast reports can be huge.
    $Blast->to_html($file);  

} elsif($opt_method == 2) {

    ## EXAMPLE 2: 
    ##
    print STDERR "Providing a special header.\n";
    $Blast->to_html( -FILE  =>$file, 
		     -HEADER=>qq|<html><head><title>BLAST results</title></head><center><h1>BLAST Results</h1><a href="our_home.html">Provided by us</a></center><hr>| );

} elsif($opt_method == 3) {
    ## EXAMPLE 3: 
    ##
    print STDERR "Trapping the HTML output for further processing:\n";
    my (@out);
    $Blast->to_html( -file  =>$file, -out=> \@out );
    printf "%d lines of blast html read.\n", scalar(@out);

} else {
    ## example 4: 
    ## using 'read' mode since we don't need to parse it if we are just
    ## dumping out html.
    print STDERR "Building a Blast.pm object and then generating HTML.\n";
    my $blast = new Bio::Tools::Blast -FILE => $file, -READ => 1;
    $blast->to_html;
}

$opt_compress && $Blast->compress_file($file);

exit(0);

