#!/usr/bin/perl -w

# Installation: 
# 1. Move this script into the same location as run_blast_remote.pl.
# 2. Edit the use lib string to be the directory containing your Bio/ directory.
use lib "/home/steve/perl/lib";

#----------------------------------------
# retrieve_blast.pl
# Retrieves a BLAST report from the NCBI queueing system.
# Author: Steve A. Chervitz, sac@neomorphic.com
# Created: Mon Oct 25 20:54:44 1999
# Revision: $Id$
#
# Input: STDIN containing the results from an NCBI BLAST request.
#        This is HTML data that contains information for
#        accessing the completed BLAST report from the NCBI
#        queueing system. The run_blast_remote.pl Bioperl script
#        now generates these files.
#
# Output: Writes two files, one containing the HTML formatted BLAST
#         report, and another containing a plain text version.
#         (Also generates some monitoring text on STDOUT).
#
# You can then feed the plain text output file into the Bioperl Blast 
# object using the parse_blast.pl script
#
#
# Example usage:
#
# 1. run_blast_remote.pl seq/yel009c.fasta -prog blastp -db swissprot            
#  
#      This generates a temp file containing HTML data from the NCBI
#      server. The file is named after the ID of the sequence in
#      the indicated fasta file. 
#
#      This temp file contains information about how 
#      to retrieve the report from the NCBI queueing system. 
#      The Blast object isn't yet savvy enough to automatically
#      retrieve the report and parse it. 
#
# 2. retrieve_blast.pl < YEL009C.blastp2.swissprot.temp.html
#       
#      This extracts the request ID from the temp.html file and 
#      requests the corresponding report from the NCBI BLAST queue.
#
#      Sample Output:
#      Obtained request ID: 940912366-18156-27559
#      GET http://www.ncbi.nlm.nih.gov/blast/blast.cgi?RID=940912366-18156-27559
#      Retrieving BLAST report from the NCBI BLAST queue...
#    
#      BLAST Search Results (786 lines, 48292 bytes)
#      Wrote HTML BLAST report to YEL009C_BLASTP.html
#      Wrote plain text BLAST report to YEL009C_BLASTP.txt
#   
# 3. parse_blast.pl 940905064-15626-17267.txt -table 1 
#
#      Parse the plain text report.
#----------------------------------------

require LWP::UserAgent;
require HTTP::Request;
require HTTP::Response;
require URI::Heuristic;

use Bio::Tools::Blast::HTML qw(&strip_html);

my $CGI = 'http://www.ncbi.nlm.nih.gov/blast/blast.cgi';

undef($/);
$data = <>;

# RID is the request ID for retrieving the BLAST report
# from the new NCBI queueing system. It Typically looks something 
# like: 940905064-15626-17267

my $rid;

if($data =~ /RID" VALUE="(\S+)"/s) {  #"
   $rid = $1;
   print "Obtained request ID: $rid\n";
}
$rid or die "Can't get RID from the input data.\n";


my $blast_url = $CGI . "?RID=$rid";

my $url = URI::Heuristic::uf_urlstr($blast_url);

$|=1;

my $ua = new LWP::UserAgent;    

my $request = new HTTP::Request(GET => $url);

print "GET $url\n", "Retrieving BLAST report from the NCBI BLAST queue...";

my $response = $ua->request($request);
my ($filestem, $content);

if($response->is_error) {
  printf "\n\nERROR: %s\n", $response->status_line();
} else {
  $content = $response->content();
  my $count;
  my $bytes;
  $bytes = length $content;
  $count = ($content =~ tr/\n/\n/);
  printf "\n\n%s (%d lines, %d bytes)\n", $response->title, $count, $bytes;

  $filestem = &get_filestem; 

  $htmlfile = "$filestem.html";
  open(OUT, ">$htmlfile") or die "Can't open output file $htmlfile: $!\n";
  print OUT $content;
  close OUT;
  print "Wrote HTML BLAST report to $htmlfile\n";

  # Now strip the HTML. Using the utility function in Bio::Tools::Blast::HTML.
  if( strip_html(\$content)) {
    my $txtfile = "$filestem.txt";
    open(OUT, ">$txtfile") or die "Can't open output file $txtfile: $!\n";
    print OUT $content;
    close OUT;
    print "Wrote plain text BLAST report to $txtfile\n";
  } else {
    print "Can't strip HTML. Plain text version not saved. Sorry.\n";
  }
}

exit(0);


# Attempts to construct a reasonable filename from the HTML BLAST report.
# Looks for the query and program name. Defaults to the RID.
sub get_filestem {

  my ($query, $prog);

  if($content =~ /Query=<\S+>\s+(\S+)/) {
    $query = $1;
  } else { 
    $query = $rid;   # Use RID as default query name (not so informative)
  }
  if($content =~ /<\S+>(T?BLAST[NPX])\s/) {
    $prog = $1;
  } else { 
    $prog = 'BLAST';
  }

  return $query . "_" . $prog;
}
