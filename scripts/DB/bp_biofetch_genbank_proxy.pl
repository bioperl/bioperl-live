#!perl

# dbfetch style caching proxy for GenBank
use strict;
use warnings;
use CGI qw(:standard);
use HTTP::Request::Common;
use LWP::UserAgent;
use Cache::FileCache;

use vars qw(%GOT $BUFFER %MAPPING $CACHE);

use constant CACHE_LOCATION => '/usr/tmp/dbfetch_cache';
use constant MAX_SIZE   => 100_000_000;  # 100 megs, roughly
use constant CACHE_DEPTH => 4;
use constant EXPIRATION => "1 week";
use constant PURGE      => "1 hour";

%MAPPING = (genbank => {db=>'nucleotide',
			rettype => 'gb'},
	    genpep  => {db=>'protein',
			rettype => 'gp'});
# we're doing everything in callbacks, so initialize globals.
$BUFFER = '';
%GOT    = ();

print header('text/plain');

param() or print_usage();

my $db     = param('db');
my $style  = param('style');
my $format = param('format');
my $id     = param('id');
my @ids    = split /\s+/,$id;

$format = 'genbank' if $format eq 'default';  #h'mmmph

$MAPPING{$db}        or error(1=>"Unknown database [$db]");
$style  eq 'raw'     or error(2=>"Unknown style [$style]");
$format eq 'genbank' or error(3=>"Format [$format] not known for database [$db]");

$CACHE = Cache::FileCache->new({cache_root          => CACHE_LOCATION,
				default_expires_in  => EXPIRATION,
				cache_DEPTH         => CACHE_DEPTH,
				namespace           => 'dbfetch',
				auto_purge_interval => PURGE});

# handle cached entries
foreach (@ids) {
  if (my $obj = $CACHE->get($_)) {
    $GOT{$_}++;
    print $obj,"//\n";
  }
}

# handle the remainder
@ids = grep {!$GOT{$_}} @ids;
if (@ids) {
  my $request = POST('http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',
		     [rettype    => $MAPPING{$db}{rettype},
		      db         => $MAPPING{$db}{db},
		      tool       => 'bioperl',
		      retmode    => 'text',
		      usehistory => 'n',
		      id         => join(',',@ids),
		     ]
		    );

  my $ua = LWP::UserAgent->new;
  my $response = $ua->request($request,\&callback);

  if ($response->is_error) {
    my $status = $response->status_line;
    error(6 => "HTTP error from GenBank [$status]");
  }
}

my @missing_ids = grep {!$GOT{$_}} @ids;
foreach (@missing_ids) {
  error(4=>"ID [$_] not found in database [$db]",1);
}

# my $response = $response->content;

sub process_record {
  my $record = shift;
  print "$record//\n";
  my ($locus)       = $record =~ /^LOCUS\s+(\S+)/m;
  my ($accession)   = $record =~ /^ACCESSION\s+(\S+)/m;
  my ($version,$gi) = $record =~ /^VERSION\s+(\S+)\s+GI:(\d+)/m;
  foreach ($locus,$accession,$version,$gi) {
    $GOT{$_}++;
    $CACHE->set($_,$record);
  }
}

sub callback {
  my $data = shift;
  $BUFFER .= $data;
  my $index = 0;
  while (($index = index($BUFFER,"//\n\n",$index))>=0) {
    my $record = substr($BUFFER,0,$index);
    $index += length("//\n\n");
    substr($BUFFER,0,$index) = '';
    process_record($record);
  }
}



sub print_usage {
  print <<'END';
This script is intended to be used non-interactively.

Brief summary of arguments:
URL

This interface does not specify what happens when biofetch is called
in interactive context. The implementations can return the entries
decorated with HTML tags and hypertext links.

A URL for biofetch consists of four sections:

			e.g.
1. protocol		http://
2. host			www.ebi.ac.uk
3. path to program	/Tools/dbfetch/dbfetch
4. query string		?style=raw;format=embl;db=embl;id=J00231


QUERY STRING

The query string options are separated from the base URL (protocol +
host + path) by a question mark (?) and from each other by a semicolon
';' (or by ampersand '&'). See CGI GET documents at
http://www.w3.org/CGI/). The order of options is not critical. It is
recommended to leave the ID to be the last item.

Input for options should be case insensitive.


option: db

  Option  : db
  Descr   : database name
  Type    : required
  Usage   : db=genpep | db=genbank
  Arg     : string 

Currently this server accepts "genbank" and "genpep"

option: style

  Option  : style
  Descr   : +/- HTML tags
  Type    : required
  Usage   : style=raw | db=html
  Arg     : enum (raw|html)

In non-interactive context, always give "style=raw". This uses
"Content-Type: text/plain". If other content types are needed (XML),
this part of the spesifications can be extended to accommodate them.

This server only accepts "raw".


option: format

  Option  : format
  Descr   : format of the database entries returned
  Type    : optional
  Usage   : format=genbank
  Arg     : enum

Format defaults to the distribution format of the database (embl for
EMBL database). If some other supported format is needed this option
is needed (E.g. formats for EMBL: fasta, bsml, agave).

This server only accepts "genbank" format.

option: id

  Option  : id
  Descr   : unique database identifier(s)
  Type    : required
  Usage   : db=J00231 | id=J00231+BUM
  Arg     : string 

The ID option should be able to process all UIDS in a database. It
should not be necessary to know if the UID is an ID, accession number
or accession.version.

The number of entry UIDs allowed is implementation specific. If the
limit is exceeded, the the program reports an error. The UIDs should
be separated by spaces (use '+' in a GET method string).


ERROR MESSAGES

The following standardized one line messages should be printed out in
case of an error.

ERROR 1 Unknown database [$db].
ERROR 2 Unknown style [$style].
ERROR 3 Format [$format] not known for database [$db].
ERROR 4 ID [$id] not found in database [$db].
ERROR 5 Too many IDs [$count]. Max [$MAXIDS] allowed.

END
;

exit 0;
}

sub error {
  my ($code,$message,$noexit) = @_;
  print "ERROR $code $message\n";
  exit 0 unless $noexit;
}

__END__

=head1 NAME

bp_biofetch_genbank_proxy.pl - Caching BioFetch-compatible web proxy for GenBank

=head1 SYNOPSIS

  Install in cgi-bin directory of a Web server.  Stand back.

=head1 DESCRIPTION

This CGI script acts as the server side of the BioFetch protocol as
described in http://obda.open-bio.org/Specs/.  It provides two
database access services, one for data source "genbank" (nucleotide
entries) and the other for data source "genpep" (protein entries).

This script works by forwarding its requests to NCBI's eutils script,
which lives at http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi.
It then reformats the output according to the BioFetch format so the
sequences can be processed and returned by the Bio::DB::BioFetch
module.  Returned entries are temporarily cached on the Web server's
file system, allowing frequently-accessed entries to be retrieved
without another round trip to NCBI.

=head2 INSTALLATION

You must have the following installed in order to run this script:

   1) perl
   2) the perl modules LWP and Cache::FileCache
   3) a web server (Apache recommended)

To install this script, copy it into the web server's cgi-bin
directory.  You might want to shorten its name; "dbfetch" is
recommended.

There are several constants located at the top of the script that you
may want to adjust.  These are:

CACHE_LOCATION

This is the location on the filesystem where the cached files will be
located.  The default is /usr/tmp/dbfetch_cache.

MAX_SIZE

This is the maximum size that the cache can grow to.  When the cache
exceeds this size older entries will be deleted automatically.  The
default setting is 100,000,000 bytes (100 MB).

EXPIRATION

Entries that haven't been accessed in this length of time will be
removed from the cache.  The default is 1 week.

PURGE

This constant specifies how often the cache will be purged for older
entries.  The default is 1 hour.

=head1 TESTING

To see if this script is performing as expected, you may test it with
this script:

 use Bio::DB::BioFetch;
 my $db = Bio::DB::BioFetch->new(-baseaddress=>'http://localhost/cgi-bin/dbfetch',
	 			 -format     =>'genbank',
				 -db         =>'genbank');
 my $seq = $db->get_Seq_by_id('DDU63596');
 print $seq->seq,"\n";

This should print out a DNA sequence.

=head1 SEE ALSO

L<Bio::DB::BioFetch>, L<Bio::DB::Registry>

=head1 AUTHOR

Lincoln Stein, E<lt>lstein-at-cshl.orgE<gt>

Copyright (c) 2003 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

