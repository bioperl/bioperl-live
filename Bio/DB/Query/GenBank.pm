#
# BioPerl module for Bio::DB::Query::GenBank.pm
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Lincoln Stein <lstein@cshl.org>
#
# Copyright Lincoln Stein
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
#

=head1 NAME

Bio::DB::Query::GenBank - Build a GenBank Entrez Query

=head1 SYNOPSIS

   use Bio::DB::Query::GenBank;
   use Bio::DB::GenBank;

   my $query_string = 'Oryza[Organism] AND EST[Keyword]';
   my $query = Bio::DB::Query::GenBank->new(-db => 'nucleotide',
                                            -query => $query_string,
                                            -mindate => '2001',
                                            -maxdate => '2002');

   print $query->count,"\n";

   # get a Genbank database handle
   my $gb = Bio::DB::GenBank->new();
   my $stream = $gb->get_Stream_by_query($query);
   while (my $seq = $stream->next_seq) {
      # do something with the sequence object
   }

   # initialize the list yourself
   my $query = Bio::DB::Query::GenBank->new(-ids=>[195052,2981014,11127914]);


=head1 DESCRIPTION

This class encapsulates NCBI Entrez queries.  It can be used to store
a list of GI numbers, to translate an Entrez query expression into a
list of GI numbers, or to count the number of terms that would be
returned by a query.  Once created, the query object can be passed to
a Bio::DB::GenBank object in order to retrieve the entries
corresponding to the query.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Lincoln Stein

Email lstein@cshl.org

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::Query::GenBank;
use strict;
use URI::Escape 'uri_unescape';
use Bio::DB::NCBIHelper;


#use constant EPOST       => $Bio::DB::NCBIHelper::HOSTBASE . '/entrez/eutils/epost.fcgi';
#use constant ESEARCH     => $Bio::DB::NCBIHelper::HOSTBASE . '/entrez/eutils/esearch.fcgi';
# the reference to the our variable of the $Bio::DB::NCBIHelper::HOSTBASE doesn't seem to work in 
# the constant definition in perl 5.10.1 or 5.16.3
use constant EPOST       => '/entrez/eutils/epost.fcgi';
use constant ESEARCH     => '/entrez/eutils/esearch.fcgi';
use constant DEFAULT_DB  => 'protein';
use constant MAXENTRY    => 100;

use vars qw(@ATTRIBUTES);

use base qw(Bio::DB::Query::WebQuery);

BEGIN {
  @ATTRIBUTES = qw(db reldate mindate maxdate datetype maxids);
  for my $method (@ATTRIBUTES) {
    eval <<END;
sub $method {
   my \$self = shift;
   my \$d    = \$self->{'_$method'};
   \$self->{'_$method'} = shift if \@_;
   \$d;
}
END
  }
}

=head2 new

 Title   : new
 Usage   : $db = Bio::DB::Query::GenBank->new(@args)
 Function: create new query object
 Returns : new query object
 Args    : -db       database (see below for allowable values)
           -query    query string
           -mindate  minimum date to retrieve from (YYYY/MM/DD)
           -maxdate  maximum date to retrieve from (YYYY/MM/DD)
           -reldate  relative date to retrieve from (days)
           -datetype date field to use ('edat' or 'mdat')
           -ids      array ref of gids (overrides query)
           -maxids   the maximum number of IDs you wish to collect
                     (defaults to 100)

This method creates a new query object.  Typically you will specify a
-db and a -query argument, possibly modified by -mindate, -maxdate, or
-reldate.  -mindate and -maxdate specify minimum and maximum dates for
entries you are interested in retrieving, expressed in the form
YYYY/MM/DD.  -reldate is used to fetch entries that are more recent
than the indicated number of days.

If you provide an array reference of IDs in -ids, the query will be
ignored and the list of IDs will be used when the query is passed to a
Bio::DB::GenBank object's get_Stream_by_query() method.  A variety of
IDs are automatically recognized, including GI numbers, Accession
numbers, Accession.version numbers and locus names.

By default, the query will collect only the first 100 IDs and will
generate an exception if you call the ids() method and the query
returned more than that number.  To increase this maximum, set -maxids
to a number larger than the number of IDs you expect to obtain.  This
only affects the list of IDs you obtain when you call the ids()
method, and does not affect in any way the number of entries you
receive when you generate a SeqIO stream from the query.

-db option values:

  The most commonly used databases are:

      protein
      nucleotide
      nuccore
      nucgss
      nucest
      unigene

  An up to date list of database names supported by NCBI eUtils is
  always available at:
  http://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?

  However, note that not all of these databases return datatypes that
  are parsable by Bio::DB::GenBank

=cut

sub new {
  my $class = shift;
  my $self  = $class->SUPER::new(@_);
  my ($query,$db,$reldate,$mindate,$maxdate,$datetype,$ids,$maxids)
    = $self->_rearrange([qw(QUERY DB RELDATE MINDATE MAXDATE DATETYPE IDS MAXIDS)],@_);
  $self->db($db || DEFAULT_DB);
  $reldate  && $self->reldate($reldate);
  $mindate  && $self->mindate($mindate);
  $maxdate  && $self->maxdate($maxdate);
  $maxids   && $self->maxids($maxids);
  $datetype ||= 'mdat';
  $datetype && $self->datetype($datetype);
  $self;
}

=head2 cookie

 Title   : cookie
 Usage   : ($cookie,$querynum) = $db->cookie
 Function: return the NCBI query cookie
 Returns : list of (cookie,querynum)
 Args    : none

NOTE: this information is used by Bio::DB::GenBank in
conjunction with efetch.

=cut

sub cookie {
  my $self = shift;
  if (@_) {
    $self->{'_cookie'}   = shift;
    $self->{'_querynum'} = shift;
  }

  else {
    $self->_run_query;
    @{$self}{qw(_cookie _querynum)};
  }
}

=head2 _request_parameters

 Title   : _request_parameters
 Usage   : ($method,$base,@params = $db->_request_parameters
 Function: return information needed to construct the request
 Returns : list of method, url base and key=>value pairs
 Args    : none

=cut

sub _request_parameters {
  my $self = shift;
  my ($method,$base);
  my @params = map {eval("\$self->$_") ? ($_ => eval("\$self->$_")) : () } @ATTRIBUTES;
  push @params,('usehistory'=>'y','tool'=>'bioperl');
  $method = 'get';
  
  $base   = $Bio::DB::NCBIHelper::HOSTBASE.ESEARCH; # this seems to need to be dynamic
  push @params,('term'   => $self->query);
  # Providing 'retmax' limits queries to 500 sequences  ?? I don't think so LS
  push @params,('retmax' => $self->maxids || MAXENTRY);

  # And actually, it seems that we need 'retstart' equal to 0 ?? I don't think so LS
  # push @params, ('retstart' => 0);

  ($method,$base,@params);
}


=head2 count

 Title   : count
 Usage   : $count = $db->count;
 Function: return count of number of entries retrieved by query
 Returns : integer
 Args    : none

Returns the number of entries that are matched by the query.

=cut

sub count   {
  my $self = shift;
  if (@_) {
    my $d = $self->{'_count'};
    $self->{'_count'}   = shift;
    return $d;
  }
  else {
    $self->_run_query;
    return $self->{'_count'};
  }
}

=head2 ids

 Title   : ids
 Usage   : @ids = $db->ids([@ids])
 Function: get/set matching ids
 Returns : array of sequence ids
 Args    : (optional) array ref with new set of ids

=cut

=head2 query

 Title   : query
 Usage   : $query = $db->query([$query])
 Function: get/set query string
 Returns : string
 Args    : (optional) new query string

=cut

=head2 _parse_response

 Title   : _parse_response
 Usage   : $db->_parse_response($content)
 Function: parse out response
 Returns : empty
 Args    : none
 Throws  : 'unparseable output exception'

=cut

sub _parse_response {
  my $self    = shift;
  my $content = shift;
  if (my ($warning) = $content =~ m!<ErrorList>(.+)</ErrorList>!s) {
    $self->warn("Warning(s) from GenBank: $warning\n");
  }
  if (my ($error) = $content =~ /<OutputMessage>([^<]+)/) {
    $self->throw("Error from Genbank: $error");
  }

  my ($count) = $content =~  /<Count>(\d+)/;
  my ($max)   = $content =~  /<RetMax>(\d+)/;
  my $truncated = $count > $max;
  $self->count($count);
  if (!$truncated) {
    my @ids = $content =~ /<Id>(\d+)/g;
    $self->ids(\@ids);
  } else {
    $self->debug("ids truncated at $max\n");
  }
  $self->_truncated($truncated);
  my ($cookie)    = $content =~ m!<WebEnv>(\S+)</WebEnv>!;
  my ($querykey)  = $content =~ m!<QueryKey>(\d+)!;
  $self->cookie(uri_unescape($cookie),$querykey);
}

=head2 _generate_id_string

 Title   : _generate_id_string
 Usage   : $string = $db->_generate_id_string
 Function: joins IDs together in string (possibly implementation-dependent)
 Returns : string of concatenated IDs
 Args    : array ref of ids (normally passed into the constructor)

=cut

sub _generate_id_string {
    my ($self, $ids) = @_;
    # this attempts to separate out accs (alphanumeric) from UIDs (numeric only)
    # recent changes to esearch has wrought this upon us.. cjf 4/19/07
    return sprintf('%s',join('|',map {
      ($_ =~ m{^\d+$}) ? $_.'[UID]' : $_.'[PACC]'
    } @$ids));
}

1;
