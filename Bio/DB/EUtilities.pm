# $Id$
#
# BioPerl module for Bio::DB::EUtilities
#
# Cared for by Chris Fields <cjfields at uiuc dot edu>
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
# 
# Interfaces with new GenericWebDBI interface 

=head1 NAME

Bio::DB::EUtilities - interface for handling web queries and data
retrieval from NCBI's Entrez Utilities.

=head1 SYNOPSIS

use Bio::DB::EUtilities;

  my $esearch = Bio::DB::EUtilities->new(-eutil      => 'esearch',
                                         -db         => 'pubmed',
                                         -term       => 'hutP',
                                         -usehistory => 'y');
  
  $esearch->get_response; # parse the response, fetch a cookie
  
  my $elink = Bio::DB::EUtilities->new(-eutil        => 'elink',
                                       -db           => 'protein',
                                       -dbfrom       => 'pubmed',
                                       -cookie       => $esearch->next_cookie,
                                       -cmd          => 'neighbor_history');
  
  $elink->get_response; # parse the response, fetch the next cookie
  
  my $efetch = Bio::DB::EUtilities->new(-db           => 'nucleotide',
                                        -cookie       => $elink->next_cookie,
                                        -retmax       => 10,
                                        -rettype      => 'fasta');
  
  print $efetch->get_response->content;

=head1 DESCRIPTION

WARNING: Please do B<NOT> spam the Entrez web server with multiple requests.
NCBI offers Batch Entrez for this purpose, now accessible here via epost!

This is a test interface to NCBI's Entrez Utilities.  The main purpose of this
is to enable access to all of NCBI's databases available through Entrez and
allow for more complex queries.  It is likely that the API for this module as
well as the documentation will change dramatically over time, so novice
users and neophytes beware!

The experimental base class is Bio::DB::GenericWebDBI, which as the name
implies enables access to any web database which will accept parameters.  This
was originally born from an idea to replace WebDBSeqI/NCBIHelper with a more
general web database accession tool so one could access sequence information,
taxonomy, SNP, PubMed, and so on.  However, this may prove to be better used
as a replacement for LWP::UserAgent when accessing NCBI-related web tools
(Entrez Utilitites, or EUtilities).  Using the base class GenericWebDBI, one
could also build web interfaces to other databases to access anything via
CGI parameters.

Currently, you can access any database available through the NCBI interface:

  http://eutils.ncbi.nlm.nih.gov/

At this point, Bio::DB::EUtilities uses the EUtilities plugin modules somewhat
like Bio::SeqIO.  So, one would call the particular EUtility (epost, efetch,
and so forth) upon instantiating the object:

  my $esearch = Bio::DB::EUtilities->new(-eutil      => 'esearch',
                                         -db         => 'pubmed',
                                         -term       => 'dihydroorotase',
                                         -usehistory => 'y');

The default EUtility (when -eutil is left out) is 'efetch'.  For specifics on
each EUtility, see their respective POD (**these are incomplete**) or
the NCBI Entrez Utilities page:

  http://eutils.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html

At this time, retrieving the response is accomplished by using the method
get_response (which also parses for cookie information, see below).  This
method returns an HTTP::Response object.  The raw data is accessed by using
the object method content, like so:

  my $efetch = Bio::DB::EUtilities->new(-db           => 'nucleotide',
                                        -cookie       => $elink->next_cookie,
                                        -retmax       => 10,
                                        -rettype      => 'fasta');
  
  print $efetch->get_response->content;

Based on this, if one wanted to retrieve sequences or other raw data but did not
want to directly use Bio* objects (such as if genome sequences were to be
retrieved) one could do so by using the proper EUtility object(s) and query(ies)
and get the raw response back from NCBI through 'efetch'.  

B<A Note on Cookies:>

Certain EUtilities (epost, esearch, or elink) are able to retain information on
the NCBI server under certain settings.  This information can be retrieved by
using a 'cookie.'  Here, the idea of the 'cookie' is similar to the 'cookie' set
on a user's computer when browsing the Web.  XML data returned by these
EUtilities, when applicable, is parsed for the cookie information (the 'WebEnv'
and 'query_key' tags to be specific)  The information along with other identifying
data, such as the calling eutility, description of query, etc.) is stored as a
Bio::DB::EUtilities::Cookie object in an internal queue.  These can be retrieved
one at a time by using the next_cookie method or all at once in an array using
get_all_cookies.  Each cookie can then be 'fed', one at a time, to another
EUtility object, thus enabling chained queries as demonstrated in the synopsis.

By default, a EUtilities object will retrieve records using a cookie if the
cookie parameter is set:

my $efetch = Bio::DB::EUtilities->new(-db           => 'taxonomy',
                                      -cookie       => $elink->next_cookie);

ELink, in particular, is capable of returning multiple cookies based on the
setting for the database; if db is set to 'all', you will retrieve a cookie for
each database with related records.  

=head1 TODO

At this time no additional parsing of the returned response is enabled, but it
is anticipated that parsing of XML for ID's and other commonly requested 
information will be added in the very near future.  Resetting internal parameters
is also planned so one could feasibly reuse the objects once instantiated, such
as if one were to use this as a replacement for LWP::UserAgent when retrieving
responses i.e. many of the Bio::DB NCBI-related modules.

Any feedback is welcome.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the 
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@lists.open-bio.org               - General discussion
  http://www.bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via the web.

  http://bugzilla.bioperl.org/

=head1 AUTHOR 

Email cjfields at uiuc dot edu

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::EUtilities;
use strict;

use vars qw(@ISA $HOSTBASE %CGILOCATION $MAX_ENTRIES %DATABASE @PARAMS
            $DEFAULT_TOOL);
use Bio::DB::GenericWebDBI;
use HTTP::Request::Common;
use URI;

@ISA = qw(Bio::DB::GenericWebDBI);

BEGIN {
    $MAX_ENTRIES = 19000;
    $DEFAULT_TOOL = 'bioperl';
    # default host base
    $HOSTBASE = 'http://eutils.ncbi.nlm.nih.gov';
    # map eutility to location
    %CGILOCATION = (
            'einfo'     => ['get'  => '/entrez/eutils/einfo.fcgi', 'xml'],
            'epost'     => ['post' => '/entrez/eutils/epost.fcgi', 'xml'],
            'efetch'    => ['get'  => '/entrez/eutils/efetch.fcgi', 'dbspec'],
            'esearch'   => ['get'  => '/entrez/eutils/esearch.fcgi', 'xml'],
            'esummary'  => ['get'  => '/entrez/eutils/esummary.fcgi', 'xml'],
            'elink'     => ['get'  => '/entrez/eutils/elink.fcgi', 'xml'],
            'egquery'   => ['get'  => '/entrez/eutils/egquery.fcgi', 'xml']
             );
    # map database to return mode
    %DATABASE = ('pubmed'           => 'xml',
                 'protein'          => 'text',
                 'nucleotide'       => 'text',
                 'nuccore'          => 'text',
                 'nucgss'           => 'text',
                 'nucest'           => 'text',
                 'structure'        => 'text',
                 'genome'           => 'text',
                 'books'            => 'xml',
                 'cancerchromosomes'=> 'xml',
                 'cdd'              => 'xml',
                 'domains'          => 'xml',
                 'gene'             => 'asn1',
                 'genomeprj'        => 'xml',
                 'gensat'           => 'xml',
                 'geo'              => 'xml',
                 'gds'              => 'xml',
                 'homologene'       => 'xml',
                 'journals'         => 'text',
                 'mesh'             => 'xml',
                 'ncbisearch'       => 'xml',
                 'nlmcatalog'       => 'xml',
                 'omia'             => 'xml',
                 'omim'             => 'xml',
                 'pmc'              => 'xml',
                 'popset'           => 'xml',
                 'probe'            => 'xml',
                 'pcassay'          => 'xml',
                 'pccompound'       => 'xml',
                 'pcsubstance'      => 'xml',
                 'snp'              => 'xml',
                 'taxonomy'         => 'xml',
                 'unigene'          => 'xml',
                 'unists'           => 'xml',
                 );
    @PARAMS = qw(rettype usehistory term field tool reldate mindate
        maxdate datetype retstart retmax sort_results seq_start seq_stop strand
        complexity report dbfrom cmd holding version linkname);
	for my $method (@PARAMS) {
		eval <<END;
sub $method {
	my \$self = shift;
    return \$self->{'_$method'} = shift if \@_;
    return \$self->{'_$method'};
}
END
	}
}

sub new {
    my($class,@args) = @_;
    if( $class =~ /Bio::DB::EUtilities::(\S+)/ ) {
        my ($self) = $class->SUPER::new(@args);
        $self->_initialize(@args);
        return $self;
    } else { 
        my %param = @args;
        @param{ map { lc $_ } keys %param } = values %param; # lowercase keys
        my $eutil = $param{'-eutil'} || 'efetch';
        return unless ($class->_load_eutil_module($eutil));
        return "Bio::DB::EUtilities::$eutil"->new(@args);
    }
}

sub _initialize {
    my ($self, @args) = @_;
    my ( $tool, $ids, $retmode, $verbose, $cookie) =
      $self->_rearrange([qw(TOOL ID RETMODE VERBOSE COOKIE)],  @args);
        # hard code the base address
    $self->url_base_address($HOSTBASE);
    $tool ||= $DEFAULT_TOOL;
    $self->tool($tool);
    $ids            && $self->id($ids);
    $verbose        && $self->verbose($verbose);
    $retmode        && $self->return_mode($retmode);
    if ($cookie) {
        $self->add_cookie($cookie);
    }
}

=head2 delay_policy

  Title   : delay_policy
  Usage   : $secs = $self->delay_policy
  Function: return number of seconds to delay between calls to remote db
  Returns : number of seconds to delay
  Args    : none

  NOTE: NCBI requests a delay of 3 seconds between requests.  This method
        implements that policy.

=cut

sub delay_policy {
  my $self = shift;
  return 3;
}

=head2 add_cookie

 Title   : cookie
 Usage   : $db->add_cookie($cookie)
 Function: adds an NCBI query cookie to the internal cookie queue
 Returns : none
 Args    : a Bio::DB::EUtilities::Cookie object

=cut

sub add_cookie {
    my $self = shift;
    if (@_) {
        my $cookie = shift;
        $self->throw("Expecting a Bio::DB::EUtilities::Cookie, got $cookie.")
          unless $cookie->isa("Bio::DB::EUtilities::Cookie");
        push @{$self->{'_cookie'}}, $cookie;
    }
}

=head2 next_cookie

 Title   : next_cookie
 Usage   : $cookie = $db->next_cookie
 Function: return a cookie from the internal cookie queue
 Returns : a Bio::DB::EUtilities::Cookie object
 Args    : none

=cut

sub next_cookie {
    my $self = shift;
    my $total = $#{$self->{'_cookie'} } + 1;
    $self->debug('Cookie Jar Total: '.$total."\n" );
    if ($self->{'_cookie'}) {
        return shift @{ $self->{'_cookie'} };
    } else {
        $self->warn("No cookies left in the jar!");
    }
}

=head2 reset_cookies

 Title   : reset_cookie
 Usage   : $db->reset_cookie
 Function: resets the internal cookie queue
 Returns : none
 Args    : none

=cut

sub reset_cookies {
    my $self = shift;
    $self->{'_cookie'} = [];
}

=head2 get_all_cookies

 Title   : get_all_cookies
 Usage   : @cookies = $db->get_all_cookies
 Function: retrieves all cookies from the internal cookie queue; this leaves
           the cookies in the queue intact 
 Returns : none
 Args    : none

=cut

sub get_all_cookies {
    my $self = shift;
    return @{ $self->{'_cookie'} };
}

=head2 parse_response

 Title   : parse_response
 Usage   : $db->_parse_response($content)
 Function: parse out response for cookies and other goodies
 Returns : empty
 Args    : none
 Throws  : Not implemented (implemented in plugin classes)

=cut

sub parse_response {
  my $self = shift;
  $self->throw_not_implemented;
}

=head2 get_response

 Title   : get_response
 Usage   : $db->get_response($content)
 Function: main method to retrieve data stream; parses out response for cookie
 Returns : HTTP::Response object
 Args    : optional : Bio::DB::EUtilities::cookie from a previous search
 Throws  : 'not a cookie' exception, response errors (via HTTP::Response)

=cut

sub get_response {
    my $self = shift;
    my $response = $self->_submit_request;
    if ($response->is_error) {
        $self->throw(ref($self)." Request Error:".$response->as_string);
    }
    $self->parse_response($response);  # grab cookies and what not
    return $response;
}

# experimental; so you can send or post more that one request with the same object

sub reset_parameters {
    my $self = shift;
    my @args = @_;
    for my $method (@PARAMS) {
        if ($self->$method) {
            $self->$method(undef);
        }
    }
    $self->reset_cookies;
    $self->_initialize(@args);
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
    return $self->{'_count'} = shift if @_;
    return $self->{'_count'};
}

=head2 _eutil

 Title   : _eutil
 Usage   : $db->_eutil;
 Function: sets eutils (will make private)
 Returns : eutil
 Args    : eutil

Returns the number of entries that are matched by the query.

=cut

sub _eutil   {
    my $self = shift;
    return $self->{'_eutil'} = shift if @_;
    return $self->{'_eutil'};
}

=head1 Private methods

=head2 _submit_request

 Title   : _submit_request
 Usage   : my $url = $self->get_request
 Function: builds request object based on set parameters
 Returns : HTTP::Request
 Args    : optional : Bio::DB::EUtilities cookie

=cut

sub _submit_request {
	my $self = shift;
    my %params = $self->_get_params;
    my $eutil = $self->_eutil;
    # build id list here
    if ($self->id) {
        my @ids = @{$self->id};
        $params{'id'} = join ',', @ids;
        $self->debug
    }
    my $url = URI->new($HOSTBASE . $CGILOCATION{$eutil}[1]);
    $url->query_form(%params);
    $self->debug("The web address:\n".$url->as_string."\n");
	if ($CGILOCATION{$eutil}[0] eq 'post') {    # epost request
		return $self->post($url);
    } else {                                    # other requests
		return $self->get($url);
    }
}

=head2 _get_params

 Title   : _get_params
 Usage   : my $url = $self->_get_params
 Function: builds parameter list for web request
 Returns : hash of parameter-value paris
 Args    : optional : Bio::DB::EUtilities cookie

=cut

sub _get_params {
    my $self = shift;
    my $cookie = $self->next_cookie;
    my @final;
    # add tests for WebEnv/query_key and id (don't need both)
    my %params;
    @final =  ($cookie && $cookie->isa("Bio::DB::EUtilities::Cookie")) ?
      qw(db sort_results seq_start seq_stop strand complexity rettype
        retstart retmax cmd) :
              @PARAMS;
    for my $method (@final) {
        if ($self->$method) {
            $params{$method} = $self->$method;
        }
    }
    if ($cookie) {
        my ($webenv, $qkey) = @{$cookie->cookie};
        $self->debug("WebEnv:$webenv\tQKey:$qkey\n");
        ($params{'WebEnv'}, $params{'query_key'}) = ($webenv, $qkey);
    }
    # to get around main function sort
    if ($params{'sort_results'}) {
        $params{'sort'} = $params{'sort_results'};
        delete $params{'sort_results'};
        # sort is broken with 'pub+date', interface doesn't like escaped '+'
    }
    $params{'db'} = $self->db if $self->db;
    unless ($self->return_mode) { # set by user
        my $format = $CGILOCATION{ $self->_eutil }[2];  # set by eutil 
        if ($format eq 'dbspec') {  # database-specific
            $format = $DATABASE{$self->db};
        }
        $params{'retmode'}=$format;
    }
    $self->debug("Param: $_\tValue: $params{$_}\n") for keys %params;
    return %params;
}

sub _load_eutil_module {
  my ($self,$eutil) = @_;
  my $module = "Bio::DB::EUtilities::" . $eutil;
  my $ok;
  
  eval {
      $ok = $self->_load_module($module);
  };
  if ( $@ ) {
      print STDERR <<END;
$self: $eutil cannot be found
Exception $@
For more information about the EUtilities system please see the EUtilities docs.
This includes ways of checking for formats at compile time, not run time
END
  ;
  }
  return $ok;
}

1;
__END__