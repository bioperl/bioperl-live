#
# BioPerl module for Bio::DB::WebQuery.pm
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

Bio::DB::Query::WebQuery - Helper class for web-based sequence queryies

=head1 SYNOPSIS

  # Do not use this class directly.  See Bio::DB::QueryI and one of
  # the implementor classes (such as Bio::DB::GenBankQuery) for
  # information.

See L<Bio::DB::QueryI>, L<Bio::DB::GenBankQuery>


=head1 DESCRIPTION

Do not use this class directly.  See Bio::DB::QueryI and one of the
implementor classes (such as Bio::DB::Query::GenBank) for information.

Those writing subclasses must define _get_params() and
_parse_response(), and possibly override _request_method().

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

package Bio::DB::Query::WebQuery;
use strict;
use URI;
use LWP::UserAgent;
use HTTP::Request::Common;

use base qw(Bio::Root::Root Bio::DB::QueryI);

=head2 new

 Title   : new
 Usage   : $db = Bio::DB::WebQuery->new(@args)
 Function: create new query object
 Returns : new query object
 Args    : -db       database (e.g. 'protein')
           -ids      array ref of ids (overrides query)
           -verbose  turn on verbose debugging

This method creates a new query object.  Typically you will specify a
-db and a -query argument.  The value of -query is a database-specific
string.

If you provide an array reference of IDs in -ids, the query will be
ignored and the list of IDs will be used when the query is passed to
the database.

=cut

# Borrowed shamelessly from WebDBSeqI.  Some of this code should be
# refactored.
sub new {
  my $class = shift;
  my $self  = $class->SUPER::new(@_);

  my ($query,$ids,$verbose) = $self->_rearrange(['QUERY','IDS','VERBOSE'],@_);
  $self->throw('must provide one of the the -query or -ids arguments')
    unless defined($query) || defined($ids);
  if ($ids) {
    $query = $self->_generate_id_string($ids);
  }
  $self->query($query);
  $verbose && $self->verbose($verbose);

  my $ua = LWP::UserAgent->new(env_proxy => 1);
  $ua->agent(ref($self) ."/".($Bio::DB::Query::WebQuery::VERSION || '0.1'));
  $self->ua($ua);
  $self->{'_authentication'} = [];
  $self;
}

=head2 ua

 Title   : ua
 Usage   : my $ua = $self->ua or 
           $self->ua($ua)
 Function: Get/Set a LWP::UserAgent for use
 Returns : reference to LWP::UserAgent Object
 Args    : $ua - must be a LWP::UserAgent

=cut

sub ua {
   my ($self, $ua) = @_;
   my $d = $self->{'_ua'};
   if( defined $ua && $ua->isa("LWP::UserAgent") ) {
      $self->{'_ua'} = $ua;
   }
   $d;
}

=head2 proxy

 Title   : proxy
 Usage   : $httpproxy = $db->proxy('http')  or 
           $db->proxy(['http','ftp'], 'http://myproxy' )
 Function: Get/Set a proxy for use of proxy
 Returns : a string indicating the proxy
 Args    : $protocol : an array ref of the protocol(s) to set/get
           $proxyurl : url of the proxy to use for the specified protocol
           $username : username (if proxy requires authentication)
           $password : password (if proxy requires authentication)

=cut

sub proxy {
    my ($self,$protocol,$proxy,$username,$password) = @_;
    return undef if ( !defined $self->ua || !defined $protocol 
		      || !defined $proxy );
    $self->authentication($username, $password) 	
	if ($username && $password);
    return $self->ua->proxy($protocol,$proxy);
}

=head2 authentication

 Title   : authentication
 Usage   : $db->authentication($user,$pass)
 Function: Get/Set authentication credentials
 Returns : Array of user/pass 
 Args    : Array or user/pass


=cut

sub authentication{
   my ($self,$u,$p) = @_;

   if( defined $u && defined $p ) {
       $self->{'_authentication'} = [ $u,$p];
   }
   return @{$self->{'_authentication'}};
}

=head2 ids

 Title   : ids
 Usage   : @ids = $db->ids([@ids])
 Function: get/set matching ids
 Returns : array of sequence ids
 Args    : (optional) array ref with new set of ids

=cut

sub ids     {
  my $self = shift;
  if (@_) {
    my $d = $self->{'_ids'};
    my $arg = shift;
    $self->{'_ids'} = ref $arg ? $arg : [$arg];
    return $d ? @$d : ();
  } else {
    $self->_fetch_ids;
    return @{$self->{'_ids'} || []};
  }
}

=head2 query

 Title   : query
 Usage   : $query = $db->query([$query])
 Function: get/set query string
 Returns : string
 Args    : (optional) new query string

=cut

sub query   {
  my $self = shift;
  my $d    = $self->{'_query'};
  $self->{'_query'} = shift if @_;
  $d;
}

=head2 _fetch_ids

 Title   : _fetch_ids
 Usage   : @ids = $db->_fetch_ids
 Function: run query, get ids
 Returns : array of sequence ids
 Args    : none

=cut

sub _fetch_ids     {
  my $self = shift;
  $self->_run_query;
  $self->_run_query(1) if $self->_truncated;
  $self->throw('Id list has been truncated even after maxids requested')
    if $self->_truncated;
  return @{$self->{'_ids'}} if $self->{'_ids'};
}

=head2 _run_query

 Title   : _run_query
 Usage   : $success = $db->_run_query
 Function: run query, parse results
 Returns : true if successful
 Args    : none

=cut

sub _run_query {
  my $self   = shift;
  my $force  = shift;

  # allow the query to be run one extra time if truncated
  return $self->{'_ran_query'} if $self->{'_ran_query'}++ && !$force;

  my $request = $self->_get_request;
  $self->debug("request is ".$request->url."\n");
  my $response = $self->ua->request($request);
  return unless $response->is_success;
  $self->debug("response is ".$response->content."\n");
  $self->_parse_response($response->content);
  1;
}

=head2 _truncated

 Title   : _truncated
 Usage   : $flag = $db->_truncated([$newflag])
 Function: get/set truncation flag
 Returns : boolean
 Args    : new flag

Some databases will truncate output unless explicitly asked
not to.  This flag allows a "two probe" attempt.

=cut

sub _truncated {
  my $self = shift;
  my $d = $self->{'_truncated'};
  $self->{'_truncated'} = shift if @_;
  $d;
}

=head2 _get_request

 Title   : _get_request
 Usage   : $http_request = $db->_get_request(@params)
 Function: create an HTTP::Request with indicated parameters
 Returns : HTTP::Request object
 Args    : CGI parameter list

=cut

sub _get_request {
  my $self   = shift;
  my ($method,$base,@params) = $self->_request_parameters;
  my $uri = URI->new($base);
  my $request;
  if ($method eq 'get') {
    $uri->query_form(@params);
    $request = GET $uri;
  } else {
    $request = POST $uri,\@params;
  }

  $request->proxy_authorization_basic($self->authentication)
	if $self->authentication;
  $request;
}

=head2 _parse_response

 Title   : _parse_response
 Usage   : $db->_parse_response($content)
 Function: parse out response
 Returns : empty
 Args    : none
 Throws  : 'unparseable output exception'

NOTE: This method must be implemented by subclass.

=cut

sub _parse_response {
  my $self    = shift;
  my $content = shift;
  $self->throw_not_implemented;
}

=head2 _request_parameters

 Title   : _request_parameters
 Usage   : ($method,$base,@params = $db->_request_parameters
 Function: return information needed to construct the request
 Returns : list of method, url base and key=>value pairs
 Args    : none

NOTE: This method must be implemented by subclass.

=cut

sub _request_parameters {
  my $self = shift;
  $self->throw_not_implemented;
}

=head2 _generate_id_string

 Title   : _generate_id_string
 Usage   : $string = $db->_generate_id_string
 Function: joins IDs together in string (implementation-dependent)
 Returns : string of concatenated IDs
 Args    : array ref of ids (normally passed into the constructor)

NOTE: This method must be implemented by subclass.

=cut

sub _generate_id_string {
  my $self = shift;
  $self->throw_not_implemented;
}

1;
