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

Bio::DB::GenericWebDBI - abstract interface for parameter-based remote
database access

=head1 SYNOPSIS

  #
  # grab data from HTTP::Response object using concrete class
  #

  $data = $db->get_response->content;

  #
  # $data is the raw data output from the HTTP::Response object;
  # this data may be preparsed using the private method _parse_response

=head1 DESCRIPTION

WARNING: Please do B<NOT> spam the web servers with multiple requests.

This class acts as a user agent interface for any generic web database, but
is specifically geared towards CGI-based databases which accept parameters.

=head1 TODO

File and filehandle support to be added

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

  http://bugzilla.open-bio.org/

=head1 AUTHOR

Email cjfields at uiuc dot edu

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::GenericWebDBI;
use strict;
use warnings;
use vars qw($MODVERSION %RETRIEVAL_TYPES $DEFAULT_RETRIEVAL_TYPE
         $DEFAULT_RETURN_FORMAT $LAST_INVOCATION_TIME);

use base qw(Bio::Root::Root LWP::UserAgent);

BEGIN {
    $MODVERSION = '0.8';
    %RETRIEVAL_TYPES = ('io_string' => 1,
                'tempfile'  => 1,
                'pipeline'  => 1,
                );
    $DEFAULT_RETRIEVAL_TYPE = 'pipeline';
    $DEFAULT_RETURN_FORMAT = 'text';
    $LAST_INVOCATION_TIME = 0;
}

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args, env_proxy => 1);
    my ($url_base, $retmode, $delay, $db) =
        $self->_rearrange([qw(URL_BASE RETMODE DELAY DB)],
        @args);
    # from LWP::UserAgent; set agent and env proxy
    $self->agent(ref($self)."/$Bio::Root::Root::VERSION");;
    $db             && $self->db($db);
    # these will likely be overridden in base classes
    $retmode        && $self->retmode($retmode);
    $url_base       && $self->url_base_address($url_base);
    # delay policy needs to be worked out; not set up correctly
    $delay = defined($delay) ? $delay: $self->delay_policy;
    $self->delay($delay);
    return $self;
}

=head2 url_base_address

 Title   : url_base_address
 Usage   : my $address = $self->url_base_address or
           $self->url_base_address($address)
 Function: Get/Set the base URL for the Web Database
 Returns : Base URL for the Web Database
 Args    : $address - URL for the WebDatabase

=cut

sub url_base_address {
    my $self = shift;
    return $self->{'_baseaddress'} = shift if @_;
    return $self->{'_baseaddress'};
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
    return undef if ( !defined $protocol || !defined $proxy );
    $self->authentication($username, $password)
    if ($username && $password);
    return $self->SUPER::proxy($protocol,$proxy);
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

=head2 db

 Title   : db
 Usage   : $db->db
 Function: Get/Set database parameter
 Returns : string
 Args    : optional string

=cut

sub db {
	my $self = shift;
	return $self->{'_db'} = shift if @_;
	return $self->{'_db'};
}

=head2 id

 Title   : id
 Usage   : $agent->id($id)
           $agent->id(\@id)
 Function: Get/Set id(s)
 Returns : reference to id(s)
 Args    : a single id or reference to array of id(s)

=cut

sub id {
	my $self = shift;
    if (@_) {
        my $id = shift;
        if (ref($id) !~ /ARRAY/) { # single ID
            $self->{'_ids'} = [$id];
        }
        else {
            $self->{'_ids'} = $id;
        }
    }
	return $self->{'_ids'};
}

=head2 retmode

 Title   : retmode
 Usage   : $agent->retmode($mode)
 Function: Get/Set return mode for query (text, xml, html, asn.1, etc)
 Returns : string for return mode
 Args    : optional string

=cut

sub retmode {
	my $self = shift;
	return $self->{'_retmode'} = shift if @_;
	return $self->{'_retmode'};
}

=head2 get_response

 Title   : get_response
 Usage   : $agent->get_response;
 Function: get the request based on set object parameters, retrieved using
           the private method _get_params
 Returns : HTTP::Response object
 Args    : none

 This is implemented by the derived class

=cut

sub get_response {
    my ($self) = @_;
    my $msg = "Implementing class must define method get_response in class GenericWebDBI";
    $self->throw($msg);
}

=head2 delay

 Title   : delay
 Usage   : $secs = $self->delay([$secs])
 Function: get/set number of seconds to delay between fetches
 Returns : number of seconds to delay
 Args    : new value

NOTE: the default is to use the value specified by delay_policy().
This can be overridden by calling this method, or by passing the
-delay argument to new().

=cut

sub delay {
   my $self = shift;
   return $self->{'_delay'} = shift if @_;
   return $self->{'_delay'};
}

=head2 delay_policy

 Title   : delay_policy
 Usage   : $secs = $self->delay_policy
 Function: return number of seconds to delay between calls to remote db
 Returns : number of seconds to delay
 Args    : none

NOTE: The default delay policy is 0s.  Override in subclasses to
implement delays.  The timer has only second resolution, so the delay
will actually be +/- 1s.

=cut

sub delay_policy {
   my $self = shift;
   return 0;
}

=head2 _submit_request

  Title   : _submit_request
  Usage   : my $url = $self->get_request
  Function: builds request object based on set parameters
  Returns : HTTP::Request
  Args    : optional : Bio::DB::EUtilities cookie

=cut

sub _submit_request {
    my ($self) = @_;
    my $msg = "Implementing class must define method _submit_request in class GenericWebDBI";
    $self->throw($msg);
}

=head2 _get_params

  Title   : _get_params
  Usage   : my $url = $self->_get_params
  Function: builds parameter list for web request
  Returns : hash of parameter-value paris
  Args    : optional : Bio::DB::EUtilities cookie

=cut

# these get sorted out in a hash originally but end up in an array to
# deal with multiple id parameters (hash values would kill that)

sub _get_params {
    my ($self) = @_;
    my $msg = "Implementing class must define method _get_params in class GenericWebDBI";
    $self->throw($msg);
}

=head2 _sleep

 Title   : _sleep
 Usage   : $self->_sleep
 Function: sleep for a number of seconds indicated by the delay policy
 Returns : none
 Args    : none

NOTE: This method keeps track of the last time it was called and only
imposes a sleep if it was called more recently than the delay_policy()
allows.

=cut

sub _sleep {
   my $self = shift;
   my $last_invocation = $LAST_INVOCATION_TIME;
   if (time - $LAST_INVOCATION_TIME < $self->delay) {
      my $delay = $self->delay - (time - $LAST_INVOCATION_TIME);
      $self->debug("sleeping for $delay seconds\n");
      sleep $delay;
   }
   $LAST_INVOCATION_TIME = time;
}

1;
__END__
