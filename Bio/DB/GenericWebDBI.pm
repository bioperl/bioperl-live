# $Id$

package Bio::DB::GenericWebDBI;
use strict;
use warnings;
use vars qw(@ISA $MODVERSION %RETRIEVAL_TYPES $DEFAULT_RETRIEVAL_TYPE
         $DEFAULT_RETURN_FORMAT $LAST_INVOCATION_TIME);
use LWP::UserAgent;
use Bio::Root::Root;

@ISA = qw(Bio::Root::Root LWP::UserAgent);

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
    my $self = $class->SUPER::new(@args);
    my ($url_base, $retrieval, $retmode, $delay, $db) =
        $self->_rearrange([qw(URL_BASE RETRIEVALTYPE RETMODE DELAY DB)],
        @args);
    # from LWP::UserAgent; set agent and env proxy
    $self->agent(ref($self)."/$Bio::Root::Root::VERSION");
    $self->env_proxy;
    $db             && $self->db($db);
    $retrieval = $DEFAULT_RETRIEVAL_TYPE unless ($retrieval);
    $retrieval      && $self->retrieval_type($retrieval);
    # these will likely be overridden in base classes
    $retmode        && $self->return_mode($retmode);
    $url_base       && $self->url_base_address($url_base);
    $delay = $self->delay_policy unless defined $delay;
    $self->delay_policy($delay);
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
    return $self->proxy($protocol,$proxy);
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

=head2 return_mode

 Title   : return_mode
 Usage   : $agent->return_mode($mode)
 Function: Get/Set return mode for query (text, xml, html, asn.1, etc)
 Returns : string for return mode
 Args    : optional string 

=cut

sub return_mode {
	my $self = shift;
	return $self->{'_retmode'} = shift if @_;
	return $self->{'_retmode'};
}

=head2 get_request

 Title   : get_request
 Usage   : $agent->get_request;
 Function: get the request based on set object parameters
 Returns : HTTP::Response object
 Args    : none
 
 This is implemented by the derived class

=cut

sub get_request {
    my ($self) = @_;
    my $msg = "Implementing class must define method get_request in class GenericWebDBI";
    $self->throw($msg);
}

=head2 request_format

 Title   : request_format
 Usage   : my ($req_format, $ioformat) = $self->request_format;
           $self->request_format("genbank");
           $self->request_format("fasta");
 Function: Get/Set sequence format retrieval. The get-form will normally not
           be used outside of this and derived modules.
 Returns : Array of two strings, the first representing the format for
           retrieval, and the second specifying the corresponding SeqIO format.
 Args    : $format = sequence format

=cut

sub request_format {
    my ($self, $value) = @_;
    if( defined $value ) {
	$self->{'_format'} = [ $value, $value];
    }
    return @{$self->{'_format'}};
}

=head2 retrieval_type

 Title   : retrieval_type
 Usage   : $self->retrieval_type($type);
           my $type = $self->retrieval_type
 Function: Get/Set a proxy for retrieval_type (pipeline, io_string or tempfile)
 Returns : string representing retrieval type
 Args    : $value - the value to store

This setting affects how the data stream from the remote web server is
processed and passed to the Bio::SeqIO layer. Three types of retrieval
types are currently allowed:

   pipeline  Perform a fork in an attempt to begin streaming
             while the data is still downloading from the remote
             server.  Disk, memory and speed efficient, but will
             not work on Windows or MacOS 9 platforms.

   io_string Store downloaded database entry(s) in memory.  Can be
             problematic for batch downloads because entire set
             of entries must fit in memory.  Alll entries must be
             downloaded before processing can begin.

   tempfile  Store downloaded database entry(s) in a temporary file.
             All entries must be downloaded before processing can
             begin.

The default is pipeline, with automatic fallback to io_string if
pipelining is not available.

=cut

sub retrieval_type {
    my ($self, $value) = @_;
    if( defined $value ) {
    $value = lc $value;
    if( ! $RETRIEVAL_TYPES{$value} ) {
        $self->warn("invalid retrieval type $value must be one of (" . 
            join(",", keys %RETRIEVAL_TYPES), ")"); 
        $value = $DEFAULT_RETRIEVAL_TYPE;
    }
    $self->{'_retrieval_type'} = $value;
    }
    return $self->{'_retrieval_type'};
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
   my $d = $self->{'_delay'};
   $self->{'_delay'} = shift if @_;
   $d;
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
      warn "sleeping for $delay seconds\n" if $self->verbose;
      sleep $delay;
   }
   $LAST_INVOCATION_TIME = time;
}

1;
__END__