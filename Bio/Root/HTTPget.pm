# $Id$
#
# BioPerl module for fallback HTTP get operations.
# Module is proxy-aware 
#
#  Cared for by Chris Dagdigian <dag@sonsorol.org>
#  but all of the good stuff was written by
#  Lincoln Stein.
# 
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Root::HTTPget - module for fallback HTTP get operations when 
LWP:: is unavailable

=head1 SYNOPSIS

 Use Bio::Root::HTTPget;

 my $response = get('http://localhost');

 $response    = get('http://localhost/images');

 $response    = eval { get('http://fred:secret@localhost/ladies_only/') 
                     } or warn $@;

 $response    = eval { get('http://jeff:secret@localhost/ladies_only/')  
                     } or warn $@;

 $response    = get('http://localhost/images/navauthors.gif');

 $response    = get(-url=>'http://www.google.com',
 		    -proxy=>'http://www.modperl.com');

=head1 DESCRIPTION

This is basically an last-chance module for doing network HTTP get requests in
situations where more advanced external CPAN modules such as LWP:: are not
installed. 

The particular reason this module was developed was so that the Open Bio
Database Access code can fallback to fetching the default registry files
from http://open-bio.org/registry/ without having to depend on
external dependencies like Bundle::LWP for network HTTP access. 

The core of this module was written by Lincoln Stein. It can handle proxies
and HTTP-based authentication.

proxy authentication.
 

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Lincoln Stein

 Cared for by Chris Dagdigian <dag@sonsorol.org>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Root::HTTPget;

use strict;
use Bio::Root::Root;
use IO::Socket qw(:DEFAULT :crlf);
use vars '@ISA';

@ISA = qw(Bio::Root::Root);


=head2 get

 Title   : get
 Usage   : 
 Function:
 Example :
 Returns : string
 Args    : 

=cut

sub get {
  my ($url,$proxy,$timeout,$auth_user,$auth_pass) = 
    __PACKAGE__->_rearrange([qw(URL PROXY TIMEOUT USER PASS)],@_);
  my $dest  = $proxy || $url;

  my ($host,$port,$path,$user,$pass) 
    = _http_parse_url($dest) or __PACKAGE__->throw("invalid URL $url");
  $auth_user ||= $user;
  $auth_pass ||= $pass;
  $path = $url if $proxy;

  # set up the connection
  my $socket = _http_connect($host,$port) or __PACKAGE__->throw("can't connect: $@");

  # the request
  print $socket "GET $path HTTP/1.0$CRLF";
  print $socket "User-Agent: Bioperl fallback fetcher/1.0$CRLF";
  # Support virtual hosts
  print $socket "HOST: $host$CRLF";

  if ($auth_user && $auth_pass) {  # authentication information
    my $token = _encode_base64("$auth_user:$auth_pass");
    print $socket "Authorization: Basic $token$CRLF";
  }
  print $socket "$CRLF";

  # read the response
  my $response;
  {
    local $/ = "$CRLF$CRLF";
    $response = <$socket>;
  }

  my ($status_line,@other_lines) = split $CRLF,$response;
  my ($stat_code,$stat_msg) = $status_line =~ m!^HTTP/1\.[01] (\d+) (.+)!
    or __PACKAGE__->throw("invalid response from web server: got $response");

  my %headers = map {/^(\S+): (.+)/} @other_lines;
  if ($stat_code == 302 || $stat_code == 301) {  # redirect
    my $location = $headers{Location} or __PACKAGE__->throw("invalid redirect: no Location header");
    return get($location,$proxy,$timeout);  # recursive call
  }

  elsif ($stat_code == 401) { # auth required
    my $auth_required = $headers{'WWW-Authenticate'};
    $auth_required =~ /^Basic realm="([^\"]+)"/
      or __PACKAGE__->throw("server requires unknown type of authentication: $auth_required");
    __PACKAGE__->throw("request failed: $status_line, realm = $1");
  }

  elsif ($stat_code != 200) {
    __PACKAGE__->throw("request failed: $status_line");
  }

  $response = '';
  while (1) {
    my $bytes = read($socket,$response,2048,length $response);
    last unless $bytes > 0;
  }

  $response;
}


=head2 _http_parse_url

 Title   :
 Usage   : 
 Function:
 Example :
 Returns :
 Args    :

=cut

sub _http_parse_url {
  my $url = shift;
  my ($user,$pass,$hostent,$path) = 
    $url =~ m!^http://(?:([^:]+):([^:]+)@)?([^/]+)(/?[^\#]*)! or return;
  $path ||= '/';
  my ($host,$port) = split(':',$hostent);
  return ($host,$port||80,$path,$user,$pass);
}

=head2 _http_connect

 Title   :
 Usage   : 
 Function:
 Example :
 Returns :
 Args    :

=cut

sub _http_connect {
  my ($host,$port,$timeout) = @_;
  my $sock = IO::Socket::INET->new(Proto     => 'tcp',
                                   Type      => SOCK_STREAM,
				   PeerHost  => $host,
				   PeerPort  => $port,
				   Timeout   => $timeout,
				  );
  $sock;
}


=head2 _encode_base64

 Title   :
 Usage   : 
 Function:
 Example :
 Returns :
 Args    :

=cut

sub _encode_base64 {
    my $res = "";
    my $eol = $_[1];
    $eol = "\n" unless defined $eol;
    pos($_[0]) = 0;                          # ensure start at the beginning

    $res = join '', map( pack('u',$_)=~ /^.(\S*)/, ($_[0]=~/(.{1,45})/gs));

    $res =~ tr|` -_|AA-Za-z0-9+/|;               # `# help emacs
    # fix padding at the end
    my $padding = (3 - length($_[0]) % 3) % 3;
    $res =~ s/.{$padding}$/'=' x $padding/e if $padding;
    # break encoded string into lines of no more than 76 characters each
    if (length $eol) {
        $res =~ s/(.{1,76})/$1$eol/g;
    }
    return $res;
}

1;
