#
# BioPerl module for Bio::DB::WebDBSeqI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
#

=head1 NAME

Bio::DB::WebDBSeqI - Object Interface to generalize Web Databases
for retrieving sequences

=head1 SYNOPSIS

   # get a WebDBSeqI object somehow
   # assuming it is a nucleotide db
   my $seq = $db->get_Seq_by_id('ROA1_HUMAN')

=head1 DESCRIPTION

Provides core set of functionality for connecting to a web based
database for retriving sequences.

Users wishing to add another Web Based Sequence Dabatase will need to
extend this class (see L<Bio::DB::SwissProt> or L<Bio::DB::NCBIHelper> for
examples) and implement the get_request method which returns a
HTTP::Request for the specified uids (accessions, ids, etc depending
on what query types the database accepts).

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

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via the web.

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email E<lt> jason@bioperl.org E<gt>

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::WebDBSeqI;
use strict;
use vars qw($MODVERSION %RETRIEVAL_TYPES $DEFAULT_RETRIEVAL_TYPE
	    $DEFAULTFORMAT $LAST_INVOCATION_TIME @ATTRIBUTES);

use Bio::SeqIO;
use Bio::Root::IO;
use LWP::UserAgent;
use POSIX 'setsid';
use HTTP::Request::Common;
use HTTP::Response;
use File::Spec;
use IO::Pipe;
use IO::String;
use Bio::Root::Root;

use base qw(Bio::DB::RandomAccessI);

BEGIN {
	$MODVERSION = '0.8';
	%RETRIEVAL_TYPES = ('io_string' => 1,
			    'tempfile'  => 1,
			    'pipeline'  => 1,
			  );
	$DEFAULT_RETRIEVAL_TYPE = 'pipeline';
	$DEFAULTFORMAT = 'fasta';
	$LAST_INVOCATION_TIME = 0;
}

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($baseaddress, $params, $ret_type, $format,$delay,$db) =
	 $self->_rearrange([qw(BASEADDRESS PARAMS RETRIEVALTYPE FORMAT DELAY DB)],
							 @args);

    $ret_type = $DEFAULT_RETRIEVAL_TYPE unless ( $ret_type);
    $baseaddress   && $self->url_base_address($baseaddress);
    $params        && $self->url_params($params);
    $db            && $self->db($db);
    $ret_type      && $self->retrieval_type($ret_type);
    $delay          = $self->delay_policy unless defined $delay;
    $self->delay($delay);


    # insure we always have a default format set for retrieval
    # even though this will be immedietly overwritten by most sub classes
    $format = $self->default_format unless ( defined $format &&
					     $format ne '' );

    $self->request_format($format);
    my $ua = LWP::UserAgent->new(env_proxy => 1);
    $ua->agent(ref($self) ."/$MODVERSION");
    $self->ua($ua);
    $self->{'_authentication'} = [];
    return $self;
}

# from Bio::DB::RandomAccessI

=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : $seq = $db->get_Seq_by_id('ROA1_HUMAN')
 Function: Gets a Bio::Seq object by its name
 Returns : a Bio::Seq object
 Args    : the id (as a string) of a sequence
 Throws  : "id does not exist" exception


=cut

sub get_Seq_by_id {
    my ($self,$seqid) = @_;
    $self->_sleep;
    my $seqio = $self->get_Stream_by_id([$seqid]);
    $self->throw("id does not exist") if( !defined $seqio ) ;
    if ($self->can('complexity') &&  defined $self->complexity && $self->complexity==0) {
        $self->warn("When complexity is set to 0, use get_Stream_by_id\n".
                    "Returning Bio::SeqIO object");
        return $seqio;
    }
    my @seqs;
    while( my $seq = $seqio->next_seq() ) { push @seqs, $seq; }

    # Since $seqio will not be used anymore, explicitly close its filehandle
    # or it will cause trouble later on cleanup
    $seqio->close;

    $self->throw("id '$seqid' does not exist") unless @seqs;
    if( wantarray ) { return @seqs } else { return shift @seqs }
}

=head2 get_Seq_by_acc

 Title   : get_Seq_by_acc
 Usage   : $seq = $db->get_Seq_by_acc('X77802');
 Function: Gets a Bio::Seq object by accession number
 Returns : A Bio::Seq object
 Args    : accession number (as a string)
 Throws  : "acc does not exist" exception

=cut

sub get_Seq_by_acc {
   my ($self,$seqid) = @_;
   $self->_sleep;
   my $seqio = $self->get_Stream_by_acc($seqid);
   $self->throw("acc '$seqid' does not exist") if( ! defined $seqio );
    if ($self->can('complexity') &&  defined $self->complexity && $self->complexity==0) {
        $self->warn("When complexity is set to 0, use get_Stream_by_acc\n".
                    "Returning Bio::SeqIO object");
        return $seqio;
    }
   my @seqs;
   while( my $seq = $seqio->next_seq() ) { push @seqs, $seq; }
   $self->throw("acc $seqid does not exist") unless @seqs;
   if( wantarray ) { return @seqs } else { return shift @seqs }
}


=head2 get_Seq_by_gi

 Title   : get_Seq_by_gi
 Usage   : $seq = $db->get_Seq_by_gi('405830');
 Function: Gets a Bio::Seq object by gi number
 Returns : A Bio::Seq object
 Args    : gi number (as a string)
 Throws  : "gi does not exist" exception

=cut

sub get_Seq_by_gi {
   my ($self,$seqid) = @_;
    $self->_sleep;
   my $seqio = $self->get_Stream_by_gi($seqid);
   $self->throw("gi does not exist") if( !defined $seqio );
    if ($self->can('complexity') &&  defined $self->complexity && $self->complexity==0) {
        $self->warn("When complexity is set to 0, use get_Stream_by_gi\n".
                    "Returning Bio::SeqIO object");
        return $seqio;
    }
   my @seqs;
   while( my $seq = $seqio->next_seq() ) { push @seqs, $seq; }
   $self->throw("gi does not exist") unless @seqs;
   if( wantarray ) { return @seqs } else { return shift @seqs }
}

=head2 get_Seq_by_version

 Title   : get_Seq_by_version
 Usage   : $seq = $db->get_Seq_by_version('X77802.1');
 Function: Gets a Bio::Seq object by sequence version
 Returns : A Bio::Seq object
 Args    : accession.version (as a string)
 Throws  : "acc.version does not exist" exception

=cut

sub get_Seq_by_version {
   my ($self,$seqid) = @_;
    $self->_sleep;
   my $seqio = $self->get_Stream_by_version($seqid);
   $self->throw("accession.version does not exist") if( !defined $seqio );
    if ($self->can('complexity') &&  defined $self->complexity && $self->complexity==0) {
        $self->warn("When complexity is set to 0, use get_Stream_by_version\n".
                    "Returning Bio::SeqIO object");
        return $seqio;
    }
    my @seqs;
   while( my $seq = $seqio->next_seq() ) { push @seqs, $seq; }
   $self->throw("accession.version does not exist") unless @seqs;
   if( wantarray ) { return @seqs } else { return shift @seqs }
}

# implementing class must define these

=head2 get_request

 Title   : get_request
 Usage   : my $url = $self->get_request
 Function: returns a HTTP::Request object
 Returns :
 Args    : %qualifiers = a hash of qualifiers (ids, format, etc)

=cut

sub get_request {
    my ($self) = @_;
    my $msg = "Implementing class must define method get_request in class WebDBSeqI";
    $self->throw($msg);
}

# class methods

=head2 get_Stream_by_id

  Title   : get_Stream_by_id
  Usage   : $stream = $db->get_Stream_by_id( [$uid1, $uid2] );
  Function: Gets a series of Seq objects by unique identifiers
  Returns : a Bio::SeqIO stream object
  Args    : $ref : a reference to an array of unique identifiers for
                   the desired sequence entries


=cut

sub get_Stream_by_id {
    my ($self, $ids) = @_;
    my ($webfmt,$localfmt) = $self->request_format;
    return $self->get_seq_stream('-uids' => $ids, '-mode' => 'single',
				 '-format' => $webfmt);
}

*get_Stream_by_batch = sub {
  my $self = shift;
  $self->deprecated('get_Stream_by_batch() is deprecated; use get_Stream_by_id() instead');
  $self->get_Stream_by_id(@_)
};


=head2 get_Stream_by_acc

  Title   : get_Stream_by_acc
  Usage   : $seq = $db->get_Stream_by_acc([$acc1, $acc2]);
  Function: Gets a series of Seq objects by accession numbers
  Returns : a Bio::SeqIO stream object
  Args    : $ref : a reference to an array of accession numbers for
                   the desired sequence entries
  Note    : For GenBank, this just calls the same code for get_Stream_by_id()

=cut

sub get_Stream_by_acc {
    my ($self, $ids ) = @_;
    return $self->get_seq_stream('-uids' => $ids, '-mode' => 'single');
}


=head2 get_Stream_by_gi

  Title   : get_Stream_by_gi
  Usage   : $seq = $db->get_Stream_by_gi([$gi1, $gi2]);
  Function: Gets a series of Seq objects by gi numbers
  Returns : a Bio::SeqIO stream object
  Args    : $ref : a reference to an array of gi numbers for
                   the desired sequence entries
  Note    : For GenBank, this just calls the same code for get_Stream_by_id()

=cut

sub get_Stream_by_gi {
    my ($self, $ids ) = @_;
    return $self->get_seq_stream('-uids' => $ids, '-mode' => 'gi');
}

=head2 get_Stream_by_version

  Title   : get_Stream_by_version
  Usage   : $seq = $db->get_Stream_by_version([$version1, $version2]);
  Function: Gets a series of Seq objects by accession.versions
  Returns : a Bio::SeqIO stream object
  Args    : $ref : a reference to an array of accession.version strings for
                   the desired sequence entries
  Note    : For GenBank, this is implemeted in NCBIHelper

=cut

sub get_Stream_by_version {
    my ($self, $ids ) = @_;
#    $self->throw("Implementing class should define this method!");
    return $self->get_seq_stream('-uids' => $ids, '-mode' => 'version'); # how it should work
}

=head2 get_Stream_by_query

  Title   : get_Stream_by_query
  Usage   : $stream = $db->get_Stream_by_query($query);
  Function: Gets a series of Seq objects by way of a query string or oject
  Returns : a Bio::SeqIO stream object
  Args    : $query :   A string that uses the appropriate query language
            for the database or a Bio::DB::QueryI object.  It is suggested
            that you create the Bio::DB::Query object first and interrogate
            it for the entry count before you fetch a potentially large stream.

=cut

sub get_Stream_by_query {
    my ($self, $query ) = @_;
    return $self->get_seq_stream('-query' => $query, '-mode'=>'query');
}

=head2 default_format

 Title   : default_format
 Usage   : my $format = $self->default_format
 Function: Returns default sequence format for this module
 Returns : string
 Args    : none

=cut

sub default_format {
    return $DEFAULTFORMAT;
}

# sorry, but this is hacked in because of BioFetch problems...
sub db {
  my $self = shift;
  my $d    = $self->{_db};
  $self->{_db} = shift if @_;
  $d;
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

=head2 get_seq_stream

 Title   : get_seq_stream
 Usage   : my $seqio = $self->get_seq_stream(%qualifiers)
 Function: builds a url and queries a web db
 Returns : a Bio::SeqIO stream capable of producing sequence
 Args    : %qualifiers = a hash qualifiers that the implementing class
           will process to make a url suitable for web querying

=cut

sub get_seq_stream {
	my ($self, %qualifiers) = @_;
	my ($rformat, $ioformat) = $self->request_format();
	my $seen = 0;
	foreach my $key ( keys %qualifiers ) {
		if( $key =~ /format/i ) {
			$rformat = $qualifiers{$key};
			$seen = 1;
		}
	}
	$qualifiers{'-format'} = $rformat if( !$seen);
	($rformat, $ioformat) = $self->request_format($rformat);
	# These parameters are implemented for Bio::DB::GenBank objects only
	if($self->isa('Bio::DB::GenBank')) {
		$self->seq_start() &&  ($qualifiers{'-seq_start'} = $self->seq_start());
		$self->seq_stop() && ($qualifiers{'-seq_stop'} = $self->seq_stop());
		$self->strand() && ($qualifiers{'-strand'} = $self->strand());
		defined $self->complexity() && ($qualifiers{'-complexity'} = $self->complexity());
	}
	my $request = $self->get_request(%qualifiers);
	$request->proxy_authorization_basic($self->authentication)
	    if ( $self->authentication);
	$self->debug("request is ". $request->as_string(). "\n");

	# workaround for MSWin systems
	$self->retrieval_type('io_string') if $self->retrieval_type =~ /pipeline/ && $^O =~ /^MSWin/;

	if ($self->retrieval_type =~ /pipeline/) {
		# Try to create a stream using POSIX fork-and-pipe facility.
		# this is a *big* win when fetching thousands of sequences from
		# a web database because we can return the first entry while
		# transmission is still in progress.
		# Also, no need to keep sequence in memory or in a temporary file.
		# If this fails (Windows, MacOS 9), we fall back to non-pipelined access.

		# fork and pipe: _stream_request()=><STREAM>
		my ($result,$stream) = $self->_open_pipe();

		if (defined $result) {
			$DB::fork_TTY = File::Spec->devnull; # prevents complaints from debugger
			if (!$result) { # in child process
			    $self->_stream_request($request,$stream);
			    POSIX::_exit(0); #prevent END blocks from executing in this forked child
			}
			else {
				return Bio::SeqIO->new('-verbose' => $self->verbose,
						       '-format'  => $ioformat,
						       '-fh'      => $stream);
			}
		}
		else {
			$self->retrieval_type('io_string');
		}
	}

	if ($self->retrieval_type =~ /temp/i) {
		my $dir = $self->io->tempdir( CLEANUP => 1);
		my ( $fh, $tmpfile) = $self->io()->tempfile( DIR => $dir );
		close $fh;
		my $resp = $self->_request($request, $tmpfile);
		if( ! -e $tmpfile || -z $tmpfile || ! $resp->is_success() ) {
			$self->throw("WebDBSeqI Error - check query sequences!\n");
		}
		$self->postprocess_data('type' => 'file',
				        'location' => $tmpfile);
		# this may get reset when requesting batch mode
		($rformat,$ioformat) = $self->request_format();
		if( $self->verbose > 0 ) {
			open my $ERR, '<', $tmpfile or $self->throw("Could not read file '$tmpfile': $!");
			while(<$ERR>) { $self->debug($_);}
			close $ERR;
		}

		return Bio::SeqIO->new('-verbose' => $self->verbose,
									  '-format' => $ioformat,
									  '-file'   => $tmpfile);
	}

	if ($self->retrieval_type =~ /io_string/i ) {
		my $resp = $self->_request($request);
		my $content = $resp->content_ref;
		$self->debug( "content is $$content\n");
		if (!$resp->is_success() || length($$content) == 0) {
			$self->throw("WebDBSeqI Error - check query sequences!\n");
		}
		($rformat,$ioformat) = $self->request_format();
		$self->postprocess_data('type'=> 'string',
				        'location' => $content);
		$self->debug( "str is $$content\n");
		return Bio::SeqIO->new('-verbose' => $self->verbose,
				       '-format' => $ioformat,
				       '-fh'   => new IO::String($$content));
	}

	# if we got here, we don't know how to handle the retrieval type
	$self->throw("retrieval type " . $self->retrieval_type .
					 " unsupported\n");
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
    my $d = $self->{'_baseaddress'};
    $self->{'_baseaddress'} = shift if @_;
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
    return if ( !defined $self->ua || !defined $protocol
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

=head2 url_params

 Title   : url_params
 Usage   : my $params = $self->url_params or
           $self->url_params($params)
 Function: Get/Set the URL parameters for the Web Database
 Returns : url parameters for Web Database
 Args    : $params - parameters to be appended to the URL for the WebDatabase

=cut

sub url_params {
	my ($self, $value) = @_;
	if( defined $value ) {
		$self->{'_urlparams'} = $value;
	}
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
	if( defined $ua && $ua->isa("LWP::UserAgent") ) {
		$self->{'_ua'} = $ua;
	}
	return $self->{'_ua'};
}

=head2 postprocess_data

 Title   : postprocess_data
 Usage   : $self->postprocess_data ( 'type' => 'string',
				     'location' => \$datastr);
 Function: process downloaded data before loading into a Bio::SeqIO
 Returns : void
 Args    : hash with two keys - 'type' can be 'string' or 'file'
                              - 'location' either file location or string
                                           reference containing data

=cut

sub postprocess_data {
	my ( $self, %args) = @_;
	return;
}

# private methods
sub _request {
	my ($self, $url,$tmpfile) = @_;
	my ($resp);
	if( defined $tmpfile && $tmpfile ne '' ) {
		$resp =  $self->ua->request($url, $tmpfile);
	} else {
		$resp =  $self->ua->request($url);
	}

	if( $resp->is_error  ) {
		$self->throw("WebDBSeqI Request Error:\n".$resp->as_string);
	}
	return $resp;
}

#mod_perl-safe replacement for the open(BLEH,'-|') call.  if running
#under mod_perl, detects it and closes the child's STDIN and STDOUT
#handles
sub _open_pipe {
  my ($self) = @_;
  # is mod_perl running?  Which API?
  my $mp = $self->mod_perl_api;
  if($mp and ! our $loaded_apache_sp) {
    my $load_api = ($mp == 1) ? 'use Apache::SubProcess': 'use Apache2::SubProcess';
    eval $load_api;
    $@ and $self->throw("$@\n$load_api module required for running under mod_perl");
    $loaded_apache_sp = 1;
  }

  my $pipe = IO::Pipe->new();

  local $SIG{CHLD} = 'IGNORE';
  defined(my $pid = fork)
    or $self->throw("Couldn't fork: $!");

  unless($pid) {
    #CHILD
    $pipe->writer();

    #if we're running under mod_perl, clean up some things after this fork
    if ($ENV{MOD_PERL} and my $r = eval{Apache->request} ) {
      $r->cleanup_for_exec;
      #don't read or write the mod_perl parent's tied filehandles
      close STDIN; close STDOUT;
      setsid() or $self->throw('Could not detach from parent');
    }
  } else {
    #PARENT
    $pipe->reader();
  }
  return ( $pid, $pipe );
}

# send web request to specified filehandle, or stdout, for streaming purposes
sub _stream_request {
  my $self    = shift;
  my $request = shift;
  my $dest_fh = shift || \*STDOUT;

  # fork so as to pipe output of fetch process through to
  # postprocess_data method call.
  my ($child,$fetch) = $self->_open_pipe();

  if ($child) {
    #PARENT
    local ($/) = "//\n";  # assume genbank/swiss format
    $| = 1;
    my $records = 0;
    while (my $record = <$fetch>) {
      $records++;
      $self->postprocess_data('type'     => 'string',
			      'location' => \$record);
      print $dest_fh $record;
    }
    $/ = "\n"; # reset to be safe;
    close $dest_fh; #must explicitly close here, because the hard
                    #exits don't cloes them for us
  }
  else {
    #CHILD
    $| = 1;
    my $resp =  $self->ua->request($request,
				   sub { print $fetch $_[0] }
				   );
    if( $resp->is_error  ) {
      $self->throw("WebDBSeqI Request Error:\n".$resp->as_string);
    }
    close $fetch; #must explicitly close here, because the hard exists
                  #don't close them for us
    POSIX::_exit(0);
  }
}

sub io {
    my ($self,$io) = @_;

    if(defined($io) || (! exists($self->{'_io'}))) {
	$io = Bio::Root::IO->new() unless $io;
	$self->{'_io'} = $io;
    }
    return $self->{'_io'};
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
      warn "sleeping for $delay seconds\n" if $self->verbose > 0;
      sleep $delay;
   }
   $LAST_INVOCATION_TIME = time;
}

=head2 mod_perl_api

 Title   : mod_perl_api
 Usage   : $version = self->mod_perl_api
 Function: Returns API version of mod_perl being used based on set env. variables
 Returns : mod_perl API version; if mod_perl isn't loaded, returns 0
 Args    : none

=cut

sub mod_perl_api {
    my $self = shift;
    my $v = $ENV{MOD_PERL} ?
            ( exists $ENV{MOD_PERL_API_VERSION} && $ENV{MOD_PERL_API_VERSION} >= 2 ) ?
            2 :
            1
        : 0;
    return $v;
}

1;
