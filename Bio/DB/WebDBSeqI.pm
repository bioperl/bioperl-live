# $Id$
#
# BioPerl module for Bio::DB::WebDBSeqI
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
extend this class (see Bio::DB::SwissProt or Bio::DB::NCBIHelper for
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

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via email or the
web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email jason@chg.mc.duke.edu

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::WebDBSeqI;
use strict;
use vars qw(@ISA $MODVERSION %RETRIEVAL_TYPES $DEFAULT_RETRIEVAL_TYPE
	    $DEFAULTFORMAT);

use Bio::DB::RandomAccessI;
use Bio::SeqIO;
use Bio::Root::IO;
use LWP::UserAgent;
use HTTP::Request::Common;
use HTTP::Response;
use File::Spec;
use IO::String;
use Bio::Root::Root;

@ISA = qw(Bio::DB::RandomAccessI);

BEGIN {
    $MODVERSION = '0.8';
    %RETRIEVAL_TYPES = ( 'io_string' => 1,
			 'tempfile'  => 1);
    $DEFAULT_RETRIEVAL_TYPE = 'io_string';
    $DEFAULTFORMAT = 'fasta';
}

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($baseaddress, $params, $ret_type, $format) = 
	$self->_rearrange([qw(BASEADDRESS PARAMS RETRIEVALTYPE FORMAT)],
			  @args);
    
    $ret_type = $DEFAULT_RETRIEVAL_TYPE unless ( $ret_type);
    $baseaddress && $self->url_base_address($baseaddress);
    $params      && $self->url_params($params);
    $ret_type    && $self->retrieval_type($ret_type);

    # insure we always have a default format set for retrieval
    # even though this will be immedietly overwritten by most sub classes
    $format = $self->default_format unless ( defined $format && 
					     $format ne '' );

    $self->request_format($format);
    my $ua = new LWP::UserAgent;
    $ua->agent(ref($self) ."/$MODVERSION");
    $self->ua($ua);  
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
    my $seqio = $self->get_Stream_by_id([$seqid]);
    $self->throw("id does not exist") if( !defined $seqio ) ;
    return $seqio->next_seq();
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
   my $seqio = $self->get_Stream_by_acc($seqid);
   $self->throw("acc does not exist") if( !defined $seqio );
   return $seqio->next_seq();
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
   my $seqio = $self->get_Stream_by_gi($seqid);
   $self->throw("gi does not exist") if( !defined $seqio );
   return $seqio->next_seq();
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
   $self->throw("Implementing class should define this method!"); 
   my $seqio = $self->get_Stream_by_version($seqid);
   $self->throw("accession.version does not exist") if( !defined $seqio );
   return $seqio->next_seq();
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
    return $self->get_seq_stream('-uids' => $ids, '-mode' => 'single');
}

=head2 get_Stream_by_acc

  Title   : get_Stream_by_acc
  Usage   : $seq = $db->get_Seq_by_acc([$acc1, $acc2]);
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
  Usage   : $seq = $db->get_Seq_by_gi([$gi1, $gi2]);
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
  Usage   : $seq = $db->get_Seq_by_version([$version1, $version2]);
  Function: Gets a series of Seq objects by accession.versions
  Returns : a Bio::SeqIO stream object
  Args    : $ref : a reference to an array of accession.version strings for
                   the desired sequence entries
  Note    : For GenBank, this is implemeted in NCBIHelper

=cut

sub get_Stream_by_version {
    my ($self, $ids ) = @_;
    $self->throw("Implementing class should define this method!"); 
    return $self->get_seq_stream('-uids' => $ids, '-mode' => 'version'); # how it should work
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
 Usage   : my $seqio = $self->get_seq_sream(%qualifiers)
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
    
    my $request = $self->get_request(%qualifiers);
    my ($stream,$resp);
    if( $self->retrieval_type =~ /temp/i ) {
	my $dir = $self->io()->tempdir( CLEANUP => 1);
	my ( $fh, $tmpfile) = $self->io()->tempfile( DIR => $dir );
	close $fh;
	my ($resp) = $self->_request($request, $tmpfile);		
	if( ! -e $tmpfile || -z $tmpfile || ! $resp->is_success() ) {
            $self->throw("WebDBSeqI Error - check query sequences!\n");
	}
	$self->postprocess_data('type' => 'file',
				'location' => $tmpfile);	
	# this may get reset when requesting batch mode
	($rformat,$ioformat) = $self->request_format();
	if( $self->verbose > 0 ) {
	    open(ERR, "<$tmpfile");
	    while(<ERR>) { $self->debug($_);}
	} 
	$stream = new Bio::SeqIO('-verbose' => $self->verbose,
				 '-format' => $ioformat,
				 '-file'   => $tmpfile);
    } elsif( $self->retrieval_type =~ /io_string/i ) {
	my ($resp) = $self->_request($request);
        my $content = $resp->content_ref;
	$self->debug( "content is $$content\n");
	if( ! $resp->is_success() || length(${$resp->content_ref()}) == 0 ) {
	    $self->throw("WebDBSeqI Error - check query sequences!\n");	
        }  
	($rformat,$ioformat) = $self->request_format();
	$self->postprocess_data('type'=> 'string',
				'location' => $content);
        print STDERR "str is $$content\n" if ( $self->verbose > 0);
	$stream = new Bio::SeqIO('-verbose' => $self->verbose,
				 '-format' => $ioformat,
				 '-fh'   => new IO::String($$content));
    } else { 
	$self->throw("retrieval type " . $self->retrieval_type . 
		     " unsupported\n");
    }
    return $stream;
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

=cut

sub proxy {
    my ($self,$protocol,$proxy) = @_;
    return undef if ( !defined $self->ua || !defined $protocol 
		      || !defined $proxy );
    return $self->ua->proxy($protocol,$proxy);
}

=head2 retrieval_type

 Title   : retrieval_type
 Usage   : $self->retrieval_type($type);
           my $type = $self->retrieval_type
 Function: Get/Set a proxy for retrieval_type (io_string or tempfile)
 Returns : string representing retrieval type
 Args    : $value - the value to store

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
    } else { $resp =  $self->ua->request($url); } 

    if( $resp->is_error  ) {
	$self->warn($resp->as_string());
	$self->throw("WebDBSeqI Request Error\n");
    }
    return $resp;
}

sub io {
    my ($self,$io) = @_;

    if(defined($io) || (! exists($self->{'_io'}))) {
	$io = Bio::Root::IO->new() unless $io;
	$self->{'_io'} = $io;
    }
    return $self->{'_io'};
}

sub DESTROY {    
}

1;
