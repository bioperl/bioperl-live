# $Id$
#
# $Id$
#
# BioPerl module for Bio::DB::GenBank
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code
# 
# Added LWP support - Jason Stajich 11/6/2000

=head1 NAME

Bio::DB::GenBank - Database object interface to GenBank

=head1 SYNOPSIS

    $gb = new Bio::DB::GenBank;

    $seq = $gb->get_Seq_by_id('MUSIGHBA1'); # Unique ID

    # or ...

    $seq = $gb->get_Seq_by_acc('J00522'); # Accession Number

=head1 DESCRIPTION

Allows the dynamic retrieval of Sequence objects (Bio::Seq) from the GenBank
database at NCBI, via an Entrez query.

WARNING: Please do NOT spam the Entrez web server with multiple requests.
NCBI offers Batch Entrez for this purpose.  Batch Entrez support will likely
be supported in a future version of DB::GenBank.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  bioperl-guts-l@bioperl.org         - Technically-oriented discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via email or the
web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Aaron Mackey

Email amackey@virginia.edu

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::GenBank;
use strict;
use vars qw(@ISA $DEFAULTFORMAT $MODVERSION $HOSTBASE);
# Object preamble - inherits from Bio::DB::RandomAccessI

$MODVERSION = '0.8';

use Bio::DB::RandomAccessI;
use Bio::SeqIO;
use IO::String;
use IO::File;
use Bio::Root::RootI;
use LWP::UserAgent;
use HTTP::Request::Common;
use HTTP::Response;
	    
$HOSTBASE = 'http://www.ncbi.nlm.nih.gov';
@ISA = qw(Bio::Root::RootI Bio::DB::RandomAccessI);
$DEFAULTFORMAT = 'GenBank';

# the new way to make modules a little more lightweight
sub new {
  my($class,@args) = @_;

  my $self = bless {}, $class;
  my ($format) = $self->_rearrange([qw(FORMAT)],
				   @args);
  $format = $DEFAULTFORMAT unless $format;
  $self->request_format($format);
  my $ua = new LWP::UserAgent;
  $ua->agent("$class/$MODVERSION");
  $self->ua($ua);
  return $self;
}

sub ua {
    my ($self, $ua) = @_;
    if( defined $ua && $ua->isa("LWP::UserAgent") ) {
	$self->{_ua} = $ua;
    }
    return $self->{_ua};
}

=head2 proxy

 Title   : proxy
 Usage   : $httpproxy = $db->proxy('http')  or 
           $db->proxy(['http','ftp'], 'http://myproxy' )
 Function: Get/Set a proxy for use
 Returns : a string indicating the proxy
 Args    : $protocol : an array ref of the protocol(s) to set/get
           $proxyurl : url of the proxy to use for the specified protocol

=cut


sub proxy {
    my ($self,$protocol,$proxy) = @_;
    return undef if ( !defined $self->ua );
    return $self->ua->proxy($protocol,$proxy);
}

=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : $seq = $db->get_Seq_by_id($uid);
 Function: Gets a Bio::Seq object by its unique identifier/name
 Returns : a Bio::Seq object
 Args    : $uid : the id (as a string) of the desired sequence entry

=cut

sub get_Seq_by_id {

  my $self = shift;
  my $uid = shift or $self->throw("Must supply an identifier!\n");

  my ($fmt, $streamfmt) = $self->request_format();
  my $entrez = "db=n&form=6&dopt=$fmt&html=no&title=no&uid=$uid";

  my $stream = $self->_get_stream($entrez,$streamfmt);
  my $seq = $stream->next_seq();
  $self->throw("Unable to get seq for id $uid, is it really a genbank id?\n") 
      if ( !defined $seq );
  return $seq;
}

=head2 get_Seq_by_acc

  Title   : get_Seq_by_acc
  Usage   : $seq = $db->get_Seq_by_acc($acc);
  Function: Gets a Bio::Seq object by its accession number
  Returns : a Bio::Seq object
  Args    : $acc : the accession number of the desired sequence entry
  Note    : For GenBank, this just calls the same code for get_Seq_by_id()

=cut

sub get_Seq_by_acc {

  my $self = shift;
  my $acc = shift or $self->throw("Must supply an accesion number!\n");
  
  return $self->get_Seq_by_id($acc);
}

=head2 get_Stream_by_id

  Title   : get_Stream_by_id
  Usage   : $stream = $db->get_Stream_by_id( [$uid1, $uid2] );
  Function: Gets a series of Seq objects by unique identifiers
  Returns : a Bio::SeqIO stream object
  Args    : $ref : a reference to an array of unique identifiers for
                   the desired sequence entries


=cut

sub get_Stream_by_id {

  my $self = shift;
  my $id = shift or $self->throw("Must supply a unique identifier!\n");
  ref($id) eq "ARRAY" or $self->throw("Must supply an array ref!\n");

  my $uid = join(',', @{$id});
  my ($fmt, $streamfmt) = $self->request_format();
  my $entrez = "db=n&form=6&dopt=$fmt&html=no&title=no&uid=$uid" ;

  return $self->_get_stream($entrez, $streamfmt);

}

=head2 get_Stream_by_acc

  Title   : get_Stream_by_acc
  Usage   : $seq = $db->get_Seq_by_acc($acc);
  Function: Gets a series of Seq objects by accession numbers
  Returns : a Bio::SeqIO stream object
  Args    : $ref : a reference to an array of accession numbers for
                   the desired sequence entries
  Note    : For GenBank, this just calls the same code for get_Stream_by_id()

=cut

sub get_Stream_by_acc {

  my $self = shift;
  my $acc = shift or $self->throw("Must supply an accession number!\n");

  return $self->get_Seq_by_id($acc);
}

=head2 get_Stream_by_batch

  Title   : get_Stream_by_batch
  Usage   : $seq = $db->get_Stream_by_batch($ref);
  Function: Retrieves Seq objects from Entrez 'en masse', rather than one
            at a time.  For large numbers of sequences, this is far superior
            than get_Stream_by_[id/acc]().
  Example :
  Returns : a Bio::SeqIO stream object
  Args    : $ref : either an array reference, a filename, or a filehandle
            from which to get the list of unique id's/accession numbers. #'

=cut

sub get_Stream_by_batch {
   my $self = shift;
   my $ref = shift or $self->throw("Must supply an argument!\n");
   my $which = ref($ref);
   my $fh;
   #my $filename; # doesn't seem to be used 
   if ( $which eq 'ARRAY') { # $ref is an array reference
       $fh = new_tmpfile IO::File;
       for ( @{$ref} ) {
	   print $fh $_ . "\n";
       }
       seek $fh, 0, 0;
       #$filename = "tempfile.txt";
   } elsif ( $which eq '') { # $ref is a filename
       $fh = new IO::File $ref, "r";
       #$filename = $ref;
   } elsif ( $which eq 'GLOB' or $which eq 'IO::File') { # $ref is assumed to be a filehandle
       $fh = $ref;
       #$filename = "tempfile.txt";
   }
   my ($fmt, $streamfmt) = $self->request_format();
   # unfortunately it seems that we must recode the format code ...
   if($fmt eq 'g') {
       # genbank
       $fmt = "6";
   } elsif($fmt eq 'f') {
       $fmt = "1";
   } # else we leave it as it is
   my $wwwbuf = [ Connection => 'Keep-Alive',
		  DB => 'n',		  
		  REQUEST_TYPE =>'LIST_OF_GIS', FORMAT => $fmt,  
		  HTML=>'FALSE', 
		  SAVETO=>'FALSE', 
		  NOHEADER=>'TRUE',
		  UID => join(',', grep { chomp; } <$fh> ),
		  ];
   my $url = "$HOSTBASE/cgi-bin/Entrez/qserver.cgi";
   my $stream;
   eval { 
       my $resp = $self->ua->request(POST $url,$wwwbuf );
       my $content = $resp->content;

       if( ! $resp->is_success || $content =~ m/^ERROR/i ) {
	   $self->warn($resp->error_as_HTML());
	   $self->throw("Entrez Error - check query sequences!\n");
       }
       if(! defined($streamfmt)) {
	   (undef, $streamfmt) = $self->request_format();
       }
       $stream = IO::String->new($content);
   };
   if ( $@ ) {
       $self->throw($@);
   }
   return Bio::SeqIO->new('-fh' => $stream, '-format' => $streamfmt);
}

=head2 request_format

 Title   : request_format
 Usage   : $db->request_format("GenBank");
           ($fmt, $streamfmt) = $db->request_format();
 Function: Sets the format in which all sequences for following queries
           shall be requested. At present, only "GenBank" and "Fasta" are
           supported. Using "Fasta" when you are solely interested in the
           sequence itself may pay off in performance, because the
           feature table doesnt have to be downloaded.

           As a Get method, returns an array of 2 strings, the first one
           denoting the NCBI code for the format, and the second one
           denoting the SeqIO format name. Should normally only be of use
           internally (as a client, you already get a ready-to-read stream).

 Returns : An array of two strings.
 Args    : The desired format for sequence retrieval.

=cut

sub request_format {
    my ($self, $fmt) = @_;

    if(defined($fmt)) {
	$fmt = lc $fmt;
	# only two formats are supported at present (does anyone know how to
	# request other formats from NCBI?)
	if($fmt eq 'genbank') {
	    $self->{'_format'} = "g";
	    $self->{'_streamfmt'} = "genbank";
	} elsif($fmt eq 'fasta') {
	    $self->{'_format'} = "f";
	    $self->{'_streamfmt'} = "fasta";
	} else {
	    $self->throw("currently only GenBank and Fasta supported");
	}
    }
    return ($self->{'_format'}, $self->{'_streamfmt'});
}

sub _get_stream {

  my($self, $entrez, $streamfmt) = @_;
  
  my $url = "$HOSTBASE/htbin-post/Entrez/query?$entrez";
  my $resp =  $self->ua->request(GET $url);

  my $content = $resp->content;

  if( !$resp->is_success || $content =~ m/^ERROR/i ) {
      $self->warn($resp->error_as_HTML());
      $self->throw("Entrez Error - check query sequences!\n");
  }
  my @data;
  my $read =0;
  foreach my $line ( split(/\n+/,$content) ) {
      if( $line =~ /^-----/ ) { $read =1 ; next }
      next if( ! $read ); 
      push @data, $line;
  }
  $content = join("\n", @data);
  if(! defined($streamfmt)) {
      (undef, $streamfmt) = $self->request_format();
  }
  my $stream = IO::String->new($content);
  return Bio::SeqIO->new('-fh' => $stream, '-format' => $streamfmt);
}

1;
__END__
