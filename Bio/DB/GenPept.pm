
#
# BioPerl module for Bio::DB::GenPept
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::GenPept - Database object interface to GenPept

=head1 SYNOPSIS

    $gb = new Bio::DB::GenPept;

    $seq = $gb->get_Seq_by_id('195055'); # Unique ID

    # or ...

    $seq = $gb->get_Seq_by_acc('195055'); # Accession Number

=head1 DESCRIPTION

Allows the dynamic retrieval of Sequence objects (Bio::Seq) from the GenPept
database at NCBI, via an Entrez query.

WARNING: Please do NOT spam the Entrez web server with multiple requests.
NCBI offers Batch Entrez for this purpose.  Batch Entrez support will likely
be supported in a future version of DB::GenPept.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Aaron Mackey

Email amackey@virginia.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::GenPept;
use vars qw($AUTOLOAD @ISA @EXPORT_OK);
use strict;

# Object preamble - inherits from Bio::DB::BioSeqI

use Bio::DB::BioSeqI;

use Bio::SeqIO;
use IO::Socket;

@ISA = qw(Bio::DB::BioSeqI Exporter);
@EXPORT_OK = qw();

# new() is inherited from Bio::DB::BioSeqI

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

# set stuff in self from @args
 return $make; # success - we hope!
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

  my $entrez = "db=p&form=6&dopt=f&html=no&title=no&uid=$uid" ;

  my $stream = $self->_get_stream($entrez);
  return $stream->next_seq();
}

=head2 get_Seq_by_acc

  Title   : get_Seq_by_acc
  Usage   : $seq = $db->get_Seq_by_acc($acc);
  Function: Gets a Bio::Seq object by its accession number
  Returns : a Bio::Seq object
  Args    : $acc : the accession number of the desired sequence entry
  Note    : For GenPept, this just calls the same code for get_Seq_by_id()

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
  my $entrez = "db=p&form=6&dopt=f&html=no&title=no&uid=$uid" ;

  return $self->_get_stream($entrez);

}

=head2 get_Stream_by_acc

  Title   : get_Stream_by_acc
  Usage   : $seq = $db->get_Seq_by_acc($acc);
  Function: Gets a series of Seq objects by accession numbers
  Returns : a Bio::SeqIO stream object
  Args    : $ref : a reference to an array of accession numbers for
                   the desired sequence entries
  Note    : For GenPept, this just calls the same code for get_Stream_by_id()

=cut

sub get_Stream_by_acc {

  my $self = shift;
  my $acc = shift or $self->throw("Must supply an accesion number!\n");

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
            from which to get the list of unique id's/accession numbers.


=cut

sub get_Stream_by_batch {
   my $self = shift;
   my $ref = shift or $self->throw("Must supply an argument!\n");
   my $which = ref($ref);
   my $fh;
   my $filename;
   if ( $which eq 'ARRAY') { # $ref is an array reference
       $fh = new_tmpfile IO::File;
       for ( @{$ref} ) {
           print $fh $_ . "\n";
       }
       seek $fh, 0, 0;
       $filename = "tempfile.txt";
   } elsif ( $which eq '') { # $ref is a filename
       $fh = new IO::File $ref, "r";
       $filename = $ref;
   } elsif ( $which eq 'GLOB' or $which eq 'IO::File') { # $ref is assumed to be a filehandle
       $fh = $ref;
       $filename = "tempfile.txt";
   }

   my $wwwbuf = "DB=n&REQUEST_TYPE=LIST_OF_GIS&FORMAT=1&HTML=FALSE&SAVETO=FALSE&NOHEADER=TRUE&UID=" . join(',', grep { chomp; } <$fh> );

   my $sock = $self->_get_sock();

   select $sock;
   print "POST /cgi-bin/Entrez/qserver.cgi HTTP/1.0\015\012";
   print "Host: www.ncbi.nlm.nih.gov\015\012";
   print "User-Agent: $0::Bio::DB::GenBank\015\012";
   print "Connection: Keep-Alive\015\012";
   print "Content-type: application/x-www-form-urlencoded\015\012";
   print "Content-length: " . length($wwwbuf) . "\015\012";
   print "\015\012";
   print $wwwbuf;

   while (<$sock>) {
       if ( m,^HTTP/\d+\.\d+\s+(\d+)[^\012]\012, ) {
           my $code = $1;
           return undef unless $code =~ /^2/;
       }
       $self->throw("Entrez Error - check query sequences!\n") if m/^ERROR/i;
       last if m/Batch Entrez results/;
   }

   return Bio::SeqIO->new('-fh' => $sock, '-format' => 'Fasta');	  

}

sub _get_stream {

  my($self, $entrez) = @_;

# most of this socket stuff is borrowed heavily from LWP::Simple, by
# Gisle Aas and Martijn Koster.  They copyleft'ed it, but we should give
# them full credit for this little diddy.

  my $sock = $self->_get_sock();

  print $sock join("\015\012" =>
		   "GET /htbin-post/Entrez/query?$entrez HTTP/1.0",
		   "Host: www.ncbi.nlm.nih.gov",
		   "User-Agent: $0::Bio::DB::GenPept",
		   "", "");

  while(<$sock>) {
    if ( m,^HTTP/\d+\.\d+\s+(\d+)[^\012]\012, ) {
      my $code = $1;
      return undef unless $code =~ /^2/;
    }
    $self->throw("Entrez Error - check query sequences!\n") if m/^ERROR/i;
    last if m/^------/; # Kludgy, but it's how L. Stein does Boulder too
  }

  return Bio::SeqIO->new('-fh' => $sock, '-format' => 'Fasta');

}

sub _get_sock {
    my $self = shift;
    my $sock = IO::Socket::INET->new(PeerAddr => 'www.ncbi.nlm.nih.gov',
                                   PeerPort => 80,
                                   Proto    => 'tcp',
                                   Timeout  => 60
                                  );
    unless ($sock) {
	$@ =~ s/^.*?: //;
	$self->throw("Can't connect to GenBank ($@)\n");
    }
    $sock->autoflush(); # just for safety's sake if they have old IO::Socket

    return $sock;
}

1;
__END__
















