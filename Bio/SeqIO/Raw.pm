## $Id$
##
# BioPerl module for Bio::SeqIO::Raw
#
# Ewan Birney <birney@sanger.ac.uk> developed the SeqIO 
# schema and the first prototype modules.
#
# This code is based on his Bio::SeqIO::Fasta module with
# the necessary minor tweaks necessary to get it to read
# and write raw formatted sequences made by
# chris dagdigian <dag@sonsorol.org>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::Raw - Raw sequence input/output stream

=head1 SYNOPSIS

It is probably best not to use this object directly, but
rather go through the SeqIO handler system. Go:

    $stream = Bio::SeqIO->new(-file => $filename, -format => 'Raw');

    while $seq ( <$stream> ) {
	# $seq is a Bio::Seq object
    }


=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from raw textfiles. 
"Raw" is defined to mean 1 or more sequences separated by a newline. Raw
sequences are just lines of data with no formatting or special fields.

The write_seq() method has not been tested yet, should be simple enough.

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

=head1 AUTHOR - Ewan Birney

Email birney@sanger.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqIO::Raw;
use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::Seq;
# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use FileHandle;

@ISA = qw(Bio::Root::Object Exporter);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;
   
  my ($file,$fh) = $self->_rearrange([qw(
					 FILE
					 FH
					 )],
				     @args,
				     );
  if( $file && $fh ) {
      $self->throw("Providing both a file and a filehandle for reading from - only one please!");
  }

  if( $file ) {
      $fh = new FileHandle;
      $fh->open($file) || $self->throw("Could not open $file for Raw stream reading $!");
   }

   $self->_filehandle($fh);
      

# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq object
 Args    :


=cut

sub next_seq{
   my ($self,@args) = @_;
   my ($seq,$fh,$c,$line,$name,$desc,$sequence);

   $fh = $self->_filehandle();

   if( $fh->eof ) {
       return undef; # no throws - end of file
   }

   ## When its 1 sequence per line with no formatting at all,
   ## grabbing it should be easy :)

   $_ = <$fh>;
   $_ = uc($_);
   s/\W//g;
   $sequence = $_;

   $seq = Bio::Seq->new(-seq => $sequence);
   
   return $seq;

}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq($seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq object


=cut

sub write_seq {
   my ($self,$seq) = @_;
   my $fh = $self->_filehandle();
   my $str = $seq->seq;
   
   print $fh $str, "\n";
   return 1;
}

=head2 _filehandle

 Title   : _filehandle
 Usage   : $obj->_filehandle($newval)
 Function: 
 Example : 
 Returns : value of _filehandle
 Args    : newvalue (optional)


=cut

sub _filehandle{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_filehandle'} = $value;
    }
    return $obj->{'_filehandle'};

}

sub DESTROY {
    my $self = shift;
    my $fh;
    $fh = $self->_filehandle();

    if( defined $fh ) {
	$fh->close();
    }

    $self->{'_filehandle'} = '';
}
    






