

#
# BioPerl module for Bio::SeqIO::SCF
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::SCF - SCF sequence input/output stream

=head1 SYNOPSIS

It is probably best not to use this object directly, but
rather go through the SeqIO handler system. Go:

    $stream = Bio::SeqIO->new(-file => $filename, -format => 'SCF');

    while $seq ( <$stream> ) {
	# $seq is a Bio::Seq object
    }


=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from SCF flat
file databases.

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


package Bio::SeqIO::SCF;
use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::Seq;
# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;
use IO::File;

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
      $self->throw("Providing both a file and a filehandle for reading from - oly one please!");
  }

  if( !$file && !$fh ) {
      $self->throw("Neither a file (-file) nor a filehandle (-fh) provided to SCF opening");
  }


  if( $file ) {

      $fh = new IO::File;
      $fh->open($file) || $self->throw("Could not open $file for SCF stream reading $!");
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
   my ($seq, $seqc, $fh, $buffer, %names);

   $fh = $self->_filehandle();
   binmode $fh; # for the Win32/Mac crowds

   if( eof $fh ) {
       return undef; # no throws - end of file
   }

   read $fh, $buffer, 128;

   my ($scf,
       $sample,
       $sample_offset,
       $bases,
       $bases_left_clip,
       $bases_right_clip,
       $bases_offset,
       $comment_size,
       $comments_offset,
       $version,
       $sample_size,
       $code_set,
       @header_spare ) = unpack "a4 NNNNNNNN a4 NN N20", $buffer;

   read STDIN, $buffer, $bases_offset - 128; # skip over sample info

   $seqc = '';
   for (0 .. ($bases - 1) ) {
       read STDIN, $buffer, 12;
       my ($index,
	   $prob_A,
	   $prob_C,
	   $prob_G,
	   $prob_T,
	   $base,
	   @base_spare) = unpack "N C C C C a C3", $buffer;
       $seqc .= $base;
   }

   {
       local ($/) = "\n";
       while (<STDIN>) {
	   m/([^\=]+)\=(.*)/;
	   $names{$1} = $2;
       }
   }

   $seq = Bio::Seq->new(-seq => $seqc,
			-id => $names{'NAME'},
			-names => \%names,
		       );
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

   $self->throw("Sorry, Bioperl cannot yet write SCF format!");

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

1;
