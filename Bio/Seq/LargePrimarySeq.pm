
#
# BioPerl module for Bio::LargePrimarySeq
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::LargePrimarySeq - PrimarySeq object that stores sequence as files in /tmp 

=head1 SYNOPSIS

  # normal primary seq usage

=head1 DESCRIPTION

This object stores a sequence as a series of files in a temporary
directory. The aim is to allow someone the ability to store very large
sequences (eg, > 100MBases) in a file system without running out of memory
(eg, on a 64 MB real memory machine!). 

Of course, to actually make use of this functionality, the programs
which use this object B<must> not call $primary_seq->seq otherwise the
entire sequence will come out into memory and probably paste your
machine. However, calls $primary_seq->subseq(10,100) will cause only
90 characters to be brought into real memory.

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

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::LargePrimarySeq;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::RootI;
use Bio::PrimarySeqI;

use POSIX;
use FileHandle;

my $DEFAULT_TEMP_DIR = "/tmp";
my $counter = 0;

@ISA = qw(Bio::PrimarySeqI Bio::Root::RootI);

sub new {
    my ($class) = @_;

    my $self = {};
    bless $self,$class;
   
    my $id = POSIX::cuserid() . "$$" . $counter;
    my $file = "$DEFAULT_TEMP_DIR/$id";
    system("touch $file");
    my $fh = new FileHandle $file,"+>"; 
    if( !defined $fh ) {
	$self->thro("Unable to make file $file $!");
    }

    $self->_fh($fh);
    $self->_filename($file);
    $self->length(0);

    return $self;
}

=head2 id

 Title   : id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub id{
   my ($self,@args) = @_;

   return $self->display_id(@args);
}

=head2 display_id

 Title   : display_id
 Usage   : $obj->display_id($newval)
 Function: 
 Example : 
 Returns : value of display_id
 Args    : newvalue (optional)


=cut

sub display_id{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'display_id'} = $value;
    }
    return $obj->{'display_id'};

}

=head2 accession_number

 Title   : accession_number
 Usage   : $obj->accession_number($newval)
 Function: 
 Example : 
 Returns : value of accession_number
 Args    : newvalue (optional)


=cut

sub accession_number{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'accession_number'} = $value;
    }
    return $obj->{'accession_number'};

}

=head2 length

 Title   : length
 Usage   : $obj->length($newval)
 Function: 
 Example : 
 Returns : value of length
 Args    : newvalue (optional)


=cut

sub length{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'length'} = $value;
    }
    return $obj->{'length'};

}

=head2 seq

 Title   : seq
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub seq{
   my ($self) = @_;

   return $self->subseq(1,$self->length);
}

=head2 subseq

 Title   : subseq
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub subseq{
   my ($self,$start,$end) = @_;

   if( $start < 1 || $end > $self->length ) {
       $self->throw("Attempting to get a subseq out of range $start:$end vs ",$self->length);
   }
   if( $end < $start ) {
       $self->throw("Attempting to subseq with end ($end) less than start ($start). To revcom use the revcom function with trunc");
   }
   
   my $string;
   if( !$self->_fh->seek($start-1,0) ) {
       $self->throw("Unable to seek on file $start:$end $!");
   }
   my $ret = $self->_fh->read($string,$end-$start+1);
   if( !defined $ret ) {
       $self->throw("Unable to read $start:$end $!");
   }


   return $string;
}


=head2 add_sequence_as_string

 Title   : add_sequence_as_string
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_sequence_as_string{
   my ($self,$str) = @_;

   my $len = $self->length + CORE::length($str);
   $self->_fh->seek(0,2);
   $self->_fh->print($str);
   $self->length($len);
}


=head2 _fh

 Title   : _fh
 Usage   : $obj->_fh($newval)
 Function: 
 Example : 
 Returns : value of _fh
 Args    : newvalue (optional)


=cut

sub _fh{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_fh'} = $value;
    }
    return $obj->{'_fh'};

}

=head2 _filename

 Title   : _filename
 Usage   : $obj->_filename($newval)
 Function: 
 Example : 
 Returns : value of _filename
 Args    : newvalue (optional)


=cut

sub _filename{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_filename'} = $value;
    }
    return $obj->{'_filename'};

}

sub DESTROY {
    my $self = shift;
    if( defined  $self->_fh ) {
	$self->_fh->close();
    }
    unlink($self->_filename);
}


