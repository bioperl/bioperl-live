# $Id$
#
# BioPerl module for Bio::Seq::LargePrimarySeq
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself
#
# updated to utilize File::Temp - jason 2000-12-12
# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::LargePrimarySeq - PrimarySeq object that stores sequence as
files in the tempdir (as found by File::Temp) or the default method in
Bio::Root::RootI

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

  bioperl-l@bioperl.org               - General discussion
  http://bio.perl.org/MailList.html   - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney, Jason Stajich

Email birney@ebi.ac.uk
Email jason@chg.mc.duke.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::LargePrimarySeq;
use vars qw($AUTOLOAD @ISA);
use strict;

use Bio::PrimarySeq;
use Bio::Root::IO;

@ISA = qw(Bio::PrimarySeq Bio::Root::IO);

sub new {
    my ($class, %params) = @_;
    
    # don't let PrimarySeq set seq until we have 
    # opened filehandle

    my $seq = $params{'-seq'} || $params{'-SEQ'};
    if($seq ) {
	delete $params{'-seq'};
	delete $params{'-SEQ'};
    }
    my $self = $class->SUPER::new(%params);
    $self->_initialize_io(%params);
    my $tempdir = $self->tempdir( CLEANUP => 1);
    my ($tfh,$file) = $self->tempfile( DIR => $tempdir );

    $tfh     && $self->_fh($tfh);
    $file    && $self->_filename($file);    
    $self->length(0);
    $seq && $self->seq($seq); 

    return $self;
}


sub length {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'length'} = $value;
    }
   
    return (defined $obj->{'length'}) ? $obj->{'length'} : 0;
}

=head2 seq

 Title   : seq
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub seq {
   my ($self, $data) = @_;   
   if( defined $data ) {
       if( $self->length() == 0) {
	   $self->add_sequence_as_string($data);
       } else { 
	   $self->warn("Trying to reset the seq string, cannot do this with a LargePrimarySeq - must allocate a new object");
       }
   } 
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
   my $fh = $self->_fh();
   if(! seek($fh,$start-1,0)) {
       $self->throw("Unable to seek on file $start:$end $!");
   }
   my $ret = read($fh, $string, $end-$start+1);
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
   my $fh = $self->_fh();
   if(! seek($fh,0,2)) {
       $self->throw("Unable to seek end of file: $!");
   }
   $self->_print($str);
   $self->length($len);
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
    my $fh = $self->_fh();
    close($fh) if( defined $fh );
    # this should be handled by Tempfile removal, but we'll unlink anyways.
    unlink $self->_filename()
    $self->SUPER::DESTROY();
}

1;
