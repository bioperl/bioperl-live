# BioPerl module for Bio::Seq::LargeLocatableSeq
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Albert Vilella
#
#	based on the Bio::LargePrimarySeq module
#       by Ewan Birney <birney@sanger.ac.uk>
#
#       and the Bio::LocatableSeq module
#       by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::LargeLocatableSeq - LocatableSeq object that stores sequence as
files in the tempdir

=head1 SYNOPSIS

  # normal primary seq usage
    use Bio::Seq::LargeLocatableSeq;
    my $seq = Bio::Seq::LargeLocatableSeq->new(-seq => "CAGT-GGT",
				              -id  => "seq1",
				              -start => 1,
				              -end   => 7);

=head1 DESCRIPTION

Bio::Seq::LargeLocatableSeq - object with start/end points on it that
can be projected into a MSA or have coordinates relative to another
seq.

This object, unlike Bio::LocatableSeq, stores a sequence as a series
of files in a temporary directory. The aim is to allow someone the
ability to store very large sequences (eg, E<gt> 100MBases) in a file
system without running out of memory (eg, on a 64 MB real memory
machine!).

Of course, to actually make use of this functionality, the programs
which use this object B<must> not call $primary_seq-E<gt>seq otherwise
the entire sequence will come out into memory and probably crash your
machine. However, calls like $primary_seq-E<gt>subseq(10,100) will cause
only 90 characters to be brought into real memory.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Albert Vilella

Email avilella-AT-gmail-DOT-com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Seq::LargeLocatableSeq;
use vars qw($AUTOLOAD);
use strict;


use base qw(Bio::Seq::LargePrimarySeq Bio::LocatableSeq Bio::Root::IO);


=head2 new

 Title   : new
 Usage   : my $obj = Bio::Seq::LargeLocatableSeq->new();
 Function: Builds a new Bio::Seq::LargeLocatableSeq object
 Returns : an instance of Bio::Seq::LargeLocatableSeq
 Args    :


=cut

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
    my $mapping = exists $params{'-mapping'} ? $params{'-mapping'} : [1,1];
    $self->mapping($mapping);
    $self->_initialize_io(%params);
    my $tempdir = $self->tempdir( CLEANUP => 1);
    my ($tfh,$file) = $self->tempfile( DIR => $tempdir );

    $tfh     && $self->_fh($tfh);
    $file    && $self->_filename($file);
    $self->length(0);
    $seq && $self->seq($seq);

    return $self;
}


=head2 length

 Title   : length
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


=cut

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
	   $self->warn("Trying to reset the seq string, cannot do this with a LargeLocatableSeq - must allocate a new object");
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
   my $string;
   my $fh = $self->_fh();

   if( ref($start) && $start->isa('Bio::LocationI') ) {
       my $loc = $start;
       if( $loc->length == 0 ) {
	   $self->warn("Expect location lengths to be > 0");
	   return '';
       } elsif( $loc->end < $loc->start ) {
	   # what about circular seqs
	   $self->warn("Expect location start to come before location end");
       }
       my $seq = '';
       if( $loc->isa('Bio::Location::SplitLocationI') ) {
	   foreach my $subloc ( $loc->sub_Location ) {
	       if(! seek($fh,$subloc->start() - 1,0)) {
		   $self->throw("Unable to seek on file $start:$end $!");
	       }
	       my $ret = read($fh, $string, $subloc->length());
	       if( !defined $ret ) {
		   $self->throw("Unable to read $start:$end $!");
	       }
	       if( $subloc->strand < 0 ) {
# 		   $string = Bio::PrimarySeq->new(-seq => $string)->revcom()->seq();
		   $string = Bio::Seq::LargePrimarySeq->new(-seq => $string)->revcom()->seq();
	       }
	       $seq .= $string;		
	   }
       } else {
	   if(! seek($fh,$loc->start()-1,0)) {
	       $self->throw("Unable to seek on file ".$loc->start.":".
			    $loc->end ." $!");
	   }
	   my $ret = read($fh, $string, $loc->length());
	   if( !defined $ret ) {
	       $self->throw("Unable to read ".$loc->start.":".
			    $loc->end ." $!");
	   }
	   $seq = $string;
       }
       if( defined $loc->strand &&
	   $loc->strand < 0 ) {
# 	   $seq = Bio::PrimarySeq->new(-seq => $string)->revcom()->seq();
	   $seq = Bio::Seq::LargePrimarySeq->new(-seq => $seq)->revcom()->seq();
       }
       return $seq;
   }
   if( $start <= 0 || $end > $self->length ) {
       $self->throw("Attempting to get a subseq out of range $start:$end vs ".
		    $self->length);
   }
   if( $end < $start ) {
       $self->throw("Attempting to subseq with end ($end) less than start ($start). To revcom use the revcom function with trunc");
   }

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
 Usage   : $seq->add_sequence_as_string("CATGAT");
 Function: Appends additional residues to an existing LargeLocatableSeq object.
           This allows one to build up a large sequence without storing
           entire object in memory.
 Returns : Current length of sequence
 Args    : string to append

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


=head2 alphabet

 Title   : alphabet
 Usage   : $obj->alphabet($newval)
 Function:
 Example :
 Returns : value of alphabet
 Args    : newvalue (optional)


=cut

sub alphabet{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->SUPER::alphabet($value);
    }
    return $self->SUPER::alphabet() || 'dna';

}

sub DESTROY {
    my $self = shift;
    my $fh = $self->_fh();
    close($fh) if( defined $fh );
    # this should be handled by Tempfile removal, but we'll unlink anyways.
    unlink $self->_filename() if defined $self->_filename() && -e $self->_filename;
    $self->SUPER::DESTROY();
}

1;
