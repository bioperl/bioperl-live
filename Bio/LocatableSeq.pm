#
# BioPerl module for Bio::LocatableSeq
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::LocatableSeq - A Sequence object with start/end points on it

=head1 SYNOPSIS

    # a normal sequence object
    $locseq->seq();
    $locseq->id();

    # has start,end points
    $locseq->start();
    $locseq->end();
    
    # inheriets off RangeI, so range operations possible

    $locseq->overlaps($seqfeature);

=head1 DESCRIPTION

The locatable sequence object was developed mainly because the 
SimpleAlign object requires this functionality, and in the rewrite
of the Sequence object we had to decide what to do with this.

It is, to be honest, not well integrated with the rest of bioperl, for
example, the ->trunc function does not return a LocatableSeq object,
as some might have thought. There are all sorts of nasty gotcha's about
interactions between coordinate systems when these sort of objects are
used. 

Please post to the guts list for more about this, or contact Ewan Birney.

=head1 CONTACT

Ewan Birney <birney@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::LocatableSeq;
use vars qw(@ISA);
use strict;

use Bio::Seq;
use Bio::RangeI;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;

@ISA = qw(Bio::Seq Bio::RangeI);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my ($start,$end,$strand) = $self->_rearrange( [qw(START END STRAND)],@args);
  my $make = $self->SUPER::_initialize(@args);
 
  defined $start && $self->start($start);
  defined $end   && $self->end($end);
  defined $strand && $self->strand($strand);

# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 start

 Title   : start
 Usage   : $obj->start($newval)
 Function: 
 Returns : value of start
 Args    : newvalue (optional)

=cut

sub start{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'start'} = $value;
    }
    return $obj->{'start'};

}

=head2 end

 Title   : end
 Usage   : $obj->end($newval)
 Function: 
 Returns : value of end
 Args    : newvalue (optional)

=cut

sub end{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'end'} = $value;
    }
    return $obj->{'end'};

}

=head2 strand

 Title   : strand
 Usage   : $obj->strand($newval)
 Function: 
 Returns : value of strand
 Args    : newvalue (optional)

=cut

sub strand{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'strand'} = $value;
    }
    return $obj->{'strand'};

}

=head2 get_nse

 Title   : get_nse
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub get_nse{
   my ($self,$char1,$char2) = @_;
  
   if( !defined $char1 ) { $char1 = "/"; }
   if( !defined $char2 ) { $char2 = "-"; }

   return $self->id() . $char1 . $self->start . $char2 . $self->end ;

}

