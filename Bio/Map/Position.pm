# $Id$
#
# BioPerl module for Bio::Map::Position
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::Position - Abstracts the notion of multiple positions for a marker within a Map

=head1 SYNOPSIS

    use Bio::Map::Position;
    my $position = new Bio::Map::Position(-positions => [ 100 ] );
    
    # or can add items to 
    $position->add_position(105);
    
    # can get listing of positions
    my @positions = $position->each_position;


=head1 DESCRIPTION

This object is an implementation of the PositionI interface that
handles the specific values of a position.  This allows an element
(e.g. Marker) to have multiple positions within a map and still be
treated as a single entity.  This does not handle the concept of a
relative map in which no known exact positions exist but markers are
just ordered relative to one another - in that case arbitrary values
must be assigned for positions.

No units are assumed here - units are handled by context of which Map
a position is placed in.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Lincoln Stein, lstein@cshl.org
Heikki Lehvaslaiho, heikki@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Map::Position;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::Map::PositionI;

@ISA = qw(Bio::Root::Root Bio::Map::PositionI );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Map::Position();
 Function: Builds a new Bio::Map::Position object 
 Returns : Bio::Map::Position
 Args    : -positions  ArrayRef or single item to be initialized as
                       positions within this object

=cut

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  $self->{'_positions'} = [];
  my ($positions) = $self->_rearrange([qw(POSITIONS)], @args);
  
  if( ref($positions) =~ /array/i ) { 
      foreach my $p ( @$positions ) {
	  $self->add_position($p);
      }
  } else { 
      $self->add_position($positions);
  }
  return $self;
}


=head2 Bio::Map::PositionI methods

=head2 each_position

 Title   : positions
 Usage   : my @positions = $position->positions();
 Function: Retrieve a list of positions coded as strings or ints 
 Returns : Array of position values 
 Args    : none

=cut

sub each_position{
   my ($self) = @_;
   return @{$self->{'_positions'}};
}

=head2 add_position

 Title   : add_position
 Usage   : $position->add_position('100')
 Function: Add a numeric or string position to the PositionI container
 Returns : none
 Args    : String or Numeric coding for a position on a map

=cut

sub add_position{
   my ($self,$value) = @_;
   if( ! $value ) { 
       $self->warn("Attempting to add a position with a null value");
       return;
   }
   push @{$self->{'_positions'}}, $value;
   return;
}

=head2 purge

 Title   : purge
 Usage   : $position->purge
 Function: Remove all the position values stored for a position
 Returns : none
 Args    : none

=cut

sub purge_positions{
   my ($self) = @_;
   $self->{'_positions'} = [];
}

1;
