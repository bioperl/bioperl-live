# $Id$
#
# BioPerl module for Bio::Map::MappableI
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::MappableI - An object that can be placed in a map

=head1 SYNOPSIS

    # get a Bio::Map::MappableI somehow
    my $position = $element->position();
    # these methods will be important for building sorted lists
    if( $position->equals($p2) ) {
	# do something
    } elsif( $position->less_tha($p2) ) {} 
      elsif( $position->greater_than($p2) ) { }    


=head1 DESCRIPTION

This object handles the generic notion of an element placed on a
(linear) Map. A Mappable can have multiple positions in multiple maps, such as
is the case of Restriction enzyme cut sites on sequence maps. For exact
information about a mappable's position in a map one must query the associate
PositionI objects which are accessible through the each_position() method.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 CONTRIBUTORS

Heikki Lehvaslaiho heikki-at-bioperl-dot-org
Lincoln Stein      lstein@cshl.org
Sendu Bala         bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::MappableI;
use vars qw(@ISA);
use strict;
use Bio::Map::EntityI;
use Bio::Map::PositionHandler;

@ISA = qw(Bio::Map::EntityI);

=head2 EntityI methods

 These are fundamental to coordination of Mappables and other entities, so are
 implemented at the interface level

=cut

=head2 get_position_handler

 Title   : get_position_handler
 Usage   : my $position_handler = $entity->get_position_handler();
 Function: Gets a PositionHandlerI that $entity is registered with.
 Returns : Bio::Map::PositionHandlerI object
 Args    : none

=cut

sub get_position_handler {
    my $self = shift;
    unless (defined $self->{_eh}) {
        my $ph = Bio::Map::PositionHandler->new($self);
        $self->{_eh} = $ph;
        $ph->register;
    }
    return $self->{_eh};
}

=head2 PositionHandlerI-related methods

 These are fundamental to coordination of Mappables and other entities, so are
 implemented at the interface level

=cut

=head2 add_position

 Title   : add_position
 Usage   : $mappable->add_position($position);
 Function: Add a position to this mappable (defining where on which map it is).
 Returns : n/a
 Args    : L<Bio::Map::PositionI> object

=cut

sub add_position {
    my $self = shift;
	# actually, we allow multiple positions to be set at once
    $self->get_position_handler->add_positions(@_);
}

=head2 get_positions

 Title   : get_positions
 Usage   : my @positions = $mappable->get_positions();
 Function: Get all the Positions of this Mappable (sorted).
 Returns : Array of L<Bio::Map::PositionI> objects
 Args    : none for all, OR
           L<Bio::Map::MapI> object for positions on the given map

=cut

sub get_positions {
    my ($self, $map) = @_;
    my @positions = $self->get_position_handler->get_positions($map);
	return sort { $a->sortable <=> $b->sortable } @positions;
}

=head2 each_position

 Title   : each_position
 Function: Synonym of the get_positions() method.
 Status  : deprecated, will be removed in next version

=cut

*each_position = \&get_positions;

=head2 purge_positions

 Title   : purge_positions
 Usage   : $mappable->purge_positions();
 Function: Remove positions from this mappable.
 Returns : n/a
 Args    : none to remove all positions, OR
           L<Bio::Map::PositionI> object to remove just that Position, OR
		   L<Bio::Map::MapI> object to remove only those positions on the given
		   map

=cut

sub purge_positions {
    my ($self, $thing) = @_;
    $self->get_position_handler->purge_positions($thing);
}

=head2 known_maps

 Title   : known_maps
 Usage   : my @maps = $marker->known_maps()
 Function: Returns the maps that this mappable is found on
 Returns : Array of L<Bio::Map::MapI> objects
 Args    : none

=cut

sub known_maps {
	my $self = shift;
	return $self->get_position_handler->get_other_entities;
}

=head2 MappableI-specific methods

=cut

=head2 in_map

 Title   : in_map
 Usage   : if ($marker->in_map($map)) {...}
 Function: Tests if this mappable is found on a specific map
 Returns : boolean
 Args    : L<Bio::Map::MapI>

=cut

sub in_map {
	my $self = shift;
	$self->throw_not_implemented();
}

=head1 RangeI-related Methods

They throw an error if start and end are not defined in the Positions of the
Mappables supplied.

=cut

=head2 equals

 Title   : equals
 Usage   : if ($mappable->equals($other_mappable)) {...}
           my @equal_positions = $mappable->equals($other_mappable);
 Function: Finds the positions in this mappable that are equal to any
           comparison positions, optionally only considering a particular map.
 Returns : array of L<Bio::Map::PositionI> objects
 Args    : L<Bio::Map::MappableI> OR L<Bio::Map::PositionI> [REQUIRED], and
           optionally a L<Bio::Map::MapI> object to only consider positions on
		   the given map.
		   [*** relative handling ***]
		   
=cut

sub equals{
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 overlaps

 Title   : overlaps
 Usage   : if ($mappable->overlaps($other_mappable)) {...}
           my @overlapping_positions = $mappable->overlaps($other_mappable);
 Function: Finds the positions in this mappable that overlap with any
           comparison positions, optionally only considering a particular map.
 Returns : array of L<Bio::Map::PositionI> objects
 Args    : L<Bio::Map::MappableI> OR L<Bio::Map::PositionI> [REQUIRED], and
           optionally a L<Bio::Map::MapI> object to only consider positions on
		   the given map.
		   [*** relative handling ***]

=cut

sub overlaps {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 contains

 Title   : contains
 Usage   : if ($mappable->contains($other_mappable)) {...}
           my @container_positions = $mappable->contains($other_mappable);
 Function: Finds the positions in this mappable that contain any comparison
           positions, optionally only considering a particular map.
 Returns : array of L<Bio::Map::PositionI> objects
 Args    : L<Bio::Map::MappableI> OR L<Bio::Map::PositionI> [REQUIRED], and
           optionally a L<Bio::Map::MapI> object to only consider positions on
		   the given map.
		   [*** relative handling ***]

=cut

sub contains {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 less_than

 Title   : less_than
 Usage   : if ($mappable->less_than($other_mappable)) {...}
           my @lower_positions = $mappable->less_than($other_mappable);
 Function: Finds the positions in this mappable that are less than all
           comparison positions, optionally only considering a particular map.
 Returns : array of L<Bio::Map::PositionI> objects
 Args    : L<Bio::Map::MappableI> OR L<Bio::Map::PositionI> [REQUIRED], and
           optionally a L<Bio::Map::MapI> object to only consider positions on
		   the given map.
		   [*** relative handling ***]

=cut

sub less_than {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 greater_than

 Title   : greater_than
 Usage   : if ($mappable->greater_than($other_mappable)) {...}
           my @higher_positions = $mappable->greater_than($other_mappable);
 Function: Finds the positions in this mappable that are greater than all
           comparison positions, optionally only considering a particular map.
 Returns : array of L<Bio::Map::PositionI> objects
 Args    : L<Bio::Map::MappableI> OR L<Bio::Map::PositionI> [REQUIRED], and
           optionally a L<Bio::Map::MapI> object to only consider positions on
		   the given map.
		   [*** relative handling ***]

=cut

sub greater_than {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 intersection

 Title   : intersection
 Usage   : my $position = $mappable->intersection($other_mappable);
           my $position = Bio::Map::Mappable->intersection(\@mappables);
 Function: Make the position that is at the intersection of all positions of all
           supplied mappables.
 Returns : L<Bio::Map::PositionI> object or undef (if not all positions overlap)
 Args    : arg #1 = [REQUIRED] a L<Bio::Map::MappableI> OR
					L<Bio::Map::PositionI> to compare this one to, or an array
					ref of such objects
           arg #2 = [*** relative handling ***]

=cut

sub intersection {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 union

 Title   : union
 Usage   : my $position = $mappable->union($other_mappable);
           my $position = Bio::Map::Mappable->union(@mappables);
 Function: Make the minimal position that contains all of the positions of all
           supplied mappables.
 Returns : L<Bio::Map::PositionI> object or undef (if not all positions overlap)
 Args    : arg #1 = [REQUIRED] a L<Bio::Map::MappableI> OR
					L<Bio::Map::PositionI> to compare this one to, or an array
					of such objects
           arg #2 = [*** relative handling ***]

=cut

sub union {
    my $self = shift;
    $self->throw_not_implemented();
}

1;
