#
# BioPerl module for Bio::Map::MappableI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

    # do not use this module directly
    # See Bio::Map::Mappable for an example of
    # implementation.

=head1 DESCRIPTION

This object handles the generic notion of an element placed on a
(linear) Map. A Mappable can have multiple positions in multiple maps, such as
is the case of Restriction enzyme cut sites on sequence maps. For exact
information about a mappable's position in a map one must query the associate
PositionI objects which are accessible through the get_positions() method.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

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
use strict;
use Bio::Map::PositionHandler;

use base qw(Bio::Map::EntityI Bio::AnnotatableI);

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
        my $ph = Bio::Map::PositionHandler->new(-self => $self);
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
           L<Bio::Map::MapI> object for positions on the given map, AND/OR some
           other true value to avoid sorting

=cut

sub get_positions {
    my ($self, $thing, $no_sort) = @_;
    my $map;
    if (ref($thing) && $thing->isa('Bio::Map::MapI')) {
        $map = $thing;
    }
    else {
        $no_sort = $thing;
    }
    my @positions = $self->get_position_handler->get_positions($map);
    return @positions if @positions == 1;
    
    unless ($no_sort) {
        # don't do
        # @positions = sort { $a->sortable <=> $b->sortable } @positions;
        # directly since sortable() can result in the call of another sort
        # routine and cause problems; pre-compute sortable values instead
        # (which is also more efficient)
        @positions = map { $_->[1] }
                     sort { $a->[0] <=> $b->[0] }
                     map  { [$_->sortable, $_] }
                     @positions;
    }
    return @positions;
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

=head2 name

 Title   : name
 Usage   : my $name = $marker->name();
           $marker->name($new_name);
 Function: Get/Set the name for this Mappable.
 Returns : A scalar representing the current name of this Mappable
 Args    : none to get
           string to set

=cut

sub name {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 id

 Title   : id
 Usage   : my $id = $marker->id();
           $marker->id($new_id);
 Function: Get/Set the id for this Mappable.
 Returns : A scalar representing the current id of this Mappable
 Args    : none to get
           string to set

=cut

sub id {
    my $self = shift;
    $self->throw_not_implemented();
}

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
           comparison positions.
 Returns : array of L<Bio::Map::PositionI> objects
 Args    : arg #1 = L<Bio::Map::MappableI> OR L<Bio::Map::PositionI> to compare
                    this one to (mandatory)
           arg #2 = optionally, the key => value pairs below
		   -map => Bio::Map::MapI           : optionally a Map to only consider
		                                      positions on the given map
		   -relative => Bio::Map::RelativeI : optionally a Relative to ask if
											  the Positions equal in terms of
											  their relative position to the
											  thing described by that Relative

=cut

sub equals {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 overlaps

 Title   : overlaps
 Usage   : if ($mappable->overlaps($other_mappable)) {...}
           my @overlapping_positions = $mappable->overlaps($other_mappable);
 Function: Finds the positions in this mappable that overlap with any
           comparison positions.
 Returns : array of L<Bio::Map::PositionI> objects
 Args    : arg #1 = L<Bio::Map::MappableI> OR L<Bio::Map::PositionI> to compare
                    this one to (mandatory)
           arg #2 = optionally, the key => value pairs below
		   -map => Bio::Map::MapI           : optionally a Map to only consider
		                                      positions on the given map
		   -relative => Bio::Map::RelativeI : optionally a Relative to ask if
                                              the Positions overlap in terms of
                                              their relative position to the
                                              thing described by that Relative

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
           positions.
 Returns : array of L<Bio::Map::PositionI> objects
 Args    : arg #1 = L<Bio::Map::MappableI> OR L<Bio::Map::PositionI> to compare
                    this one to (mandatory)
           arg #2 = optionally, the key => value pairs below
		   -map => Bio::Map::MapI           : optionally a Map to only consider
		                                      positions on the given map
		   -relative => Bio::Map::RelativeI : optionally a Relative to ask if
                                              the Positions contains in terms of
                                              their relative position to the
                                              thing described by that Relative

=cut

sub contains {
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
 Args    : arg #1 = L<Bio::Map::MappableI> OR L<Bio::Map::PositionI> to  compare
                    this one to, or an array ref of such objects (mandatory)
           arg #2 = optionally, the key => value pairs below
		   -map => Bio::Map::MapI           : optionally a Map to only consider
		                                      positions on the given map
		   -relative => Bio::Map::RelativeI : optionally a Relative to to ask
											  how the Positions intersect in
											  terms of their relative position
											  to the thing described by that
											  Relative

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
 Args    : arg #1 = L<Bio::Map::MappableI> OR L<Bio::Map::PositionI> to  compare
                    this one to, or an array ref of such objects (mandatory)
           arg #2 = optionally, the key => value pairs below
		   -map => Bio::Map::MapI           : optionally a Map to only consider
		                                      positions on the given map
		   -relative => Bio::Map::RelativeI : optionally a Relative to to ask
											  if the union of the Positions in
											  terms of their relative position
											  to the thing described by that
											  Relative

=cut

sub union {
    my $self = shift;
    $self->throw_not_implemented();
}

1;
