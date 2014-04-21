#
# BioPerl module for Bio::Map::PositionHandler
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Sendu Bala
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::PositionHandler - A Position Handler Implementation

=head1 SYNOPSIS

    # This is used by modules when they want to implement being a
    # Position or being something that has Positions (when they are
    # a L<Bio::Map::EntityI>)

    # Make a PositionHandler that knows about you
    my $ph = Bio::Map::PositionHandler->new($self);

    # Register with it so that it handles your Position-related needs
    $ph->register;

    # If you are a position, get/set the map you are on and the marker you are
    # for
    $ph->map($map);
    $ph->element($marker);
    my $map = $ph->map;
    my $marker = $ph->element;

    # If you are a marker, add a new position to yourself
    $ph->add_positions($pos);

    # And then get all your positions on a particular map
    foreach my $pos ($ph->get_positions($map)) {
        # do something with this Bio::Map::PositionI
    }

    # Or find out what maps you exist on
    my @maps = $ph->get_other_entities;

    # The same applies if you were a map

=head1 DESCRIPTION

A Position Handler copes with the coordination of different Bio::Map::EntityI
objects, adding and removing them from each other and knowning who belongs to
who. These relationships between objects are based around shared Positions,
hence PositionHandler.

This PositionHandler is able to cope with Bio::Map::PositionI objects,
Bio::Map::MappableI objects and Bio::Map::MapI objects.

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

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::PositionHandler;
use strict;

use base qw(Bio::Root::Root Bio::Map::PositionHandlerI);

# globally accessible hash, via private instance methods
my $RELATIONS = {};

=head2 General methods

=cut

=head2 new

 Title   : new
 Usage   : my $position_handler = Bio::Map::PositionHandler->new(-self => $self);
 Function: Get a Bio::Map::PositionHandler that knows who you are.
 Returns : Bio::Map::PositionHandler object
 Args    : -self => Bio::Map::EntityI that is you

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    my ($you) = $self->_rearrange([qw(SELF)], @args);
    
    $self->throw('Must supply -self') unless $you;
    $self->throw('-self must be a reference (object)') unless ref($you);
    $self->throw('This is [$you], not a Bio::Map::EntityI object') unless $you->isa('Bio::Map::EntityI');
    $self->{_who} = $you;
    $self->{_rel} = $RELATIONS;
    return $self;
}

=head2 register

 Title   : register
 Usage   : $position_handler->register();
 Function: Ask this Position Handler to look after your entity relationships.
 Returns : n/a
 Args    : none

=cut

sub register {
    my $self = shift;
    my $you = $self->{_who};
    
    $self->throw("Trying to re-register [$you], which could be bad") if $you->get_position_handler->index;
    
    $self->{_index} = ++$self->{_rel}->{assigned_indices};
    $self->{_rel}->{registered}->{$self->{_index}} = $you;
}

=head2 index

 Title   : index
 Usage   : my $index = $position_handler->index();
 Function: Get the unique registry index for yourself, generated during the
           resistration process.
 Returns : int
 Args    : none

=cut

sub index {
    my $self = shift;
    return $self->{_index};
}

=head2 get_entity

 Title   : get_entity
 Usage   : my $entity = $position_handler->get_entity($index);
 Function: Get the entity that corresponds to the supplied registry index.
 Returns : Bio::Map::EntityI object
 Args    : int

=cut

sub get_entity {
    my ($self, $index) = @_;
    return $self->{_rel}->{registered}->{$index} || $self->throw("Requested registy index '$index' but that index isn't in the registry");
}

=head2 Methods for Bio::Map::PositionI objects

=cut

=head2 map

 Title   : map
 Usage   : my $map = $position_handler->map();
           $position_handler->map($map);
 Function: Get/Set the map you are on. You must be a Position.
 Returns : L<Bio::Map::MapI>
 Args    : none to get, OR
           new L<Bio::Map::MapI> to set

=cut

sub map {
    my ($self, $entity) = @_;
    return $self->_pos_get_set($entity, 'position_maps', 'Bio::Map::MapI');
}

=head2 element

 Title   : element
 Usage   : my $element = $position_handler->element();
           $position_handler->element($element);
 Function: Get/Set the map element you are for. You must be a Position.
 Returns : L<Bio::Map::MappableI>
 Args    : none to get, OR
           new L<Bio::Map::MappableI> to set

=cut

sub element {
    my ($self, $entity) = @_;
    return $self->_pos_get_set($entity, 'position_elements', 'Bio::Map::MappableI');
}

=head2 Methods for all other Bio::Map::EntityI objects

=cut

=head2 add_positions

 Title   : add_positions
 Usage   : $position_handler->add_positions($pos1, $pos2, ...);
 Function: Add some positions to yourself. You can't be a position.
 Returns : n/a
 Args    : Array of Bio::Map::PositionI objects

=cut

sub add_positions {
    my $self = shift;
    $self->throw('Must supply at least one Bio::Map::EntityI') unless @_ > 0;
    my $you_index = $self->_get_you_index(0);
    my $kind = $self->_get_kind;
    
    foreach my $pos (@_) {
        $self->_check_object($pos, 'Bio::Map::PositionI');
        my $pos_index = $self->_get_other_index($pos);
        
        $self->_pos_set($pos_index, $you_index, $kind);
    }
}

=head2 get_positions

 Title   : get_positions
 Usage   : my @positions = $position_handler->get_positions();
 Function: Get all your positions. You can't be a Position.
 Returns : Array of Bio::Map::PositionI objects
 Args    : none for all, OR
           Bio::Map::EntityI object to limit the Positions to those that
           are shared by you and this other entity.

=cut

sub get_positions {
    my ($self, $entity) = @_;
    my $you_index = $self->_get_you_index(0);
    
    my @positions = keys %{$self->{_rel}->{has}->{$you_index}};
    
    if ($entity) {
        my $entity_index = $self->_get_other_index($entity);
        my $pos_ref = $self->{_rel}->{has}->{$entity_index};
        @positions = grep { $pos_ref->{$_} } @positions;
    }
    
    return map { $self->get_entity($_) } @positions;
}

=head2 purge_positions

 Title   : purge_positions
 Usage   : $position_handler->purge_positions();
 Function: Remove all positions from yourself. You can't be a Position.
 Returns : n/a
 Args    : none to remove all, OR
           Bio::Map::PositionI object to remove only that entity, OR
           Bio::Map::EntityI object to limit the removal to those Positions that
           are shared by you and this other entity.

=cut

sub purge_positions {
    my ($self, $thing) = @_;
    my $you_index = $self->_get_you_index(0);
    my $kind = $self->_get_kind;
    
    my @pos_indices;
    if ($thing) {
        $self->throw("Must supply an object") unless ref($thing);
        if ($thing->isa("Bio::Map::PositionI")) {
            @pos_indices = ($self->_get_other_index($thing));
        }
        else {
            my $entity_index = $self->_get_other_index($thing);
            my $pos_ref = $self->{_rel}->{has}->{$entity_index};
            @pos_indices = grep { $pos_ref->{$_} } keys %{$self->{_rel}->{has}->{$you_index}};
        }
    }
    else {
        @pos_indices = keys %{$self->{_rel}->{has}->{$you_index}};
    }
    
    foreach my $pos_index (@pos_indices) {
        $self->_purge_pos_entity($pos_index, $you_index, $kind);
    }
}

=head2 get_other_entities

 Title   : get_other_entities
 Usage   : my @entities = $position_handler->get_other_entities();
 Function: Get all the entities that share your Positions. You can't be a
           Position.
 Returns : Array of Bio::Map::EntityI objects
 Args    : none

=cut

sub get_other_entities {
    my $self = shift;
    my $you_index = $self->_get_you_index(0);
    my $kind = $self->_get_kind;
    my $want = $kind eq 'position_elements' ? 'position_maps' : 'position_elements';
    
    my %entities;
    while (my ($pos_index) = each %{$self->{_rel}->{has}->{$you_index}}) {
        my $entity_index = $self->{_rel}->{$want}->{$pos_index};
        $entities{$entity_index} = 1 if $entity_index;
    }
    
    return map { $self->get_entity($_) } keys %entities;
}

# do basic check on an object, make sure it is the right type
sub _check_object {
    my ($self, $object, $interface) = @_;
    $self->throw("Must supply an arguement") unless $object;
    $self->throw("This is [$object], not an object") unless ref($object);
    $self->throw("This is [$object], not a $interface") unless $object->isa($interface);
}

# get the object we are the handler of, its index, and throw depending on if
# we're a Position
sub _get_you_index {
    my ($self, $should_be_pos) = @_;
    my $you = $self->{_who};
    if ($should_be_pos) {
        $self->throw("This is not a Position, method invalid") unless $you->isa('Bio::Map::PositionI');
    }
    else {
        $self->throw("This is a Position, method invalid") if $you->isa('Bio::Map::PositionI');
    }
    return $self->index;
}

# check an entity is registered and get its index
sub _get_other_index {
    my ($self, $entity) = @_;
    $self->throw("Must supply an object") unless ref($entity);
    my $index = $entity->get_position_handler->index;
    $self->throw("Entity doesn't seem like it's been registered") unless $index;
    $self->throw("Entity may have been registered with a different PositionHandler, can't deal with it") unless $entity eq $self->get_entity($index);
    return $index;
}

# which of the position hashes should we be recorded under?
sub _get_kind {
    my $self = shift;
    my $you = $self->{_who};
    return $you->isa('Bio::Map::MapI') ? 'position_maps' : $you->isa('Bio::Map::MappableI') ? 'position_elements' : $self->throw("This is [$you] which is an unsupported kind of entity");
}

# get/set position entity
sub _pos_get_set {
    my ($self, $entity, $kind, $interface) = @_;
    my $you_index = $self->_get_you_index(1);
    
    my $entity_index;
    if ($entity) {
        $self->_check_object($entity, $interface);
        my $new_entity_index = $self->_get_other_index($entity);
        $entity_index = $self->_pos_set($you_index, $new_entity_index, $kind);
    }
    
    $entity_index ||= $self->{_rel}->{$kind}->{$you_index} || 0;
    if ($entity_index) {
        return $self->get_entity($entity_index);
    }
    return;
}

# set position entity
sub _pos_set {
    my ($self, $pos_index, $new_entity_index, $kind) = @_;
    my $current_entity_index = $self->{_rel}->{$kind}->{$pos_index} || 0;
    
    if ($current_entity_index) {
        if ($current_entity_index == $new_entity_index) {
            return $current_entity_index;
        }
        
        $self->_purge_pos_entity($pos_index, $current_entity_index, $kind);
    }
    
    $self->{_rel}->{has}->{$new_entity_index}->{$pos_index} = 1;
    $self->{_rel}->{$kind}->{$pos_index} = $new_entity_index;
    return $new_entity_index;
}

# disassociate position from one of its current entities
sub _purge_pos_entity {
    my ($self, $pos_index, $entity_index, $kind) = @_;
    delete $self->{_rel}->{has}->{$entity_index}->{$pos_index};
    delete $self->{_rel}->{$kind}->{$pos_index};
}

1;
