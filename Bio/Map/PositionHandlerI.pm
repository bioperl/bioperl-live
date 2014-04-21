#
# BioPerl module for Bio::Map::PositionHandlerI
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

Bio::Map::PositionHandlerI - A Position Handler Interface

=head1 SYNOPSIS

    # do not use this module directly
    # See Bio::Map::PositionHandler for an example of
    # implementation.

=head1 DESCRIPTION

This interface describes the basic methods required for Position Handlers. A
Position Handler copes with the coordination of different Bio::Map::EntityI
objects, adding and removing them from each other and knowning who belongs to
who. These relationships between objects are based around shared Positions,
hence PositionHandler.

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

package Bio::Map::PositionHandlerI;
use strict;

use base qw(Bio::Root::RootI);

=head2 General methods

=cut

=head2 register

 Title   : register
 Usage   : $position_handler->register();
 Function: Ask this Position Handler to look after your entity relationships.
 Returns : n/a
 Args    : none

=cut

sub register {
    my $self = shift;
    $self->throw_not_implemented();
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
    $self->throw_not_implemented();
}

=head2 get_entity

 Title   : get_entity
 Usage   : my $entity = $position_handler->get_entity($index);
 Function: Get the entity that corresponds to the supplied registry index.
 Returns : Bio::Map::EntityI object
 Args    : int

=cut

sub get_entity {
        my $self = shift;
    $self->throw_not_implemented();
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
    my $self = shift;
    $self->throw_not_implemented();
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
    my $self = shift;
    $self->throw_not_implemented();
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
    $self->throw_not_implemented();
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
    my $self = shift;
    $self->throw_not_implemented();
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
    my $self = shift;
    $self->throw_not_implemented();
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
    $self->throw_not_implemented();
}

1;
