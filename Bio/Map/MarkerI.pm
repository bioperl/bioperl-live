# $Id$
#
# BioPerl module for Bio::Map::MarkerI
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::MarkerI - Interface for basic marker functionality

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the interface here

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

Heikki Lehvaslaiho heikki@ebi.ac.uk
Lincoln Stein      lstein@cshl.org
Jason Stajich      jason@bioperl.org
Chad Matsalla      bioinformatics1@dieselwurks.com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Map::MarkerI;
use vars qw(@ISA);
use strict;
use Bio::Map::MappableI;

@ISA = qw(Bio::Map::MappableI);

=head2 name($new_name)

 Title   : name($new_name)
 Usage   : my $name = $o_usat->name($new_name) _or_
	   my $name = $o_usat->name()
 Function: Get/Set the name for this Marker
 Returns : A scalar representing the current name of this Marker
 Args    : If provided, the current name of this Marker
	   will be set to $new_name.

=cut

sub name {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 position

 Title   : position
 Usage   : my position_string = $position->position('mapname');
 Function: Get/Set method for single value positions. 
           Gives a simplified interface when only one map and 
           one position per marker is used.
 Returns : a position value
 Args    : optional:
           Map - Reference to Bio::Map::MapI 
           String or Numeric coding for a position on a map

=cut

sub position{
   my ($self,$map, $value) = @_;
   $self->throw_not_implemented();
}

=head2 add_position

 Title   : add_position
 Usage   : $position->add_position($map,'100')
 Function: Add a numeric or string position to the PositionI container
 Returns : none
 Args    : Map - Reference to Bio::Map::MapI 
           String or Numeric coding for a position on a map

=cut

sub add_position{
   my ($self,$map,$value) = @_;
   $self->throw_not_implemented();
}


=head2 each_position

 Title   : positions
 Usage   : my @positions = $position->each_position_value('mapname');
 Function: Retrieve a list of positions coded as strings or ints 
 Returns : Array of position values 
 Args    : none

=cut

sub each_position {
   my ($self,$mapname) = @_;
   $self->throw_not_implemented();
}

=head2 purge_positions

 Title   : purge_positions
 Usage   : $position->purge_positions
 Function: Remove all the position values stored for a Marker
 Returns : none
 Args    : [optional] only purge values for a given map

=cut

sub purge_position_values{
   my ($self, $map) = @_;
   $self->throw_not_implemented();
}

=head2 known_maps

 Title   : known_maps
 Usage   : my @maps = $marker->known_maps
 Function: Returns the list of maps that this position has values for
 Returns : list of Bio::Map::MapI unique ids
 Args    : none

=cut

sub known_maps{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 in_map

 Title   : in_map
 Usage   : if ( $position->in_map($map) ) {}
 Function: Tests if a position has values in a specific map
 Returns : boolean
 Args    : a map unique id OR Bio::Map::MapI


=cut

sub in_map{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 get_position_object

 Title   : get_position_class
 Usage   : my $pos = $marker->get_position_object();
 Function: To get an object of the default Position class
           for this Marker. Subclasses should redefine this method.
           The Position needs to be L<Bio::Map::PositionI>.
 Returns : L<Bio::Map::PositionI>
 Args    : none

=cut

sub get_position_object {
    my ($self) = @_;
    $self->throw_not_implemented();
}


=head2 tuple

 Title   : tuple
 Usage   : ($me, $you) = $self->_tuple($compare)
 Function: Utility method to extract numbers and test for missing values.
           Makes writing subsequent tests easier.
 Returns : a tuple of values or ranges
 Args    : Bio::Map::MappableI or Bio::Map::PositionI

=cut

sub in_map{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}


=head2 Bio::Map::MappableI methods

=cut

=head2 position

 Title   : position
 Usage   : my $position = $mappable->position(); 
 Function: Get/Set the Bio::Map::PositionI for a mappable element
 Returns : Bio::Map::PositionI
 Args    : (optional) Bio::Map::PositionI

=head2 equals

 Title   : equals
 Usage   : if( $mappable->equals($mapable2)) ...
 Function: Test if a position is equal to another position
 Returns : boolean
 Args    : Bio::Map::MappableI

=head2 less_than

 Title   : less_than
 Usage   : if( $mappable->less_than($m2) ) ...
 Function: Tests if a position is less than another position
 Returns : boolean
 Args    : Bio::Map::MappableI

=head2 greater_than

 Title   : greater_than
 Usage   : if( $mappable->greater_than($m2) ) ...
 Function: Tests if position is greater than another position
 Returns : boolean
 Args    : Bio::Map::MappableI

=cut

1;
