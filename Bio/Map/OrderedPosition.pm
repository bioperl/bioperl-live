# BioPerl module for Bio::Map::OrderedPosition
#
# Cared for by Chad Matsalla <bioinformatics1@dieselwurks.com>
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::OrderedPosition - Abstracts the notion of a member
	of an ordered list of markers. Each marker is a certain distance
	from the one in the ordered list before it.

=head1 SYNOPSIS

    use Bio::Map::OrderedPosition;
	# the first marker in the sequence
    my $position = new Bio::Map::OrderedPosition(-order => 1,
			-positions => [ [ $map, 22.3] ] );
	# the second marker in the sequence, 15.6 units from the fist one
    my $position2 = new Bio::Map::OrderedPosition(-order => 2,
			-positions => [ [ $map, 37.9] ] );
	# the third marker in the sequence, coincidental with the second
	# marker
    my $position3 = new Bio::Map::OrderedPosition(-order => 3,
                        -posititions => [ [ $map, 37.9]] );

=head1 DESCRIPTION

This object is an implementation of the PositionI interface and the
Position object handles the specific values of a position.
OrderedPosition is intended to be slightly more specific then Position
but only specific enough for a parser from the MarkerIO subsystem to
create and then pass to a client application to bless into the proper
type. For an example of how this is intended to work, see the
Mapmaker.pm.

No units are assumed here - units are handled by context of which Map
a position is placed in.

Se Bio::Map::Position for additional information.

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
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Chad Matsalla

Email bioinformatics1@dieselwurks.com

=head1 CONTRIBUTORS

Lincoln Stein, lstein@cshl.org
Heikki Lehvaslaiho, heikki@ebi.ac.uk
Jason Stajich, jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Map::OrderedPosition;
use vars qw(@ISA);
use strict;

use Bio::Map::Position;

@ISA = qw(Bio::Map::Position);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Map::OrderedPosition();
 Function: Builds a new Bio::Map::OrderedPosition object 
 Returns : Bio::Map::OrderedPosition
 Args    : -order - The order of this position

=cut

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
#    $self->{'_order'} = [];
  
    my ($map, $marker, $value, $order) = 
	$self->_rearrange([qw( MAP 
			       MARKER 
			       VALUE
			       ORDER
			       )], @args);
#    print join ("|-|", ($map, $marker, $value, $order)), "\n";
    $map     && $self->map($map);
    $marker  && $self->marker($marker);
    $value   && $self->value($value);
    $order   && $self->order($order);

    return $self;
}

=head2 order

 Title   : order
 Usage   : $o_position->order($new_position) _or_
           $o_position->order()
 Function: get/set the order position of this position in a map
 Returns : 
 Args    : If $new_position is provided, the current position of this Position
           will be set to $new_position.

=cut

sub order {
    my ($self,$order) = @_;
    if ($order) {
	# no point in keeping the old ones
	$self->{'_order'} = $order;
    }
    return $self->{'_order'};
}

=head2 Bio::Map::Position functions

=cut

=head2 known_maps

 Title   : known_maps
 Usage   : my @maps = $position->known_maps
 Function: Returns the list of maps that this position has values for
 Returns : list of Bio::Map::MapI unique ids
 Args    : none

=head2 in_map

 Title   : in_map
 Usage   : if ( $position->in_map($map) ) {}
 Function: Tests if a position has values in a specific map
 Returns : boolean
 Args    : a map unique id OR Bio::Map::MapI

=head2 each_position_value

 Title   : positions
 Usage   : my @positions = $position->each_position_value($map);
 Function: Retrieve a list of positions coded as strings or ints 
 Returns : Array of position values for a Map
 Args    : Bio::Map::MapI object to get positions for

=head2 add_position_value

 Title   : add_position_value
 Usage   : $position->add_position_value($map,'100');
 Function: Add a numeric or string position to the PositionI container 
           and assoiciate it with a Bio::Map::MapI
 Returns : none
 Args    : $map - Bio::Map::MapI
           String or Numeric coding for a position on a map

=head2 purge_position_values

 Title   : purge_position_values
 Usage   : $position->purge_position_values
 Function: Remove all the position values stored for a position
 Returns : none
 Args    : [optional] only purge values for a given map

=head2 equals

 Title   : equals
 Usage   : if( $mappable->equals($mapable2)) ...
 Function: Test if a position is equal to another position.
 Returns : boolean
 Args    : Bio::Map::PositionI

=cut

sub equals{
   my ($self,$compare) = @_;
   return 0 if ( ! defined $compare || ! $compare->isa('Bio::Map::OrderedPosition'));
   return ( $compare->order == $self->order);
}

# admittedly these are really the best comparisons in the world
# but it is a first pass we'll need to refine the algorithm or not 
# provide general comparisions and require these to be implemented
# by objects closer to the specific type of data

=head2 less_than

 Title   : less_than
 Usage   : if( $mappable->less_than($m2) ) ...
 Function: Tests if a position is less than another position
           It is assumed that 2 positions are in the same map.
 Returns : boolean
 Args    : Bio::Map::PositionI

=cut


sub less_than{
   my ($self,$compare) = @_;
   return 0 if ( ! defined $compare || ! $compare->isa('Bio::Map::OrderedPosition'));
   return ( $compare->order < $self->order);
}

=head2 greater_than

 Title   : greater_than
 Usage   : if( $mappable->greater_than($m2) ) ...
 Function: Tests if position is greater than another position.
           It is assumed that 2 positions are in the same map.
 Returns : boolean
 Args    : Bio::Map::PositionI

=cut

sub greater_than{
   my ($self,$compare) = @_;
   return 0 if ( ! defined $compare || ! $compare->isa('Bio::Map::OrderedPosition'));
   return ( $compare->order > $self->order);
}

1;
