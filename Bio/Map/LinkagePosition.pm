# BioPerl module for Bio::Map::LinkagePosition
#
# Cared for by Chad Matsalla <bioinformatics1@dieselwurks.com>
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::LinkagePosition - Create a Position for a Marker that will be placed
	on a Bio::Map::LinkageMap

=head1 SYNOPSIS

    use Bio::Map::Position;
    my $position = new Bio::Map::LinkagePosition(-positions => 1,
						 -distance => 22.1 );
    
	    # can get listing of positions
    my @positions = $position->each_position;


=head1 DESCRIPTION

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

=head1 AUTHOR - Chad Matsalla

Email bioinformatics1@dieselwurks.com

=head1 CONTRIBUTORS

Lincoln Stein, lstein@cshl.org
Heikki Lehvaslaiho, heikki@ebi.ac.uk
Jason Stajich jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Map::LinkagePosition;
use vars qw(@ISA);
use strict;
require 'dumpvar.pl';

use Bio::Map::OrderedPosition;

@ISA = qw(Bio::Map::OrderedPosition);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Map::LinkagePosition(-positions => $position,
				-distance => $distance );
 Function: Builds a new Bio::Map::LinkagePosition object 
 Returns : Bio::Map::LinkagePosition
 Args    : -order => the relative order of this marker on a linkage map
 	   -positions => positions on a map
=cut

=head2 Bio::Map::OrderedPosition methods

=head2 order

 Title   : order
 Usage   : $o_position->order($new_position) _or_
           $o_position->order()
 Function: get/set the order position of this position in a map
 Returns : 
 Args    : If $new_position is provided, the current position of this Position
           will be set to $new_position.

=cut

=head2 Bio::Map::Position functions

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

=cut

sub each_position_value {
    my ($self, @args) = @_;
    my @v = $self->SUPER::each_position_value(@args);
    if ( ! @v ) { 
	return ('0.0');
    }
    return @v;
}

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


=head2 greater_than

 Title   : greater_than
 Usage   : if( $mappable->greater_than($m2) ) ...
 Function: Tests if position is greater than another position.
           It is assumed that 2 positions are in the same map.
 Returns : boolean
 Args    : Bio::Map::PositionI

=cut


1;
