# $Id$
#
# BioPerl module for Bio::Map::PositionI
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::PositionI - Abstracts the notion of an element having
  multiple positions for a marker within a Map

=head1 SYNOPSIS

    # get a Bio::Map::PositionI object somehow
    print "The marker maps to the following positions: ", 
    join(",", $position->each_position_value), "\n";

=head1 DESCRIPTION

This object handles the concept that a mappable object (e.g. Marker)
may have multiple positions in a map (e.g. restriction enzymes).  

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


package Bio::Map::PositionI;
use vars qw(@ISA);
use strict;
use Bio::Root::RootI;
use Carp;

@ISA = qw(Bio::Root::RootI);

=head2 each_position_value

 Title   : positions
 Usage   : my @positions = $position->each_position_value('mapname');
 Function: Retrieve a list of positions coded as strings or ints 
 Returns : Array of position values 
 Args    : none

=cut

sub each_position_value{
   my ($self,$mapname) = @_;
   $self->_abstractDeath('each_position_value');
}

=head2 add_position_value

 Title   : add_position_value
 Usage   : $position->add_position_value($map,'100')
 Function: Add a numeric or string position to the PositionI container
 Returns : none
 Args    : Map - Reference to Bio::Map::MapI 
           String or Numeric coding for a position on a map

=cut

sub add_position_value{
   my ($self,$map,$value) = @_;
   $self->_abstractDeath('add_position_value');
}

=head2 purge_position_values

 Title   : purge_position_values
 Usage   : $position->purge_position_values
 Function: Remove all the position values stored for a position
 Returns : none
 Args    : [optional] only purge values for a given map

=cut

sub purge_position_values{
   my ($self) = @_;
   $self->_abstractDeath('purge_position_values');
}

=head2 known_maps

 Title   : known_maps
 Usage   : my @maps = $position->known_maps
 Function: Returns the list of maps that this position has values for
 Returns : list of Bio::Map::MapI unique ids
 Args    : none

=cut

sub known_maps{
   my ($self) = @_;
   $self->_abstractDeath('known_maps');
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
   $self->_abstractDeath('in_map');
}

=head2 equals

 Title   : equals
 Usage   : if( $mappable->equals($mapable2)) ...
 Function: Test if a position is equal to another position.
           It is assumed that 2 positions are in the same map.
 Returns : boolean
 Args    : Bio::Map::PositionI

=cut

sub equals{
   my ($self,$compare) = @_;
   $self->_abstractDeath('less_than');
}

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
   $self->_abstractDeath('less_than');
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
   $self->_abstractDeath('greater_than');
}

1;
