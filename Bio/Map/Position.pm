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

Bio::Map::Position - Abstracts the notion of multiple positions 
  for a marker within a Map

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
Chad Matsalla, bioinformatics1@dieselwurks.com

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
 Args    : -positions  ArrayRef of tuples, 
                       each tuple contains a Bio::Map::MapI reference
                       and a position value to associate with this map.

=cut

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  $self->{'_positions'} = {};
  
  my ($positions) = $self->_rearrange([qw(POSITIONS)], @args);
  if( defined $positions ) {
      if( ref($positions) =~ /array/i ) { 
	  foreach my $p ( @$positions ) {
	      if( ref($p) !~ /array/i ) {
		  $self->warn("Must specify a tuple array ref of a Bio::Map::MapI and position value for each position when initializing a $class");
		  next;
	      }	      
	      $self->add_position_value(@$p);
	  }
      } else { 
	  $self->warn("Must specify a valid ArrayReference for -positions in $class");
      }
  }
  return $self;
}

=head2 Bio::Map::PositionI methods

=cut

=head2 known_maps

 Title   : known_maps
 Usage   : my @maps = $position->known_maps
 Function: Returns the list of maps that this position has values for
 Returns : list of Bio::Map::MapI unique ids
 Args    : none

=cut

sub known_maps{
   my ($self) = @_;
   return keys %{$self->{'_positions'}};
}

=head2 in_map

 Title   : in_map
 Usage   : if ( $position->in_map($map) ) {}
 Function: Tests if a position has values in a specific map
 Returns : boolean
 Args    : a map unique id OR Bio::Map::MapI


=cut

sub in_map{
   my ($self,$val) = @_;
   return undef if ! defined $val;
   if( ref($val) ) {
       if( ! $val->isa('Bio::Map::MapI') ) {
	   $self->warn('Must either specify a map unique Id or a Bio::Map::MapI when calling in_map');
	   return undef;
       }
       return scalar @{$self->{'_positions'}->{$val->unique_id()}};
   }
   return scalar @{$self->{'_positions'}->{$val}};
}

=head2 each_position_value

 Title   : positions
 Usage   : my @positions = $position->each_position_value($map);
 Function: Retrieve a list of positions coded as strings or ints 
 Returns : Array of position values for a Map
 Args    : Bio::Map::MapI object to get positions for

=cut

sub each_position_value {
   my ($self,$map) = @_;
   if( ! defined $map || ( ref($map) && ! $map->isa('Bio::Map::MapI'))) {
       $self->warn("Must specify a valid Bio::Map::MapI when retrieving each position value");
       return undef;
   }
   my $id = ref($map) ? $map->unique_id : $map;
   if( ! defined $self->{'_positions'}->{$id} ) {
       return ();
   }   
   return @{$self->{'_positions'}->{$id}};
}

=head2 add_position_value

 Title   : add_position_value
 Usage   : $position->add_position_value($map,'100');
 Function: Add a numeric or string position to the PositionI container 
           and assoiciate it with a Bio::Map::MapI
 Returns : none
 Args    : $map - Bio::Map::MapI
           String or Numeric coding for a position on a map

=cut

sub add_position_value{
   my ($self,$map, $value) = @_;
   if( ! defined $map || (ref($map) && ! $map->isa('Bio::Map::MapI')) ) {
       $self->warn("Must specify a valid Bio::Map::MapI when adding a position value");
       return undef;
   }
   if( ref($map) ) {
       push @{$self->{'_positions'}->{$map->unique_id()}}, $value;
   } else {
       push @{$self->{'_positions'}->{$map}}, $value;
   }
   return;
}


=head2 purge_position_values

 Title   : purge_position_values
 Usage   : $position->purge_position_values
 Function: Remove all the position values stored for a position
 Returns : none
 Args    : [optional] only purge values for a given map

=cut

sub purge_position_values{
   my ($self,$map) = @_;
   if( defined $map ) { 
       if( ref($map) && ! $map->isa('Bio::Map::MapI') ) {
	   $self->warn("Must specify a valid Bio::Map::MapI when calling purge_position_values with a parameter");
	   return undef;
       }
       if( ! ref($map) ) {
	   $self->{'_positions'}->{$map} = [];
       } else {
	   $self->{'_positions'}->{$map->unique_id()} = [];
       }
   }
   $self->{'_positions'} = {};
}

=head2 equals

 Title   : equals
 Usage   : if( $mappable->equals($mapable2)) ...
 Function: Test if a position is equal to another position.
 Returns : boolean
 Args    : Bio::Map::PositionI

=cut

sub equals{
   my ($self,$compare) = @_;
   return 0 unless ( defined $compare && ref($compare) &&
		     $compare->isa('Bio::Map::PositionI') );

   foreach my $map ( $self->known_maps ) {
       if( $compare->in_map($map) ) {
	   my @p = $self->each_position_value($map);
	   my @p2 = $compare->each_position_value($map);
   
	   while ( @p && @p2) {
	       return 0 if( pop @p != pop @p2 );
	   }	   
	   return 1;
       }
   }
   return 0;
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
   return 0 unless ( defined $compare && ref($compare) &&
		     $compare->isa('Bio::Map::PositionI') );

   foreach my $map ( $self->known_maps ) {
       if( $compare->in_map($map) ) {
	   my @p = $self->each_position_value($map);
	   my @p2 = $compare->each_position_value($map);
   
	   while ( @p && @p2) {
	       return 0 if( pop @p > pop @p2 );
	   }	   
	   return 1;
       }
   }
   return 0;
}

=head2 greater_than

 Title   : greater_than
 Usage   : if( $mappable->greater_than($m2) ) ...
 Function: Tests if position is greater than another position.
           It is assumed that 2 positions are in the same map.
 Returns : boolean
 Args    : Bio::Map::PositionI

=cut

sub greater_than {
    my ($self,$compare) = @_;
    return 0 unless ( defined $compare && ref($compare) &&
		      $compare->isa('Bio::Map::PositionI') );
    foreach my $map ( $self->known_maps ) {
	if( $compare->in_map($map) ) {
	    my @p = $self->each_position_value($map);
	    my @p2 = $compare->each_position_value($map);

	    while ( @p && @p2) {
		return 0 if( pop @p < pop @p2 );
	    }	   
	    return 1;
	}
    }
    return 0;
}

1;
