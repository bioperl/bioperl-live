# BioPerl module for Bio::Map::Marker
#
# Cared for by Chad Matsalla <bioinformatics1@dieselwurks.com>
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::Marker - An central map object representing a generic marker
that can have multiple location in several maps.

=head1 SYNOPSIS

  # get map objects somehow

  # a marker with complex localisation
  $o_usat = new Bio::Map::Marker (-name=>'Chad Super Marker 2',
				  -positions => [ [$map1, $position1], 
                                                  [$map1, $position2] 
						] );

  # The markers deal with Bio::Map::Position objects which can also
  # be explicitely created and passed on to markers as an array ref:
  $o_usat2 = new Bio::Map::Marker (-name=>'Chad Super Marker 3',
				  -positions => [ $pos1, 
                                                  $pos2 
						] );

  # a marker with unique position in a map
  $marker1 = new Bio::Map::Marker (-name=>'hypervariable1',
				   -map => $map1,
				   -position => 100
				   )

  # an other way of creating a marker with unique position in a map:
  $marker2 = new Bio::Map::Marker (-name=>'hypervariable2');
  $map1->add_marker($marker2);
  $marker2->position(100);

  # position method is a short cut for get/set'ing  unigue positions
  # which overwrites previous values
  # to place a marker to other maps or to have multiple positions 
  # for a map within the same map use add_position()

  $marker2->add_position(200);	# new position in the same map
  $marker2->add_position($map2,200); # new map

  # setting a map() in a marker or adding a marker into a map are
  # identical mathods. Both set the bidirectional connection which is
  # used by the marker to remember its latest, default map.

  # Regardes of how marker positions are created, they are stored and
  # returned as Bio::Map::PositionI objects:

  # unique position
  print $marker1->position->value, "\n";
  # several positions
  foreach $pos ($marker2->each_position($map1)) {
     print $pos->value, "\n";
  }

See L<Bio::Map::Position> and L<Bio::Map::PositionI> for more information.

=head1 DESCRIPTION

This object handles the notion of a generic marker. This marker will
have a name and a position in a map.

This object is intended to be used by a marker parser like Mapmaker.pm
and then blessed into the proper type of marker (ie Microsatellite) by
the calling script.

=head2 Design principles

A Marker is a central object in Bio::Map name space. A Map is a holder
class for objects. A Marker has a Position in a Map.  A Marker can be
compared to an other Markers using boolean methods. Positions can have
non-numeric values or other methods to store the locations, so they
have a method numeric() which does the conversion. 

A Marker has a convinience method position() which is able to create
Positions of required class from scalars by calling method
get_position_object().

For more complex situations, a Marker can have multiple positions in
multiple Maps. It is therefore possible to extract Positions (all or
belonging to certain Map) and compare Markers to them. It is up to the
programmer to make sure position values and Maps they belong to can be
sensibly compared.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing list.  Your participation is much appreciated.

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

Heikki Lehvaslaiho heikki@ebi.ac.uk
Lincoln Stein      lstein@cshl.org
Jason Stajich      jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::Marker;
use vars qw(@ISA);
use strict;
use Bio::Root::Root;
use Bio::Map::MarkerI;
use Bio::Map::Position;

@ISA = qw(Bio::Root::Root Bio::Map::MarkerI);

=head2 new

 Title   : new
 Usage   : $o_marker = new Bio::Map::Marker( -name => 'Whizzy marker',
	                                     -position => $position);
 Function: Builds a new Bio::Map::Marker object
 Returns : Bio::Map::Marker
 Args    :
           -name    => name of this microsatellite 
                       [optional], string,default 'Unknown'

           -positions => map position for this marker, [optional]
                Bio::Map::PositionI-inherited obj, no default)

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    $self->{'_positions'} = [];
    my ($name, $map, $position, $positions) = 
	$self->_rearrange([qw(NAME 
			      MAP
			      POSITION
			      POSITIONS
			      )], @args);
    if ($name) { $self->name($name); } 
    else {$self->name('Unnamed marker'); }
    $position  && $self->position($position); 
    $positions && $self->positions($positions); 
    $map       && $self->map($map);
 
    return $self;
}

=head2 name

 Title   : name
 Usage   : $o_usat->name($new_name) _or_
	   my $name = $o_usat->name()
 Function: Get/Set the name for this Microsatellite
 Returns : A scalar representing the current name of this marker
 Args    : If provided, the current name of this marker
	   will be set to $new_name.

=cut

sub name {
    my ($self,$name) = @_;
    my $last = $self->{'_name'};
    if ($name) {
	$self->{'_name'} = $name;
    }
    return $last;
}


=head2 map

 Title   : map
 Usage   : my $mymap = $marker->map();
 Function: Get/Set the default map for the marker.
           This is set by L<Bio::Map::CytoMap::add_element> method
 Returns : L<Bio::Map::MapI>
 Args    : [optional] new L<Bio::Map::MapI>

=cut

sub map {
   my ($self,$map) = @_;
   if( defined $map ) {
       $self->thow('This is [$map], not Bio::Map::MapI object')
	   unless $map->isa('Bio::Map::MapI');
       $self->{'_default_map'} = $map;
   }
   return $self->{'_default_map'};
}



=head2 get_position_object

 Title   : get_position_class
 Usage   : my $pos = $marker->get_position_object();
 Function: To get an object of the default Position class
           for this Marker. Subclasses should redefine this method.
           The Position needs to be Bio::Map::PositionI.
 Returns : Bio::Map::Position
 Args    : none

See L<Bio::Map::Position> and L<Bio::Map::PositionI> for more information.

=cut

sub get_position_object {
   my ($self) = @_;
   return new Bio::Map::Position();
}


=head2 position

 Title   : position
 Usage   : $position = $mappable->position($map); OR
           $mappable->position($position); # $position can be  Bio::Map::PositionI
           $mappable->position(100); # or scalar if the marker has a default map
           $mappable->position($map, 100); # if not give explicit $map
 Function: Get/Set the Bio::Map::PositionI for a mappable element
           in a specific Map
           Adds the marker to a map automatically if Map is given. 
           Altenaitvely, you can add the merker to the map first 
           (L<Bio::Map::Map::add_element>) to set the default map
 Returns : Bio::Map::PositionI
 Args    : $position - Bio::Map::PositionI # Position we want to set
            OR
           $map  - Bio::Map::MapI AND
           scalar
            OR
           scalar, but only if the marker has been added to a map

=cut

sub position {
    my ($self, $pos, $secondary_pos) = @_;
    my ($map);
  POS: {
      if ($pos) {
	  if (ref($pos) eq 'SCALAR' || ref($pos) eq '') {
	      $map = $self->map;
	  }
	  elsif (ref($pos) eq 'ARRAY') {
	      $map = $pos->[0];
	      $pos = $pos->[1];
	  }
	  elsif ($pos->isa('Bio::Map::PositionI')) {
	      $pos->marker($self);

	      $self->purge_positions;
	      $self->add_position($pos);
	      $map = $pos->map;
	      $map->add_element($self) unless defined($self->map) && $self->map eq $map;
	      last POS;
	  }

	  elsif ($pos->isa('Bio::Map::MapI')) {
	      $map = $pos;
	      $pos = $secondary_pos;
	  } else {
	      $map = $self->map;
	  }
	  $self->throw("You need to add a marker to a map before ". 
		       "you can set positions without explicit map!" )
	      unless $map;
	  $self->throw("Position better be scalar, not [$pos=". ref($pos)  ."]")
	      unless ref($pos) eq 'SCALAR' || ref($pos) eq ''; 

	  my $newpos = $self->get_position_object;
	  $newpos->map($map);
	  $newpos->value($pos);
	  $newpos->marker($self);

	  $map->add_element($self) unless defined($self->map) && $self->map eq $map;
	  $self->purge_positions;
	  $self->add_position($newpos)
	  }
  }
    my @array = $self->each_position();
    $self->warn('More than one value is associated with this position')
	if scalar @array > 1;
    return $array[0];
}

=head2 add_position

 Title   : add_position
 Usage   : $position->add_position($position)
 Function: Add the Position to the Marker container.
           If you are using this method, you need to 
           add the Marker to the Map yourself
 Returns : none
 Args    : Position - Reference to Bio::Map::PositionI 

=cut

sub add_position{
   my ($self, $pos) = @_;
   $self->throw("Must give a Position") unless defined $pos;

   $self->throw("Must give a Bio::Map::PositionI, not [". ref($pos) ."]") 
       unless $pos->isa('Bio::Map::PositionI');

   my $map = $pos->map;
   $map->add_element($self) unless defined($self->map) && $self->map eq $map;

   push @{$self->{'_positions'}}, $pos;

}

=head2 positions

 Title   : positions
 Usage   : $mappable->positions([$pos1, $pos2, $pos3]); 
 Function: Add multiple Bio::Map::PositionI for a mappable element
           in a Map.
 Returns : boolean
 Args    : array ref of $map/value tuples or array ref of Positions

=cut

sub positions {
    my ($self, $arrayref) = @_;
    my ($map);
   $self->throw_not_implemented();
}

=head2 each_position

 Title   : each_position
 Usage   : my @positions = $position->each_position('mapname');
 Function: Retrieve a list of Positions 
 Returns : Array of L<Bio::Map::PositionI>
 Args    : none

=cut

sub each_position {
   my ($self,$mapname) = @_;
   $self->warn("Retrieving positions in a named map only is  ".
	       "not implemented. Getting all.") if $mapname;
   return @{$self->{'_positions'}}; 
}

=head2 purge_positions

 Title   : purge_positions
 Usage   : $marker->purge_positions
 Function: Remove all the position values stored for a Marker
 Returns : none
 Args    : [optional] only purge values for a given map

=cut

sub purge_positions{
    my ($self, $map) = @_;
    $self->warn("Retrieving positions in a named map only, not implemented ") if $map;
    $self->{'_positions'} = [];
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
   my %hash;
   foreach my $pos ($self->each_position) {
       $hash{$pos->map->unique_id} = 1; 
   }
   return keys %hash;
}

=head2 in_map

 Title   : in_map
 Usage   : if ( $position->in_map($map) ) {}
 Function: Tests if a position has values in a specific map
 Returns : boolean
 Args    : a map unique id OR Bio::Map::MapI

=cut

sub in_map{
   my ($self,$map) = @_;

   $self->throw("Need  an argument") unless $map;

   if (ref($map) && $map->isa('Bio::Map::MapI')) {
       foreach my $pos ($self->each_position) {
	   return 1 if $pos->map eq $map;
       }
   } else { # assuming a scalar
       foreach my $pos ($self->each_position) {
	   return 1 if $pos->map->unique_id eq $map;
       }
   }
   return 0;
}

=head2 Comparison methods

=cut

=head2 tuple

 Title   : tuple
 Usage   : ($me, $you) = $self->_tuple($compare)
 Function: Utility ethod to extract numbers and test for missing values.
 Returns : tuple values
 Args    : Bio::Map::MappableI or Bio::Map::PositionI

=cut

sub tuple {
    my ($self,$compare) = @_;
    my ($me, $you) = (-1, -1);

    $self->warn("Trying to compare [". $self->name. "] to nothing.") &&
	return ($me, $you) unless defined $compare;
    $self->warn("[". $self->name. "] has no position.") &&
	return ($me, $you) unless $self->position;

    $me  = $self->position->numeric;

    if( $compare->isa('Bio::Map::MappableI') ){
	$self->warn("[". $compare->name. "] has no position.") &&
	    return ($me, $you) unless $compare->position;

	$you = $compare->position->numeric;
	return ($me, $you);

    } elsif( $compare->isa('Bio::Map::PositionI') ) {

	$you = $compare->numeric;
	return ($me, $you);

    } else { 
	$self->warn("Can only run equals with Bio::Map::MappableI or ".
		    "Bio::Map::PositionI not [$compare]"); 
    }
    return ($me, $you);
}


=head2 equals

 Title   : equals
 Usage   : if( $mappable->equals($mapable2)) ...
 Function: Test if a position is equal to another position
 Returns : boolean
 Args    : Bio::Map::MappableI or Bio::Map::PositionI

=cut

sub equals {
    my ($self,$compare) = @_;

    my ($me, $you) = $self->tuple($compare);
    return 0 if $me == -1 or $you == -1 ;
    return $me == $you;
}

=head2 less_than

 Title   : less_than
 Usage   : if( $mappable->less_than($m2) ) ...
 Function: Tests if a position is less than another position
 Returns : boolean
 Args    : Bio::Map::MappableI  or Bio::Map::PositionI

=cut

sub less_than {
    my ($self,$compare) = @_;

    my ($me, $you) = $self->tuple($compare);
    return 0 if $me == -1 or $you == -1 ;
    return $me < $you;
}

=head2 greater_than

 Title   : greater_than
 Usage   : if( $mappable->greater_than($m2) ) ...
 Function: Tests if position is greater than another position
 Returns : boolean
 Args    : Bio::Map::MappableI or Bio::Map::PositionI

=cut

sub greater_than {
    my ($self,$compare) = @_;


    my ($me, $you) = $self->tuple($compare);
    return 0 if $me == -1 or $you == -1 ;
    return $me > $you;
}

1;
