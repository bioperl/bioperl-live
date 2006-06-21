# $Id$
#
# BioPerl module for Bio::Map::Marker
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
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
				   );

  # an other way of creating a marker with unique position in a map:
  $marker2 = new Bio::Map::Marker (-name=>'hypervariable2');
  $map1->add_element($marker2);
  $marker2->position(100);

  # position method is a short cut for get/setting unique positions
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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Chad Matsalla

Email bioinformatics1@dieselwurks.com

=head1 CONTRIBUTORS

Heikki Lehvaslaiho heikki-at-bioperl-dot-org
Lincoln Stein      lstein@cshl.org
Jason Stajich      jason@bioperl.org
Sendu Bala         bix@sendu.me.uk

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
           -map => the default map for this marker, a Bio::Map::MapI
           -position => map position for this marker, [optional]
                Bio::Map::PositionI-inherited obj, no default)
           -positions => array ref of map positions, as above
           
           position and positions can also take as values anything the
           corresponding methods can take

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
    
    $map       && $self->map($map);
    $position  && $self->position($position); 
    $positions && $self->positions($positions);
 
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
           Alternatively, you can add the marker to the map first 
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
            }
            else {
                $map = $self->map;
            }
            
            $self->throw("You need to add a marker to a map before you can set positions without explicit map!" ) unless $map;
            $self->throw("Position better be scalar, not [$pos=". ref($pos)  ."]") unless ref($pos) eq 'SCALAR' || ref($pos) eq ''; 
            
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
    $self->warn('More than one value is associated with this position') if scalar @array > 1;
    return $array[0];
}

=head2 add_position

 Title   : add_position
 Usage   : $marker->add_position($position);
 Function: Add a Position to this marker
 Returns : n/a
 Args    : Bio::Map::PositionI
           * or *
           {
             Bio::Map::MapI
             * AND *
             scalar
           }
           * or *
           scalar, but only if the marker has been added to a map

=cut
sub add_position  {
    my ($self, $pos, $pos_actual) = @_;
    $self->throw("Must give a Position") unless defined $pos;
    
    if (ref($pos) && $pos->isa('Bio::Map::PositionI')) {
        my $map = $pos->map;
        $map->add_element($self) unless defined($self->map) && $self->map eq $map;
    }
    elsif ($pos_actual && ref($pos_actual) eq '' && ref($pos) && $pos->isa('Bio::Map::MapI')) {
        my $map = $pos;
        $map->add_element($self) unless defined($self->map) && $self->map eq $map;
        $pos = $self->_make_pos($map, $pos_actual);
    }
    elsif (ref($pos) eq '' && $pos =~ /^\d+$/) {
        my $newpos = $self->get_position_object;
        my $map = $self->map();
        $map or $self->throw("Must already have a default map set");
        $pos = $self->_make_pos($map, $pos);
    }
    else {
        $self->throw("Must give a Bio::Map::PositionI or an int, not [". ref($pos) ."]");
    }
    
    push @{$self->{'_positions'}}, $pos;
}

=head2 positions

 Title   : positions
 Usage   : $marker->positions([$pos1, $pos2, $pos3]);
 Function: Add multiple Bio::Map::PositionI to this marker
 Returns : n/a
 Args    : array ref of $map/value tuples or array ref of Bio::Map::PositionI

=cut
sub positions {
    my ($self, $args_ref) = @_;
    
    foreach my $arg (@{$args_ref}) {
        if (ref($arg) eq 'ARRAY') {
            $self->add_position(@{$arg});
        }
        else {
            $self->add_position($arg);
        }
    }
}

=head2 each_position

 Title   : each_position
 Usage   : my @positions = $marker->each_position('map_unique_id');
 Function: Retrieve a list of Positions
 Returns : array of TFBS::Map::Position
 Args    : nothing for all
           Bio::Map::MapI OR unique_id for positions on the given map

=cut
sub each_position {
    my ($self, $map) = @_;
    
    if ($map) {
        my @positions;
        if (ref($map) && $map->isa('Bio::Map::MapI')) {
            foreach my $pos (@{$self->{'_positions'}}) {
                if ($pos->map && $pos->map eq $map) {
                    push(@positions, $pos);
                }
            }
        }
        else {
            foreach my $pos (@{$self->{'_positions'}}) {
                if ($pos->map && $pos->map()->unique_id eq $map) {
                    push(@positions, $pos);
                }
            }
        }
        
        return @positions;
    }
    
    return @{$self->{'_positions'}};
}

=head2 purge_positions

 Title   : purge_positions
 Usage   : $marker->purge_positions
 Function: remove all the position values stored for a Marker
 Returns : n/a
 Args    : nothing for all
           Bio::Map::MapI OR unique_id to only purge positions on the given map

=cut
sub purge_positions  {
    my ($self, $map) = @_;
    
    if ($map) {
        my @dontwant = $self->each_position($map);
        my %dontwant;
        foreach my $pos (@dontwant) {
            $dontwant{$pos} = 1;
        }
        
        my @wanted;
        foreach my $pos (@{$self->{'_positions'}}) {
            exists $dontwant{$pos} and next;
            push(@wanted, $pos);
        }
        
        $self->{'_positions'} = \@wanted;
        return;
    }
    
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
 Usage   : Do Not Use!
 Function: tuple was supposed to be a private method; this method no longer
           does anything
 Returns : warning
 Args    : none

=cut

sub tuple {
    my $self = shift;
    $self->warn("The tuple method was supposed to be a private method, don't call it!");
}

=head2 equals

 Title   : equals
 Usage   : if ($marker->equals($other_marker)) {...}
           my @equal_positions = $marker->equals($comparison_marker);
 Function: Finds the positions in this marker that are equal to any
           comparison positions, optionally only considering a particular map.
 Returns : array of those Bio::Map::PositionI objects in self that satisfy the
           above criteria
 Args    : Bio::Map::MappableI or Bio::Map::PositionI [REQUIRED]
           Bio::Map::MapI (to only consider positions on this map)

=cut
sub equals {
    my $self = shift;
    return $self->_compare('equals', @_);
}

=head2 less_than

 Title   : less_than
 Usage   : if ($marker->less_than($comparison_marker)) {...}
           my @lower_positions = $marker->less_than($comparison_marker);
 Function: Finds the positions in this marker that are less than all
           comparison positions, optionally only considering a particular map.
 Returns : array of those Bio::Map::PositionI objects in self that satisfy the
           above criteria
 Args    : Bio::Map::MappableI or Bio::Map::PositionI [REQUIRED]
           Bio::Map::MapI (to only consider positions on this map)

=cut
sub less_than {
    my ($self, $compare, $map) = @_;
    my ($mine, $yours) = $self->_pre_compare($compare, $map);
    
    (@{$mine} > 0 && @{$yours} > 0) or return ();
    my $least = ${$yours}[0]->start;
    
    my @less;
    foreach my $pos (@{$mine}) {
        $pos->end < $least or last;
        push(@less, $pos);
    }
    
    return @less;
}

=head2 greater_than

 Title   : greater_than
 Usage   : if ($marker->greater_than($comparison_marker)) {...}
           my @higher_positions = $marker->greater_than($comparison_marker);
 Function: Finds the positions in this marker that are greater than all
           comparison positions, optionally only considering a particular map.
 Returns : array of those Bio::Map::PositionI objects in self that satisfy the
           above criteria
 Args    : Bio::Map::MappableI or Bio::Map::PositionI [REQUIRED]
           Bio::Map::MapI (to only consider positions on this map)

=cut
sub greater_than {
    my ($self, $compare, $map) = @_;
    my ($mine, $yours) = $self->_pre_compare($compare, $map);
    
    (@{$mine} > 0 && @{$yours} > 0) or return ();
    
    my $greatest = ${$yours}[-1]->end;
    
    my @more;
    foreach my $pos (@{$mine}) {
        $pos->start > $greatest or next;
        push(@more, $pos);
    }
    
    return @more;
}

=head2 overlaps

 Title   : overlaps
 Usage   : if ($marker->overlaps($other_marker)) {...}
           my @overlapping_positions = $marker->overlaps($comparison_marker);
 Function: Finds the positions in this marker that overlap with any
           comparison positions, optionally only considering a particular map.
 Returns : array of those Bio::Map::PositionI objects in self that satisfy the
           above criteria
 Args    : Bio::Map::MappableI or Bio::Map::PositionI [REQUIRED]
           Bio::Map::MapI (to only consider positions on this map)

=cut
sub overlaps {
    my $self = shift;
    return $self->_compare('overlaps', @_);
}

=head2 contains

 Title   : contains
 Usage   : if ($marker->overlaps($other_marker)) {...}
           my @container_positions = $marker->contains($comparison_marker);
 Function: Finds the positions in this marker that contain any comparison
           positions, optionally only considering a particular map.
 Returns : array of those Bio::Map::PositionI objects in self that satisfy the
           above criteria
 Args    : Bio::Map::MappableI or Bio::Map::PositionI [REQUIRED]
           Bio::Map::MapI (to only consider positions on this map)

=cut
sub contains {
    my $self = shift;
    return $self->_compare('contains', @_);
}

=head2 _compare

 Title   : _compare
 Usage   : my @positions = $marker->_compare($method, $comparison_marker, $map);
 Function: do a RangeI comparison
 Returns : array of Bio::Map::PositionI
 Args    : string, a RangeI comparison method name,
           AND a Bio::Map::MappableI [REQUIRED]
           Bio::Map::MapI (to only consider positions on this map)

=cut

sub _compare {
    my ($self, $method, $compare, $map) = @_;
    my ($mine, $yours) = $self->_pre_compare($compare, $map);
    
    (@{$mine} > 0 && @{$yours} > 0) or return ();
    
    my @ok;
    foreach my $my_pos (@{$mine}) {
        foreach my $your_pos (@{$yours}) {
            if ($my_pos->$method($your_pos)) {
                push(@ok, $my_pos);
                last;
            }
        }
    }
    
    return @ok;
}

=head2 _pre_compare

 Title   : _pre_compare
 Usage   : my @positions = $marker->_compare($method, $comparison_marker, $map);
 Function: test for missing values and discover if we have multiple positions so
           that we can do some kind of comparison later
 Returns : (\@sorted_self_positions, \@sorted_compare_positions)
 Args    : Bio::Map::MappableI OR Bio::Map::PositionI

=cut

sub _pre_compare {
    my ($self, $compare, $map) = @_;
    my (@mine, @yours);
    
    $self->warn("Trying to compare [". $self->name. "] to nothing or scalar; need object.") && return (\@mine, \@yours) unless (defined $compare && ref($compare));
    
    @mine = $self->each_position($map);
    my $on_this_map = $map ? ' (on the supplied map)' : '';
    @mine > 0 or ($self->warn("[". $self->name. "] has no positions$on_this_map.") && return (\@mine, \@yours));
    
    if ($compare->isa('Bio::Map::PositionI')) {
        push(@yours, $compare);
    }
    elsif ($compare->isa('Bio::Map::MappableI')) {
        @yours = $compare->each_position($map);
        @yours > 0 or ($self->warn("[". $compare->name. "] has no positions$on_this_map.") && return (\@mine, \@yours));
    }
    else {
        $self->warn("Can only run a comparison with a Bio::Map::MappableI or Bio::Map::PositionI, not [$compare]");
        return (\@mine, \@yours);
    }
    
    @mine = sort { $a->numeric <=> $b->numeric } @mine;
    @yours = sort { $a->numeric <=> $b->numeric } @yours;
    return (\@mine, \@yours)
}

=head2 _make_pos

 Title   : _make_pos
 Usage   : my $pos = $marker->_make_pos($map, $value);
 Function: make a new Position object on a given map at the given position
 Returns : L<Bio::Map::Position>
 Args    : L<Bio::Map::MapI>, scalar (the position on the map)

=cut

sub _make_pos {
    my ($self, $map, $value) = @_;
    my $pos = $self->get_position_object;
    $pos->map($map);
    $pos->value($value);
    $pos->marker($self);
    
    return $pos;
}

1;
