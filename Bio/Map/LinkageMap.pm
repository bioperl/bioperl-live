# BioPerl module for Bio::Map::LinkageMap
#
# Cared for by Chad Matsalla <bioinformatics1@dieselwurks.com>
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::LinkageMap - A representation of a genetic linkage map.

=head1 SYNOPSIS

    use Bio::Map::LinkageMap;
	# create a new map
    my $map = new Bio::Map::LinkageMap(-name => 'Chads Superterriffic Map',
                                      -type => 'Linkage',
                                      -units=> 'cM');
	# create the location of a marker for that map
    my $position = new Bio::Map::LinkagePosition( -positions => 1,
                -distance => "22.3");
	# create a marker and place it at that position
    my $marker = new Bio::Map::Marker::Microsatellite(
			-name => 'SuuuperMarker',
			-position => $position);
	# place that marker on that map
    $map->add_element($marker);

	# done!

=head1 DESCRIPTION

This object describes the basic functionality of a genetic linkage map in
Bioperl. Each 'position' can have one or more markers that map some number of
units from the markers at the previous position.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org            - General discussion
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

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Map::LinkageMap;
use vars qw(@ISA);
use strict;
use Bio::Root::RootI;
use Carp;

@ISA = qw(Bio::Root::RootI Bio::Map::MapI);

=head2 new

 Title   : new
 Usage   : my $linkage_map = new Bio::Map::LinkageMap();
 Function: Builds a new Bio::Map::LinkageMap object
 Returns : Bio::Map::LinkageMap
 Args    : -name    => the name of the map (string) [optional]
	   -type    => the type of this map (string, defaults to Linkage) [optional]
           -species => species for this map (Bio::Species) [optional]
           -units   => the map units (string, defaults to cM) [optional]
           -elements=> elements to initialize with
                       (arrayref of Bio::Map::MappableI objects) [optional]

=cut

sub new {
  my ($class,@args) = @_;
	my $self = {};
  $self->{'_elements'} = [];
  $self->{'_name'}     = '';
  $self->{'_species'}  = '';
  $self->{'_units'}    = '';

  my $self = $class->SUPER::new(@args);
  my ($name, $species, $units,
      $elements, $type) = $self->_rearrange([qw(NAME SPECIES UNITS ELEMENTS TYPE)], @args);
  $name     && $self->name($name);
  $species  && $self->species($species);
	if ($units) { $self->units($units); }
	else { $self->units('cM'); }
	if ($type) { $self->type($type); }
	else { $self->type('Linkage'); }
  if( $elements && ref($elements) =~ /array/ ) {
      foreach my $item ( @$elements ) {
          $self->add_element($item);
      }
  }
  return $self;
}


=head2 species($new_species)

 Title   : species($new_species)
 Usage   : my $species = $map->species($new_species) _or_
	my $species = $map->species()
 Function: Get/Set Species for a map
 Returns : A Bio::Species object representing the current species of this
	LinkageMap.
 Args    : If provided, the species of this LinkageMap will be set to
	$new_species.
	

=cut

sub species {
	my ($self,$species) = @_;
	if ($species) {
		$self->{'_species'} = $species;
	}
	return $self->{'_species'};
}

=head2 units($new_units)

 Title   : units($new_units)
 Usage   : $map->units($new_units) _or_
	$map->units()
 Function: Get/Set units for a map
 Returns : A scalar representing the units for this LinkageMap
 Args    : If provided, the units for this LinkageMap will be set to
	$new_units.

=cut

sub units {
	my ($self,$units) = @_;
	if (!$units) {
		return $self->{'_units'};
	}
	$self->{'_units'} = $units;
}

=head2 type($new_type)

 Title   : type($new_type)
 Usage   : my $type = $map->type($new_type) _or_
	my $type = $map->type()
 Function: Get/Set Map type
 Returns : A scalar representing the current map type.
 Args    : If provided, the current map type will be set to $new_type.

=cut

sub type {
	my ($self,$type) = @_;
	if (!$type) {
		$self->{'_type'} = $type;
	}
	return $self->{'_type'};
}

=head2 length()

 Title   : length()
 Usage   : my $length = $map->length();
 Function: Retrieves the length of the map. In the case of a LinkageMap, the
	length is the sum of all marker distances.
 Returns : An integer representing the length of this LinkageMap. Will return
	undef if length is not calculateable
 Args    : None.


=cut

sub length {
	my ($self) = @_;
	my $total_distance;
	foreach (@{$self->{'_elements'}}) {
		if ($_) {
			$total_distance += $_->position()->distance();
		}
	}
	return $total_distance;
}

=head2 name($new_name)

 Title   : name($new_name)
 Usage   : my $name = $map->name($new_name) _or_
	my $length = $map->name()
 Function: Get/set the name of the map.
 Returns : The current name of the map.
 Args    : If provided, the name of the map is set to $new_name.

=cut

sub name {
	my ($self,$name) = @_;
	if ($name) {
		$self->{'_name'} = $name;
	}
	return $self->{'_name'};
}

=head2 add_element($marker)

 Title   : add_element($marker)
 Usage   : $map->add_element($marker)
 Function: Add a Bio::Map::MappableI object to the Map
 Returns : none
 Args    : Bio::Map::MappableI object
 Notes   : It is strongly recommended that you use a Bio::Map::LinkagePosition
	as the position in any Bio::Map::Mappable that you create to place on
	this map. Using some other Bio::Map::Position might work but might be
	unpredictable.

=cut

sub add_element {
	my ($self,$marker) = @_;
	my $o_position = $marker->position();
	if (ref($o_position) != /Linkage/) {
		$self->warn("You really should use a Linkage Position for this object. This insures that there is only one position. Trying anyway...");
	}
	my $position = $o_position->position();
	if ($self->{'_elements'}[$position]) {
		$self->warn("Replacing the marker in position $position because in a linkage map the position is a key.");
	}	
	$self->{'_elements'}[$position] = $marker;
}

=head2 each_element

 Title   : each_element
 Usage   : my @elements = $map->each_element;
 Function: Retrieves all the elements in a map
           _ordered_.
 Returns : An array containing MappableI objects.
 Args    : None.
 Notes   : This is a useless concept in the context of a linkage map but is
	included if you want a list of all of the marker names on the map.

=cut

sub each_element {
	my ($self) = @_;
	return @{$self->{'_elements'}};
}

1;
