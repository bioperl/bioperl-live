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
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Chad Matsalla

Email bioinformatics1@dieselwurks.com

=head1 CONTRIBUTORS

Lincoln Stein       lstein@cshl.org
Heikki Lehvaslaiho  heikki@ebi.ac.uk
Jason Stajich       jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Map::LinkageMap;
use vars qw(@ISA);
use strict;
use Bio::Map::SimpleMap;

@ISA = qw(Bio::Map::SimpleMap);

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

          -uid      => Unique ID of this map
=cut

# new provided by SimpleMap



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
	    $total_distance += ($_->position()->each_position_value($self))[0];
	}
    }
    return $total_distance;
}

=head2 add_element($marker)

 Title   : add_element($marker)
 Usage   : $map->add_element($marker)
 Function: Add a Bio::Map::MappableI object to the Map
 Returns : none
 Args    : Bio::Map::MappableI object
 Notes   : It is strongly recommended that you use a
	   Bio::Map::LinkagePosition as the position in any
	   Bio::Map::Mappable that you create to place on this
	   map. Using some other Bio::Map::Position might work but might
	   be unpredictable.
           N.B. I've added Bio::Map::OrderedPosition which should achieve
                similar things from LinkagePosition and will work for
                RH markers too.
=cut

#'
sub _add_element {
    my ($self,$marker) = @_;

    my $o_position = $marker->position();

    $self->debug( "marker position is ". $marker->position());
#     print("add_element: \$o_position is $o_position\n");
#     print("add_element: \$marker is $marker\n");

    my $position;
    unless ( $o_position->isa('Bio::Map::LinkagePosition') ||
	     $o_position->isa('Bio::Map::OrderedPosition')
	     ) {
	$self->warn("You really should use a Linkage Position for this object. This insures that there is only one position. Trying anyway...");	
	my @p = ( $o_position->each_position_value($self));
	$position = shift @p;
	if( ! defined $position ) {
	    $self->throw("This marker ($marker) does not have a position in this map ($self)");
	}
    } else {
	$position = $o_position->order;
    }

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

=head2 implemented by Bio::Map::SimpleMap

=cut

=head2 name($new_name)

 Title   : name($new_name)
 Usage   : my $name = $map->name($new_name) _or_
	   my $length = $map->name()
 Function: Get/set the name of the map.
 Returns : The current name of the map.
 Args    : If provided, the name of the map is set to $new_name.

=head2 species

 Title   : species
 Usage   : my $species = $map->species;
 Function: Get/Set Species for a map
 Returns : Bio::Species object
 Args    : (optional) Bio::Species


=head2 units

 Title   : units
 Usage   : $map->units('cM');
 Function: Get/Set units for a map
 Returns : units for a map
 Args    : units for a map (string)


=head2 type

 Title   : type
 Usage   : my $type = $map->type
 Function: Get/Set Map type
 Returns : String coding map type
 Args    : (optional) string

=head2 unique_id

 Title   : unique_id
 Usage   : my $id = $map->unique_id;
 Function: Get/Set the unique ID for this map
 Returns : a unique identifier
 Args    : [optional] new identifier to set

=cut

1;
