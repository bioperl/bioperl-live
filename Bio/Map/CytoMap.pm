# $Id$
#
# BioPerl module for Bio::Map::CytoMap
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::CytoMap - A Bio::MapI compliant map implementation handling cytogenic bands 

=head1 SYNOPSIS

    use Bio::Map::CytoMap;
    my $map = new Bio::Map::CytoMap(-name => 'human1',
				      -species => $human);
    foreach my $marker ( @markers ) { # get a list of markers somewhere
	$map->add_element($marker);
    }

=head1 DESCRIPTION

This is the simple implementation of cytogenetic maps based on
L<Bio::Map::MapI>.  It handles the essential storage of name, species,
type, and units as well as in memory representation of the elements of
a map.

For CytoMaps type is hard coded to be 'cytogeneticmap' and
units are set to '' but can be set to something else.

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

=head1 AUTHOR - Heikki Lehvaslaiho

Email heikki@ebi.ac.uk

=head1 CONTRIBUTORS

Jason Stajich      jason@bioperl.org
Lincoln Stein      lstein@cshl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Map::CytoMap;
use vars qw(@ISA $MAPCOUNT);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::Map::SimpleMap;

@ISA = qw(Bio::Root::Root Bio::Map::SimpleMap);
BEGIN { $MAPCOUNT = 1; }

=head2 Modified methods 

All methods present in L<Bio::Map::SimpleMap> are implemted by this
class. Most of the methods are inherited from SimpleMap. The following
methods have been modified to refelect the needs of cytogenetic maps.

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Map::CytoMap();
 Function: Builds a new Bio::Map::CytoMap object
 Returns : Bio::Map::CytoMap
 Args    : -name    => name of map (string)
           -species => species for this map (Bio::Species) [optional]
           -elements=> elements to initialize with
                       (arrayref of Bio::Map::MappableI objects) [optional]

           -uid     => Unique Id
=cut

sub new {
    my($class,@args) = @_;

    my $self = $class->SUPER::new(@args);

    $self->{'_elements'} = [];
    $self->{'_name'}     = '';
    $self->{'_species'}  = '';
    $self->{'_units'}    = '';
    $self->{'_type'}    = 'cyto';
    $self->{'_uid'} = $MAPCOUNT++;
    my ($name, $type,$species, $units,
	$elements,$uid) = $self->_rearrange([qw(NAME TYPE
						SPECIES UNITS
						ELEMENTS UID)], @args);
    defined $name     && $self->name($name);
    defined $species  && $self->species($species);
    defined $units    && $self->units($units);
    defined $type     && $self->type($type);
    defined $uid      && $self->unique_id($uid);

    if( $elements && ref($elements) =~ /array/ ) {
	foreach my $item ( @$elements ) {
	    $self->add_element($item);
	}
    }
    return $self;
}

=head2 type

 Title   : type
 Usage   : my $type = $map->type
 Function: Get hard-coded  Map type
 Returns : String coding map type
 Args    : 

=cut

sub type {
   my ($self) = @_;
   return $self->{'_type'};
}


=head2 length

 Title   : length
 Usage   : my $length = $map->length();
 Function: Retrieves the length of the map,
 Returns : undef since length is not calculatable for 
           cytogenetic maps
 Args    : none

=cut

sub length{
   my ($self,@args) = @_;
   return undef;
}

=head2 Methods inherited from L<Bio::Map::SimpleMap>

=cut

=head2 species

 Title   : species
 Usage   : my $species = $map->species;
 Function: Get/Set Species for a map
 Returns : Bio::Species object or string
 Args    : (optional) Bio::Species or string

=cut

=head2 units

 Title   : units
 Usage   : $map->units('cM');
 Function: Get/Set units for a map
 Returns : units for a map
 Args    : units for a map (string)

=cut

=head2 name

 Title   : name
 Usage   : my $name = $map->name
 Function: Get/Set Map name
 Returns : Map name
 Args    : (optional) string

=cut

=head2 unique_id

 Title   : unique_id
 Usage   : my $id = $map->unique_id;
 Function: Get/Set the unique ID for this map
 Returns : a unique identifier
 Args    : [optional] new identifier to set

=cut

=head2 each_element

 Title   : each_element
 Usage   : my @elements = $map->each_element;
 Function: Retrieves all the elements in a map
           unordered
 Returns : Array of Bio::Map::MappableI objects
 Args    : none


=cut

=head2 New methods

=cut


1;
