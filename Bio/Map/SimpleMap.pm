# $Id$
#
# BioPerl module for Bio::Map::SimpleMap
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::SimpleMap - A MapI implementation handling the basics of a Map 

=head1 SYNOPSIS

    use Bio::Map::SimpleMap;
    my $map = new Bio::Map::SimpleMap(-name => 'genethon',
				      -type => 'Genetic',
				      -units=> 'cM',
				      -species => $human);
    foreach my $marker ( @markers ) { # get a list of markers somewhere
	$map->add_element($marker);
    }
	
=head1 DESCRIPTION

This is the basic implementation of a Bio::Map::MapI.  It handles the
essential storage of name, species, type, and units as well as in
memory representation of the elements of a map.

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

=head1 CONTRIBUTORS

Heikki Lehvaslaiho heikki@ebi.ac.uk
Lincoln Stein      lstein@cshl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Map::SimpleMap;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::Map::MapI;

@ISA = qw(Bio::Root::Root Bio::Map::MapI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Map::SimpleMap();
 Function: Builds a new Bio::Map::SimpleMap object 
 Returns : Bio::Map::SimpleMap
 Args    : -name    => name of map (string)
           -species => species for this map (Bio::Species) [optional]
           -units   => map units (string)
           -elements=> elements to initialize with
                       (arrayref of Bio::Map::MappableI objects) [optional]

=cut

sub new {
  my($class,@args) = @_;  
  
  my $self = $class->SUPER::new(@args);
  
  $self->{'_elements'} = [];
  $self->{'_name'}     = '';
  $self->{'_species'}  = '';
  $self->{'_units'}    = '';
  $self->{'_type'}    = '';

  my ($name, $type,$species, $units,
      $elements) = $self->_rearrange([qw(NAME TYPE 
					 SPECIES UNITS 
					 ELEMENTS)], @args);
  defined $name     && $self->name($name);
  defined $species  && $self->species($species);
  defined $units    && $self->units($units);
  defined $type     && $self->type($type);
  if( $elements && ref($elements) =~ /array/ ) {
      foreach my $item ( @$elements ) {
	  $self->add_element($item);
      }
  } 
  return $self;
}

=head2 species

 Title   : species
 Usage   : my $species = $map->species;
 Function: Get/Set Species for a map
 Returns : Bio::Species object
 Args    : (optional) Bio::Species

=cut

sub species{
   my ($self,$value) = @_;
   if( defined $value ) {
       if( ! $value->isa('Bio::Species') ) {
	   $self->warn("Trying to set species to invalid type. $value is not a Bio::Species ");
       } else { 
	   $self->{'_species'} = $value;
       }
   }
   return $self->{'_species'};
}

=head2 units

 Title   : units
 Usage   : $map->units('cM');
 Function: Get/Set units for a map
 Returns : units for a map
 Args    : units for a map (string)

=cut

sub units{
   my ($self,$value) = @_;
   if( defined $value ) {
       $self->{'_units'} = $value;
   }
   return $self->{'_units'};
}

=head2 type

 Title   : type
 Usage   : my $type = $map->type
 Function: Get/Set Map type
 Returns : String coding map type
 Args    : (optional) string

=cut

sub type {
   my ($self,$value) = @_;
   # this may be hardcoded/overriden by subclasses

   if( defined $value ) {
       $self->{'_type'} = $value;
   }
   return $self->{'_type'};
}


=head2 name

 Title   : name
 Usage   : my $name = $map->name
 Function: Get/Set Map name
 Returns : Map name
 Args    : (optional) string

=cut

sub name {
   my ($self,$value) = @_;
   if( defined $value ) {
       $self->{'_name'} = $value;
   }
   return $self->{'_name'};
}

=head2 length

 Title   : length
 Usage   : my $length = $map->length();
 Function: Retrieves the length of the map, 
           It is possible for the length to be unknown
           for maps such as Restriction Enzyme, will return undef
           in that case
 Returns : integer representing length of map in current units
           will return undef if length is not calculateable
 Args    : none

=cut

sub length{
   my ($self,@args) = @_;
   if( defined $self->{'_lastelement'} ) {
       my @p = $self->{'_lastlement'}->position->each_position;
       return pop @p; # let's just get the last position for the last element
                      # for now
   } else { return 0; } 
}


=head2 add_element

 Title   : add_element
 Usage   : $map->add_element($marker)
 Function: Add a Bio::Map::MappableI object to the Map
 Returns : none
 Args    : Bio::Map::MappableI object

=cut

sub add_element{
    my ($self,$mapelement) = @_;
    return unless ( defined $mapelement);

    push @{$self->{'_elements'}}, $mapelement;
    if( !defined $self->{'_lastelement'} || 
	$mapelement->greater_than($self->{'_lastelement'}) ) {
	$self->{'_lastelement'} = $mapelement;    
    }
}

=head2 each_element

 Title   : each_element
 Usage   : my @elements = $map->each_element;
 Function: Retrieves all the elements in a map
           unordered
 Returns : Array of Bio::Map::MappableI objects
 Args    : none


=cut

sub each_element{
   my ($self) = @_;
   return @{$self->{'_elements'}};
}

1;
