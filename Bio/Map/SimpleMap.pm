#
# BioPerl module for Bio::Map::SimpleMap
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::SimpleMap - A MapI implementation handling the basics of a Map 

=head1 SYNOPSIS

    use Bio::Map::SimpleMap;

    my $map = Bio::Map::SimpleMap->new(-name => 'genethon',
				      -type => 'Genetic',
				      -units=> 'cM',
				      -species => $human);

    foreach my $marker ( @markers ) { # get a list of markers somewhere
		$map->add_element($marker);
    }

    foreach my $marker ($map->get_elements) {
        # do something with this Bio::Map::MappableI
    }

=head1 DESCRIPTION

This is the basic implementation of a Bio::Map::MapI. It handles the
essential storage of name, species, type, and units.

It knows which map elements (mappables) belong to it, and their
position.

Subclasses might need to redefine or hardcode type(), length() and
units().

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 CONTRIBUTORS

Heikki Lehvaslaiho heikki-at-bioperl-dot-org
Lincoln Stein      lstein@cshl.org
Sendu Bala         bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::SimpleMap;
use vars qw($MAPCOUNT);
use strict;


use base qw(Bio::Root::Root Bio::Map::MapI);
BEGIN { $MAPCOUNT = 1; }

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Map::SimpleMap->new();
 Function: Builds a new Bio::Map::SimpleMap object
 Returns : Bio::Map::SimpleMap
 Args    : -name    => name of map (string)
           -species => species for this map (Bio::Species) [optional]
           -units   => map units (string)
           -uid     => Unique Id [defaults to a unique integer]

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  $self->{'_name'}     = '';
  $self->{'_species'}  = '';
  $self->{'_units'}    = '';
  $self->{'_type'}    = '';
  $self->{'_uid'} = $MAPCOUNT++;
  my ($name, $type,$species, $units,$uid) = $self->_rearrange([qw(NAME TYPE
					      SPECIES UNITS
					      UID)], @args);
  defined $name     && $self->name($name);
  defined $species  && $self->species($species);
  defined $units    && $self->units($units);
  defined $type     && $self->type($type);
  defined $uid      && $self->unique_id($uid);
    
  return $self;
}

=head2 species

 Title   : species
 Usage   : my $species = $map->species;
 Function: Get/Set Species for a map
 Returns : Bio::Taxon object or string
 Args    : (optional) Bio::Taxon or string

=cut

sub species{
   my ($self,$value) = @_;
   if( defined $value ) {
       $self->{'_species'} = $value;
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
 Function: Retrieves the length of the map.
           It is possible for the length to be unknown for maps such as
           Restriction Enzyme, will return 0 in that case.
 Returns : integer representing length of map in current units
           will return 0 if length is not calculateable
 Args    : none

=cut

sub length {
	my $self = shift;
    
	my $len = 0;
    foreach my $element ($self->get_elements) {
        foreach my $pos ($element->get_positions($self)) {
            if ($pos->value) {
                $len = $pos->end if $pos->end > $len;
            }
        }
	}
    
	return $len;
}

=head2 unique_id

 Title   : unique_id
 Usage   : my $id = $map->unique_id;
 Function: Get/Set the unique ID for this map
 Returns : a unique identifier
 Args    : [optional] new identifier to set

=cut

sub unique_id {
   my ($self,$id) = @_;
   if( defined $id ) {
       $self->{'_uid'} = $id;
   }
   return $self->{'_uid'};
}

=head2 add_element

 Title   : add_element
 Usage   : $map->add_element($element)
 Function: Tell a Bio::Map::MappableI object its default Map is this one; same
           as calling $element->default_map($map).

		   *** does not actually add the element to this map! ***

 Returns : none
 Args    : Bio::Map::MappableI object
 Status  : Deprecated, will be removed in next version

=cut

sub add_element {
    my ($self, $element) = @_;
    return unless $element;
	
    $self->throw("This is not a Bio::Map::MappableI object but a [$element]")
	unless $element->isa('Bio::Map::MappableI');
	
    $element->default_map($self);
}

=head2 get_elements

 Title   : get_elements
 Usage   : my @elements = $map->get_elements;
 Function: Retrieves all the elements on a map (unordered unless all elements
           have just 1 position on the map, in which case sorted)
 Returns : Array of Map elements (L<Bio::Map::MappableI>)
 Args    : none

=cut

sub get_elements {
	my $self = shift;
	
    my @elements = $self->SUPER::get_elements;
    
	# for backward compatability with MapIO tests, and for 'niceness', when
	# there is only 1 position per element we will return the elements in
	# order, as long as the positions have values set
    my $only_1 = 1;
    foreach my $element (@elements) {
        my @positions = $element->get_positions($self);
        if (@positions > 1 || (@positions == 1 && ! $positions[0]->value)) {
            $only_1 = 0;
        }
    }
	if ($only_1) {
		@elements = map { $_->[1] }
                    sort { $a->[0] <=> $b->[0] }
                    map { [${[$_->get_positions($self)]}[0]->sortable, $_] }
                    @elements;
	}
	
	return @elements;
}

=head2 each_element

 Title   : each_element
 Function: Synonym of the get_elements() method.
 Status  : deprecated, will be removed in the next version

=cut

*each_element = \&get_elements;

=head2 purge_element

 Title   : purge_element
 Usage   : $map->purge_element($element)
 Function: Purge an element from the map.
 Returns : none
 Args    : Bio::Map::MappableI object

=cut

sub purge_element {
    my ($self, $element) = @_;
    $self->throw("Must supply an argument") unless $element;
    $self->throw("This is [$element], not an object") unless ref($element);
    $self->throw("This is [$element], not a Bio::Map::MappableI object") unless $element->isa('Bio::Map::MappableI');
	
	$self->purge_positions($element);
}

=head2 annotation

 Title   : annotation
 Usage   : $map->annotation($an_col);
           my $an_col = $map->annotation();
 Function: Get the annotation collection (see Bio::AnnotationCollectionI)
           for this annotatable object.
 Returns : a Bio::AnnotationCollectionI implementing object, or undef
 Args    : none to get, OR
           a Bio::AnnotationCollectionI implementing object to set

=cut

sub annotation {
    my $self = shift;
    if (@_) { $self->{_annotation} = shift }
    return $self->{_annotation} || return;
}

1;
