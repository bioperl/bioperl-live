# $Id$
#
# BioPerl module for Bio::Map::MapI
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::MapI - Interface for describing Map objects in bioperl 

=head1 SYNOPSIS

    # get a MapI somehowe
    my $name   = $map->name();     # string
    my $length = $map->length();   # integer
    my $species= $map->species;    # Bio::Species
    my $type   = $map->type();     # genetic/sts/rh/

=head1 DESCRIPTION

This object describes the basic functionality of a Map in bioperl.
Maps are anything from Genetic Map to Sequence Map to and Assembly Map
to Restriction Enzyme to FPC.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 CONTRIBUTORS

Lincoln Stein, lstein@cshl.org
Heikki Lehvaslaiho, heikki-at-bioperl-dot-org
Sendu Bala, bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::MapI;
use vars qw(@ISA);
use strict;
use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);


=head2 add_position

 Title   : add_position
 Usage   : $map->add_position($position);
 Function: Add a position to this map. This is how a mappable element is placed
           on the map. Notifies the position it is on this map - same as calling
           $position->map($map);.
 Returns : n/a
 Args    : L<Bio::Map::PositionI> object

=cut

sub add_position {
    # this is fundamental to coordination of Positions and Maps, so is
    # implemented at the interface level
    my ($self, $pos) = @_;
    $self->throw("Must supply an argument") unless $pos;
    $self->throw("This is [$pos], not an object") unless ref($pos);
    $self->throw("This is [$pos], not a Bio::Map::PositionI object") unless $pos->isa('Bio::Map::PositionI');
    
    $self->{'_positions'}->{$pos} = $pos;
    
    my $old_map = $pos->map;
    unless ($old_map && $old_map eq $self) {
        $pos->map($self);
    }
}

=head2 purge_positions

 Title   : purge_positions
 Usage   : $map->purge_position();
 Function: Remove all positions from this map. Notifies the positions they are
           no longer on this map.
 Returns : n/a
 Args    : none to remove all positions, OR
           L<Bio::Map::PositionI> object to remove just that Position, OR
		   L<Bio::Map::MappableI> object to remove only those positions of the
           given mappable

=cut

sub purge_positions {
    # this is fundamental to coordination of Positions and Maps, so is
    # implemented at the interface level
    my ($self, $thing) = @_;
    
	my @positions;
	if ($thing) {
		$self->throw("This is [$thing], not an object") unless ref($thing);
		$self->throw("This is [$thing], not a Bio::Map::PositionI or Bio::Map::MappableI object")
		unless ($thing->isa('Bio::Map::PositionI') || $thing->isa('Bio::Map::MappableI'));
		
		$thing->isa('Bio::Map::PositionI') ? push(@positions, $thing) : (@positions = $thing->each_position($self));
	}
	else {
		@positions = $self->each_position;
	}
	
	foreach my $pos (@positions) {
		defined $self->{'_positions'}->{$pos} or next;
		delete $self->{'_positions'}->{$pos};
		
        my $existing_map = $pos->map;
		if ($existing_map && $existing_map eq $self) {
			$pos->purge_map;
		}
	}
}

=head2 each_position

 Title   : each_position
 Usage   : my @positions = $map->each_position()
 Function: get all the Positions on this map (sorted)
 Returns : array of L<Bio::Map::PositionI> objects
 Args    : none for all, OR
           L<Bio::Map::MappableI> object for positions of the given mappable

=cut

sub each_position {
    # Positions were added with a method implemented in the interface, so this
    # accessor is implemented at the interface level as well
    my ($self, $element) = @_;
	
    if (defined $self->{'_positions'}) {
        my @positions;
		if ($element) {
			$self->throw("This is [$element], not an object") unless ref($element);
			$self->throw("This is [$element], not a Bio::Map::MappableI object") unless $element->isa('Bio::Map::MappableI');
			
			foreach my $pos (values %{$self->{'_positions'}}) {
                my $existing_element = $pos->element;
				if ($existing_element && $existing_element eq $element) {
					push(@positions, $pos);
				}
			}
		}
		else {
			@positions = values %{$self->{'_positions'}};
		}
		
        return sort { $a->sortable <=> $b->sortable } @positions;
    }
    return;
}

=head2 species

 Title   : species
 Usage   : my $species = $map->species;
 Function: Get/Set Species for a map
 Returns : L<Bio::Species> object
 Args    : (optional) Bio::Species

=cut

sub species{
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 units

 Title   : units
 Usage   : $map->units('cM');
 Function: Get/Set units for a map
 Returns : units for a map
 Args    : units for a map (string)

=cut

sub units{
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 type

 Title   : type
 Usage   : my $type = $map->type
 Function: Get/Set Map type
 Returns : String coding map type
 Args    : (optional) string

=cut

sub type {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 name

 Title   : name
 Usage   : my $name = $map->name
 Function: Get/Set Map name
 Returns : Map name
 Args    : (optional) string

=cut

sub name {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 length

 Title   : length
 Usage   : my $length = $map->length();
 Function: Retrieves the length of the map. 
           It is possible for the length to be unknown for maps such as
           Restriction Enzyme, will return undef in that case
 Returns : integer representing length of map in current units
           will return undef if length is not calculateable
 Args    : none

=cut

sub length{
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 unique_id

 Title   : unique_id
 Usage   : my $id = $map->unique_id;
 Function: Get/Set the unique ID for this map
 Returns : a unique identifier
 Args    : [optional] new identifier to set 

=cut

sub unique_id{
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 each_element

 Title   : each_element
 Usage   : my @elements = $map->each_element;
 Function: Retrieves all the elements on a map (unordered)
 Returns : Array of Map elements (L<Bio::Map::MappableI>)
 Args    : none

=cut

sub each_element{
    my $self = shift;
    $self->throw_not_implemented();
}

1;
