# $Id$
#
# BioPerl module for Bio::Map::Marker
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Sendu Bala
# 
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::Mappable - An object representing a generic map element
that can have multiple locations in several maps.

=head1 SYNOPSIS

  # a map element in two different positions on the same map
  $map1 = new Bio::Map::SimpleMap ();
  $position1 = new Bio::Map::Position (-map => $map1, $value => 100);
  $position2 = new Bio::Map::Position (-map => $map1, $value => 200);
  $mappable = new Bio::Map::Mappable (-positions => [$position1, $position2] );

  # add another position on a different map
  $map2 = new Bio::Map::SimpleMap ();
  $position3 = new Bio::Map::Position (-map => $map2, $value => 50);
  $mappable->add_position($position3);

  # get all the places our map element is found, on a particular map of interest
  foreach $pos ($mappable->each_position($map1)) {
     print $pos->value, "\n";
  }

=head1 DESCRIPTION

This object handles the notion of a generic map element. Mappables are
entities with one or more positions on one or more maps.

This object is a pure perl implementation of L<Bio::Map::MappableI> with a
convienience method positions() to add multiple positions at once.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 CONTRIBUTORS

Heikki Lehvaslaiho heikki-at-bioperl-dot-org
Lincoln Stein      lstein@cshl.org
Jason Stajich      jason@bioperl.org
Chad Matsalla      bioinformatics1@dieselwurks.com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Map::Mappable;
use vars qw(@ISA);
use strict;
use Bio::Root::Root;
use Bio::Map::MappableI;

@ISA = qw(Bio::Root::Root Bio::Map::MappableI);

=head2 new

 Title   : new
 Usage   : my $mappable = new Bio::Map::Mappable();
 Function: Builds a new Bio::Map::Mappable object
 Returns : Bio::Map::Mappable
 Args    : none

=cut

sub new {
    my $class = shift;
    my $self = $class->SUPER::new();
    return bless($self, ref $class || $class);
}

=head2 known_maps

 Title   : known_maps
 Usage   : my @maps = $marker->known_maps()
 Function: Returns the maps that this mappable is found on
 Returns : Array of L<Bio::Map::MapI> objects
 Args    : none

=cut

sub known_maps {
	my $self = shift;
	my %maps;
    foreach my $pos ($self->each_position) {
        my $map = $pos->map;
        if ($map) {
            $maps{$map} = $map;
        }
    }
    
    return values %maps;
}

=head2 in_map

 Title   : in_map
 Usage   : if ($marker->in_map($map)) {...}
 Function: Tests if this mappable is found on a specific map
 Returns : boolean
 Args    : L<Bio::Map::MapI>

=cut

sub in_map {
	my ($self, $query_map) = @_;
	$self->throw("Must supply an argument") unless $query_map;
    $self->throw("This is [$query_map], not an object") unless ref($query_map);
    $self->throw("This is [$query_map], not a Bio::Map::MapI object") unless $query_map->isa('Bio::Map::MapI');
    
    foreach my $map ($self->known_maps) {
        ($map eq $query_map) && return 1;
    }
    
    return 0;
}

=head2 Comparison methods

=cut

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
    
    (@{$mine} > 0 && @{$yours} > 0) or return;
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
    
    (@{$mine} > 0 && @{$yours} > 0) or return;
    
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

=head2 tuple

 Title   : tuple
 Usage   : Do Not Use!
 Function: tuple was supposed to be a private method; this method no longer
           does anything
 Returns : warning
 Args    : none
 Status  : deprecated, will be removed in next version

=cut

sub tuple {
    my $self = shift;
    $self->warn("The tuple method was supposed to be a private method, don't call it!");
}

1;
