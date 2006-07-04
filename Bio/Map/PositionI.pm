# $Id$
#
# BioPerl module for Bio::Map::PositionI
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::PositionI - Abstracts the notion of a position having a value in the context of a marker and a Map

=head1 SYNOPSIS

    # do not use directly
    my $marker_position = $position->value;
    my $map             = $position->map;
    my $marker          = $position->marker;

=head1 DESCRIPTION

This object stores one of the postions that a mappable object
(e.g. Marker) may have in a map.

The method numeric() returns the position in a form that can be
compared between other positions of the same type.

A 'position', in addition to being a single point, can also be an area and so
can be imagined as a range and compared with other positions on the basis of
overlap, intersection etc.

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

Email jason-at-bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Lincoln Stein, lstein-at-cshl.org
Heikki Lehvaslaiho, heikki-at-bioperl-dot-org
Sendu Bala, bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Map::PositionI;
use vars qw(@ISA);
use strict;
use Bio::Root::RootI;
use Bio::RangeI;
use Carp;

@ISA = qw(Bio::Root::RootI Bio::RangeI);


=head2 map

 Title   : map
 Usage   : my $map = $position->map();
           $position->map($map);
 Function: Get/Set the map the position is in. When setting, notifies the map
           in question we belong to it, and any previous map that we no longer
           belong to it.
 Returns : L<Bio::Map::MapI>
 Args    : none to get
           new L<Bio::Map::MapI> to set

=cut

sub map {
    # this is fundamental to coordination of Positions and Maps, so is
    # implemented at the interface level
    my ($self, $map) = @_;
    
    if (defined $map) {
        $self->throw("This is [$map], not an object") unless ref($map);
        $self->throw("This is [$map], not a Bio::Map::MapI object") unless $map->isa('Bio::Map::MapI');
        
        if (defined $self->{'_map'} && $self->{'_map'} ne $map) {
            $self->{'_map'}->purge_positions($self);
        }
        $self->{'_map'} = $map;
        $map->add_position($self);
    }
    
    return $self->{'_map'} || return;
}

=head2 purge_map

 Title   : purge_map
 Usage   : $position->purge_map();
 Function: Disassociate this Position from any map. Notifies any pre-existing
           map that we are no longer on it.
 Returns : n/a
 Args    : none

=cut

sub purge_map {
    # this is fundamental to coordination of Positions and Maps, so is
    # implemented at the interface level
    my $self = shift;
    
    if (defined $self->{'_map'}) {
        $self->{'_map'}->purge_positions($self);
    }
    delete $self->{'_map'};
}

=head2 element

 Title   : element
 Usage   : my $element = $position->element();
           $position->element($element);
 Function: Get/Set the element the position is for. When setting, notifies the
           element in question we belong to it, and any previous element that we
           no longer belong to it.
 Returns : L<Bio::Map::MappableI>
 Args    : none to get
           new L<Bio::Map::MappableI> to set

=cut

sub element {
    # this is fundamental to coordination of Positions and Mappables, so is
    # implemented at the interface level
    my ($self, $element) = @_;
    
    if (defined $element) {
        $self->throw("This is [$element], not an object") unless ref($element);
        $self->throw("This is [$element], not a Bio::Map::MappableI object") unless $element->isa('Bio::Map::MappableI');
        
        if (defined $self->{'_element'} && $self->{'_element'} ne $element) {
            $self->{'_element'}->purge_positions($self);
        }
        $self->{'_element'} = $element;
        $element->add_position($self);
    }
    
    return $self->{'_element'} || return;
}

=head2 purge_element

 Title   : purge_element
 Usage   : $position->purge_element();
 Function: Disassociate this Position from any element. Notifies any
           pre-existing element that we no longer belong to it.
 Returns : n/a
 Args    : none

=cut

sub purge_element {
    # this is fundamental to coordination of Positions and Mappables, so is
    # implemented at the interface level
    my $self = shift;
    
    if (defined $self->{'_element'}) {
        $self->{'_element'}->purge_positions($self);
    }
    delete $self->{'_element'};
}

=head2 marker

 Title   : marker
 Function: This is a synonym of the element() method
 Status  : deprecated, will be removed in the next version

=cut

*marker = \&element;
*marker = \&element; # avoid warning

=head2 value

 Title   : value
 Usage   : my $pos = $position->value();
 Function: Get/Set the value for this position
 Returns : scalar, value
 Args    : [optional] new value to set

=cut

sub value {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 numeric

 Title   : numeric
 Usage   : my $num = $position->numeric();
 Function: Read-only method that is guaranteed to return 
           representation for this position that can be compared with others
 Returns : numeric
 Args    : none

=cut

sub numeric {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 sortable

 Title   : sortable
 Usage   : my $num = $position->sortable();
 Function: Read-only method that is guaranteed to return a value suitable
           for correctly sorting this kind of position
 Returns : numeric
 Args    : none

=cut

sub sortable {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 start

  Title   : start
  Usage   : $start = $position->start();
  Function: get/set the start of this position
  Returns : the start of this position
  Args    : optionally allows the start to be set
            using $position->start($start)

=cut

=head2 end

  Title   : end
  Usage   : $end = $position->end();
  Function: get/set the end of this position
  Returns : the end of this position
  Args    : optionally allows the end to be set
            using $position->end($end)

=cut

=head2 length

  Title   : length
  Usage   : $length = $position->length();
  Function: get the length of this position
  Returns : the length of this position
  Args    : none

=cut

=head2 strand

  Title   : strand
  Usage   : $strand = $position->strand();
  Function: get the strand of this position; it is always 1 (the forward strand)
  Returns : 1, the strandedness
  Args    : none

=cut

sub strand {
    return 1;
}

=head1 Boolean Methods

These range-related methods return true or false. They throw an error if start
and end are not defined.


=head2 overlaps

  Title   : overlaps
  Usage   : if ($p1->overlaps($p2)) {...}
  Function: tests if $p1 overlaps $p2
  Args    : a Bio::RangeI (eg. a Bio::Map::Position) to compare this one to
  Returns : true if the positions overlap (regardless of map), false otherwise
            [*** relative handling ***]

=cut

=head2 contains

  Title   : contains
  Usage   : if ($p1->contains($p2) {...}
  Function: tests whether $p1 totally contains $p2 
  Args    : a Bio::RangeI (eg. a Bio::Map::Position) to compare this one to.
	        alternatively, integer scalar to test
  Returns : true if the argument is totally contained within this position
            (regardless of map), false otherwise
            [*** relative handling ***]

=cut

=head1 Geometrical methods

These methods do things to the geometry of ranges, and return
Bio::Map::PositionI compliant objects or triplets (start, stop, strand) from
which new positions could be built.

=head2 intersection

 Title   : intersection
 Usage   : ($start, $stop, $strand) = $r1->intersection($r2)
           ($start, $stop, $strand) = Bio::Map::Position->intersection(\@positions);
           my $containing_range = $r1->intersection($r2);
           my $containing_range = Bio::Map::Position->intersection(\@positions);
 Function: gives the range that is contained by all ranges
 Returns : undef if they do not overlap, or
           the range that they do overlap in an object like the calling one
           (with map transfered across from self or a position in the array ref
           of positions supplied, as long as the map is shared by all inputs),
           OR a three element array
 Args    : a position/range to compare this one to, or an array ref of
           positions/ranges
           [*** relative handling ***]

=cut

sub intersection {
    # overriding the RangeI implementation so we can transfer map and make the
    # result a PositionI
    my $self = shift;
    
    if (wantarray()) {
        return $self->SUPER::intersection(@_);
    }
    
    my $range = $self->SUPER::intersection(@_);
    $range or return;
    
    # transfer the map
    my @things;
    if ($self eq "Bio::Map::PositionI") {
        $self = "Bio::Map::Position";
        $self->warn("calling static methods of an interface is deprecated; use $self instead");
    }
    if (ref $self) {
        push(@things, $self);
    }
    my $given = shift;
    ref($given) eq 'ARRAY' ? push(@things, @{$given}) : push(@things, $given);
    
    my $map;
    foreach my $thing (@things) {
        if ($thing->isa("Bio::Map::PositionI")) {
            my $this_map = $thing->map || next;
            $map ||= $this_map;
            
            if ($map ne $this_map) {
                $map = undef;
                last;
            }
        }
    }
    
    return $self->new(-start => $range->start, -end => $range->end, -map => $map);
}

=head2 union

 Title   : union
 Usage   : ($start, $stop, $strand) = $r1->union($r2);
           ($start, $stop, $strand) = Bio::Map::Position->union(@positions);
           my $newpos = Bio::Map::Position->union(@positions);
 Function: finds the minimal position/range that contains all of the positions
 Returns : the position object containing all of the ranges in an object like
           the calling one (with map transfered across from self or a position
           in the array ref of positions supplied, as long as the map is shared
           by all inputs), OR a three element array
 Args    : a range or list of positions/ranges to find the union of
           [*** relative handling ***]

=cut
sub union {
    # overriding the RangeI implementation so we can transfer map and make the
    # result a PositionI
    my $self = shift;
    
    if (wantarray()) {
        return $self->SUPER::union(@_);
    }
    
    my $range = $self->SUPER::union(@_);
    $range or return;
    
    # transfer the map
    my @things;
    if ($self eq "Bio::Map::PositionI") {
        $self = "Bio::Map::Position";
        $self->warn("calling static methods of an interface is deprecated; use $self instead");
    }
    if (ref $self) {
        push(@things, $self);
    }
    my $given = shift;
    ref($given) eq 'ARRAY' ? push(@things, @{$given}) : push(@things, $given);
    
    my $map;
    foreach my $thing (@things) {
        if ($thing->isa("Bio::Map::PositionI")) {
            my $this_map = $thing->map || next;
            $map ||= $this_map;
            
            if ($map ne $this_map) {
                $map = undef;
                last;
            }
        }
    }
    
    return $self->new(-start => $range->start, -end => $range->end, -map => $map);
}

=head2 overlap_extent

 Title   : overlap_extent
 Usage   : ($a_unique,$common,$b_unique) = $a->overlap_extent($b)
 Function: Provides actual amount of overlap between two different
           positions
 Example :
 Returns : array of values containing the length unique to the calling 
           position, the length common to both, and the length unique to 
           the argument position
 Args    : a position

=cut

#*** should this be overridden from RangeI?

=head2 disconnected_ranges

    Title   : disconnected_ranges
    Usage   : my @disc_ranges = Bio::Range->disconnected_ranges(@ranges);
    Function: finds the minimal set of ranges such that each input range
              is fully contained by at least one output range, and none of
              the output ranges overlap
    Args    : a list of ranges
    Returns : a list of objects of the same type as the input 
              (conforms to RangeI)

=cut

#*** this should be overridden from RangeI

1;
