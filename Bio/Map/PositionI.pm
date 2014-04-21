#
# BioPerl module for Bio::Map::PositionI
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

Bio::Map::PositionI - Abstracts the notion of a position having a value in the context of a marker and a Map

=head1 SYNOPSIS

    # do not use this module directly
    # See Bio::Map::Position for an example of
    # implementation.

=head1 DESCRIPTION

This object stores one of the postions that a mappable object
(e.g. Marker) may have in a map.

Positions can have non-numeric values or other methods to store the locations,
so they have a method numeric() which does the conversion. numeric()
returns the position in a form that can be compared between other positions of
the same type. It is not necessarily a value suitable for sorting positions (it
may be the distance from the previous position); for that purpose the result of
sortable() should be used.

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

Email jason-at-bioperl.org

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
use strict;
use Bio::Map::PositionHandler;
use Bio::Map::Mappable;
use Scalar::Util qw(looks_like_number);

use base qw(Bio::Map::EntityI Bio::RangeI);

=head2 EntityI methods

 These are fundamental to coordination of Positions and other entities, so are
 implemented at the interface level

=cut

=head2 get_position_handler

 Title   : get_position_handler
 Usage   : my $position_handler = $entity->get_position_handler();
 Function: Gets a PositionHandlerI that $entity is registered with.
 Returns : Bio::Map::PositionHandlerI object
 Args    : none

=cut

sub get_position_handler {
    my $self = shift;
    unless (defined $self->{_eh}) {
        my $ph = Bio::Map::PositionHandler->new(-self => $self);
        $self->{_eh} = $ph;
        $ph->register;
    }
    return $self->{_eh};
}

=head2 PositionHandlerI-related methods

 These are fundamental to coordination of Positions and other entities, so are
 implemented at the interface level

=cut

=head2 map

 Title   : map
 Usage   : my $map = $position->map();
           $position->map($map);
 Function: Get/Set the map the position is in.
 Returns : L<Bio::Map::MapI>
 Args    : none to get
           new L<Bio::Map::MapI> to set

=cut

sub map {
    my ($self, $map) = @_;
    return $self->get_position_handler->map($map);
}

=head2 element

 Title   : element
 Usage   : my $element = $position->element();
           $position->element($element);
 Function: Get/Set the element the position is for.
 Returns : L<Bio::Map::MappableI>
 Args    : none to get
           new L<Bio::Map::MappableI> to set

=cut

sub element {
    my ($self, $element) = @_;
    return $self->get_position_handler->element($element);
}

=head2 marker

 Title   : marker
 Function: This is a synonym of the element() method
 Status  : deprecated, will be removed in the next version

=cut

*marker = \&element;

=head2 PositionI-specific methods

=cut

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
 Usage   : my $num = $position->numeric;
 Function: Read-only method that is guaranteed to return a numeric 
           representation of the start of this position. 
 Returns : scalar numeric
 Args    : none to get the co-ordinate normally (see absolute() method), OR
           Bio::Map::RelativeI to get the co-ordinate converted to be
           relative to what this Relative describes.

=cut

sub numeric {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 sortable

 Title   : sortable
 Usage   : my $num = $position->sortable();
 Function: Read-only method that is guaranteed to return a value suitable
           for correctly sorting this kind of position amongst other positions
           of the same kind on the same map. Note that sorting different kinds
           of position together is unlikely to give sane results.
 Returns : numeric
 Args    : none

=cut

sub sortable {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 relative

  Title   : relative
  Usage   : my $relative = $position->relative();
            $position->relative($relative);
  Function: Get/set the thing this Position's coordinates (numerical(), start(),
            end()) are relative to, as described by a Relative object.
  Returns : Bio::Map::RelativeI (default is one describing "relative to the
            start of the Position's map")
  Args    : none to get, OR
            Bio::Map::RelativeI to set

=cut

sub relative {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 absolute

  Title   : absolute
  Usage   : my $absolute = $position->absolute();
            $position->absolute($absolute);
  Function: Get/set how this Position's co-ordinates (numerical(), start(),
            end()) are reported. When absolute is off, co-ordinates are
            relative to the thing described by relative(). Ie. the value
            returned by start() will be the same as the value you set start()
            to. When absolute is on, co-ordinates are converted to be relative
            to the start of the map.

            So if relative() currently points to a Relative object describing
            "relative to another position which is 100 bp from the start of
            the map", this Position's start() had been set to 50 and absolute()
            returns 1, $position->start() will return 150. If absolute() returns
            0 in the same situation, $position->start() would return 50.

  Returns : boolean (default 0)
  Args    : none to get, OR
            boolean to set

=cut

sub absolute {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 RangeI-based methods

=cut

=head2 start

  Title   : start
  Usage   : my $start = $position->start();
            $position->start($start);
  Function: Get/set the start co-ordinate of this position.
  Returns : the start of this position
  Args    : scalar numeric to set, OR
            none to get the co-ordinate normally (see absolute() method), OR
            Bio::Map::RelativeI to get the co-ordinate converted to be
            relative to what this Relative describes.

=cut

=head2 end

  Title   : end
  Usage   : my $end = $position->end();
            $position->end($end);
  Function: Get/set the end co-ordinate of this position.
  Returns : the end of this position
  Args    : scalar numeric to set, OR
            none to get the co-ordinate normally (see absolute() method), OR
            Bio::Map::RelativeI to get the co-ordinate converted to be
            relative to what this Relative describes.

=cut

=head2 length

  Title   : length
  Usage   : $length = $position->length();
  Function: Get the length of this position.
  Returns : the length of this position
  Args    : none

=cut

=head2 strand

  Title   : strand
  Usage   : $strand = $position->strand();
  Function: Get the strand of this position; it is always 1 since maps to not
            have strands.
  Returns : 1
  Args    : none

=cut

sub strand {
    return 1;
}

=head2 toString

  Title   : toString
  Usage   : print $position->toString(), "\n";
  Function: stringifies this range
  Returns : a string representation of the range of this Position
  Args    : optional Bio::Map::RelativeI to have the co-ordinates reported
            relative to the thing described by that Relative

=cut

sub toString {
    my $self = shift;
    $self->throw_not_implemented();
}

=head1 RangeI-related methods

These methods work by considering only the values of start() and end(), as
modified by considering every such co-ordinate relative to the start of the map
(ie. absolute(1) is set temporarily during the calculation), or any supplied
Relative. For the boolean methods, when the comparison Position is on the same
map as the calling Position, there is no point supplying a Relative since the
answer will be the same as without. Relative is most useful when comparing
Positions on different maps and you have a Relative that describes some special
place on each map like 'the start of the gene', where the actual start of the
gene relative to the start of the map is different for each map.

The methods do not consider maps during their calculations - things on different
maps can overlap/contain/intersect/etc. each other.

The geometrical methods (intersect, union etc.) do things to the geometry of
ranges, and return Bio::Map::PositionI compliant objects or triplets (start,
stop, strand) from which new positions could be built. When a PositionI is made
it will have a map transferred to it if all the arguments shared the same map.
If a Relative was supplied the result will have that same Relative.

Note that the strand-testing args are there for compatibility with the RangeI
interface. They have no meaning when only using PositionI objects since maps do
not have strands. Typically you will just set the argument to undef if you want
to supply the argument after it.

=head2 equals

  Title   : equals
  Usage   : if ($p1->equals($p2)) {...}
  Function: Test whether $p1 has the same start, end, length as $p2.
  Returns : true if they are describing the same position (regardless of map)
  Args    : arg #1 = a Bio::RangeI (eg. a Bio::Map::Position) to compare this
                     one to (mandatory)
            arg #2 = optional strand-testing arg ('strong', 'weak', 'ignore')
            arg #3 = optional Bio::Map::RelativeI to ask if the Positions
                     equal in terms of their relative position to the thing
                     described by that Relative

=cut

sub equals {
    # overriding the RangeI implementation so we can handle Relative
    my ($self, $other, $so, $rel) = @_;
    
    my ($own_start, $own_end) = $self->_pre_rangei($self, $rel);
    my ($other_start, $other_end) = $self->_pre_rangei($other, $rel);
    
    return ($self->_testStrand($other, $so) and
            $own_start == $other_start and $own_end == $other_end);
}


=head2 less_than

 Title   : less_than
 Usage   : if ($position->less_than($other_position)) {...}
 Function: Ask if this Position ends before another starts.
 Returns : boolean
 Args    : arg #1 = a Bio::RangeI (eg. a Bio::Map::Position) to compare this
                    one to (mandatory)
           arg #2 = optional Bio::Map::RelativeI to ask if the Position is less
                    in terms of their relative position to the thing described
                    by that Relative

=cut

sub less_than {
    my ($self, $other, $rel) = @_;
    
    my ($own_start, $own_end) = $self->_pre_rangei($self, $rel);
    my ($other_start, $other_end) = $self->_pre_rangei($other, $rel);
    
    return $own_end < $other_start;
}

=head2 greater_than

 Title   : greater_than
 Usage   : if ($position->greater_than($other_position)) {...}
 Function: Ask if this Position starts after another ends.
 Returns : boolean
 Args    : arg #1 = a Bio::RangeI (eg. a Bio::Map::Position) to compare this
                    one to (mandatory)
           arg #2 = optional Bio::Map::RelativeI to ask if the Position is
                    greater in terms of their relative position to the thing
                    described by that Relative

=cut

sub greater_than {
    my ($self, $other, $rel) = @_;
    
    my ($own_start, $own_end) = $self->_pre_rangei($self, $rel);
    my ($other_start, $other_end) = $self->_pre_rangei($other, $rel);
    
    return $own_start > $other_end;
}

=head2 overlaps

  Title   : overlaps
  Usage   : if ($p1->overlaps($p2)) {...}
  Function: Tests if $p1 overlaps $p2.
  Returns : True if the positions overlap (regardless of map), false otherwise
  Args    : arg #1 = a Bio::RangeI (eg. a Bio::Map::Position) to compare this
                     one to (mandatory)
            arg #2 = optional strand-testing arg ('strong', 'weak', 'ignore')
            arg #3 = optional Bio::Map::RelativeI to ask if the Positions
                     overlap in terms of their relative position to the thing
                     described by that Relative
            arg #4 = optional minimum percentage length of the overlap before
                     reporting an overlap exists (default 0)

=cut

sub overlaps {
    # overriding the RangeI implementation so we can handle Relative
    my ($self, $other, $so, $rel, $min_percent) = @_;
    $min_percent ||= 0;
    
    my ($own_min, $other_min) = (0, 0);
    if ($min_percent > 0) {
        $own_min = (($self->length / 100) * $min_percent) - 1;
        $other_min = (($other->length / 100) * $min_percent) - 1;
    }
    
    my ($own_start, $own_end) = $self->_pre_rangei($self, $rel);
    my ($other_start, $other_end) = $self->_pre_rangei($other, $rel);
    
    return ($self->_testStrand($other, $so) and not
            (($own_start + $own_min > $other_end or $own_end - $own_min < $other_start) ||
             ($own_start > $other_end - $other_min or $own_end < $other_start + $other_min)));
}

=head2 contains

  Title   : contains
  Usage   : if ($p1->contains($p2)) {...}
  Function: Tests whether $p1 totally contains $p2.
  Returns : true if the argument is totally contained within this position
            (regardless of map), false otherwise
  Args    : arg #1 = a Bio::RangeI (eg. a Bio::Map::Position) to compare this
                     one to, or scalar number (mandatory)
            arg #2 = optional strand-testing arg ('strong', 'weak', 'ignore')
            arg #3 = optional Bio::Map::RelativeI to ask if the Position
                     is contained in terms of their relative position to the
                     thing described by that Relative

=cut

sub contains {
    # overriding the RangeI implementation so we can handle Relative
    my ($self, $other, $so, $rel) = @_;
    
    my ($own_start, $own_end) = $self->_pre_rangei($self, $rel);
    my ($other_start, $other_end) = $self->_pre_rangei($other, $rel);
    
    return ($self->_testStrand($other, $so) and
			$other_start >= $own_start and $other_end <= $own_end);
}

=head2 intersection

 Title   : intersection
 Usage   : ($start, $stop, $strand) = $p1->intersection($p2)
           ($start, $stop, $strand) = Bio::Map::Position->intersection(\@positions);
           $mappable = $p1->intersection($p2, undef, $relative);
           $mappable = Bio::Map::Position->intersection(\@positions);
 Function: gives the range that is contained by all ranges
 Returns : undef if they do not overlap, OR
           Bio::Map::Mappable object who's positions are the
           cross-map-calculated intersection of the input positions on all the
           maps that the input positions belong to, OR, in list context, a three
           element array (start, end, strand)
 Args    : arg #1 = [REQUIRED] a Bio::RangeI (eg. a Bio::Map::Position) to
                    compare this one to, or an array ref of Bio::RangeI
           arg #2 = optional strand-testing arg ('strong', 'weak', 'ignore')
           arg #3 = optional Bio::Map::RelativeI to ask how the Positions
                    intersect in terms of their relative position to the thing
                    described by that Relative

=cut

sub intersection {
    # overriding the RangeI implementation so we can transfer map and handle
    # Relative
    my ($self, $given, $so, $rel) = @_;
	$self->throw("missing arg: you need to pass in another argument") unless $given;
    
    my @positions;
    if ($self eq "Bio::Map::PositionI") {
		$self = "Bio::Map::Position";
		$self->warn("calling static methods of an interface is deprecated; use $self instead");
	}
	if (ref $self) {
		push(@positions, $self);
	}
    ref($given) eq 'ARRAY' ? push(@positions, @{$given}) : push(@positions, $given);
    $self->throw("Need at least 2 Positions") unless @positions >= 2;
    
    my ($intersect, $i_start, $i_end, $c_start, $c_end, %known_maps);
    while (@positions > 0) {
        unless ($intersect) {
            $intersect = shift(@positions);
            ($i_start, $i_end) = $self->_pre_rangei($intersect, $rel);
            my $map = $intersect->map;
            $known_maps{$map->unique_id} = $map;
        }
        
        my $compare = shift(@positions);
        ($c_start, $c_end) = $self->_pre_rangei($compare, $rel);
        return unless $compare->_testStrand($intersect, $so);
        if ($compare->isa('Bio::Map::PositionI')) {
            my $this_map = $compare->map;
            if ($this_map) {
                $known_maps{$this_map->unique_id} = $this_map;
            }
        }
        else {
            $self->throw("Only Bio::Map::PositionI objects are supported, not [$compare]");
        }
        
        my @starts = sort {$a <=> $b} ($i_start, $c_start);
        my @ends   = sort {$a <=> $b} ($i_end, $c_end);
        
        my $start = pop @starts; # larger of the 2 starts
        my $end = shift @ends;   # smaller of the 2 ends
        
        my $intersect_strand;    # strand for the intersection
        if (defined($intersect->strand) && defined($compare->strand) && $intersect->strand == $compare->strand) {
            $intersect_strand = $compare->strand;
        }
        else {
            $intersect_strand = 0;
        }
        
        if ($start > $end) {
            return;
        }
        else {
            $intersect = $self->new(-start  => $start,
                                    -end    => $end,
                                    -strand => $intersect_strand);
        }
    }
    
    $intersect || return;
    my ($start, $end, $strand) = ($intersect->start, $intersect->end, $intersect->strand);
    
    my @intersects;
    foreach my $known_map (values %known_maps) {
        my $new_intersect = $intersect->new(-start => $start,
                                            -end => $end,
                                            -strand => $strand,
                                            -map => $known_map);
        $new_intersect->relative($rel) if $rel;
        push(@intersects, $new_intersect);
    }
    unless (@intersects) {
        $intersect->relative($rel) if $rel;
        @intersects = ($intersect);
    }
    
    my $result = Bio::Map::Mappable->new();
    $result->add_position(@intersects); # sneaky, add_position can take a list of positions
    return $result;
}

=head2 union

 Title   : union
 Usage   : ($start, $stop, $strand) = $p1->union($p2);
           ($start, $stop, $strand) = Bio::Map::Position->union(@positions);
           my $mappable = $p1->union($p2);
           my $mappable = Bio::Map::Position->union(@positions);
 Function: finds the minimal position/range that contains all of the positions
 Returns : Bio::Map::Mappable object who's positions are the
           cross-map-calculated union of the input positions on all the maps
           that the input positions belong to, OR, in list context, a three
           element array (start, end, strand)
 Args    : a Bio::Map::PositionI to compare this one to, or a list of such
           OR
           a single Bio::Map::PositionI or array ref of such AND a
           Bio::Map::RelativeI to ask for the Position's union in terms of their
           relative position to the thing described by that Relative

=cut

sub union {
    # overriding the RangeI implementation so we can transfer map and handle
    # Relative
    my ($self, @args) = @_;
    $self->throw("Not enough arguments") unless @args >= 1;
    
    my @positions;
    my $rel;
    if ($self eq "Bio::Map::PositionI") {
		$self = "Bio::Map::Position";
		$self->warn("calling static methods of an interface is deprecated; use $self instead");
	}
	if (ref $self) {
		push(@positions, $self);
	}
    if (ref $args[0] eq 'ARRAY') {
        push(@positions, @{shift(@args)});
    }
    else {
        push(@positions, shift(@args));
    }
    if ($args[0] && $args[0]->isa('Bio::Map::RelativeI')) {
        $rel = shift(@args);
    }
    foreach my $arg (@args) {
        # avoid pushing undefined values into @positions
        push(@positions, $arg) if $arg;
    }
    $self->throw("Need at least 2 Positions") unless @positions >= 2;
    
    my (@starts, @ends, %known_maps, $union_strand);
    foreach my $compare (@positions) {
        # RangeI union allows start or end to be undefined; however _pre_rangei
        # will throw
        my ($start, $end) = $self->_pre_rangei($compare, $rel);
        
        if ($compare->isa('Bio::Map::PositionI')) {
            my $this_map = $compare->map;
            if ($this_map) {
                $known_maps{$this_map->unique_id} = $this_map;
            }
        }
        else {
            $self->throw("Only Bio::Map::PositionI objects are supported, not [$compare]");
        }
        
        if (! defined $union_strand) {
			$union_strand = $compare->strand;
		}
        else {
			if (! defined $compare->strand or $union_strand ne $compare->strand) {
				$union_strand = 0;
			}
		}
        
        push(@starts, $start);
        push(@ends, $end);
    }
    
	@starts = sort { $a <=> $b } @starts;
	@ends   = sort { $a <=> $b } @ends;
	my $start = shift @starts;
	my $end = pop @ends;
    
    my @unions;
    foreach my $known_map (values %known_maps) {
        my $new_union = $self->new(-start => $start,
                                   -end => $end,
                                   -strand => $union_strand,
                                   -map => $known_map);
        $new_union->relative($rel) if $rel;
        push(@unions, $new_union);
    }
    unless (@unions) {
        @unions = ($self->new(-start => $start,
                         -end => $end,
                         -strand => $union_strand));
        $unions[0]->relative($rel) if $rel;
    }
    
    my $result = Bio::Map::Mappable->new();
    $result->add_position(@unions); # sneaky, add_position can take a list of positions
    return $result;
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
 Usage   : my @disc_ranges = Bio::Map::Position->disconnected_ranges(@ranges);
 Function: Creates the minimal set of positions such that each input position is
           fully contained by at least one output position, and none of the
           output positions overlap.
 Returns : Bio::Map::Mappable with the calculated disconnected ranges
 Args    : a Bio::Map::PositionI to compare this one to, or a list of such,
           OR
           a single Bio::Map::PositionI or array ref of such AND a
           Bio::Map::RelativeI to consider all Position's co-ordinates in terms
           of their relative position to the thing described by that Relative,
           AND, optionally, an int for the minimum percentage of overlap that
           must be present before considering two ranges to be overlapping
           (default 0)

=cut

sub disconnected_ranges {
    # overriding the RangeI implementation so we can transfer map and handle
    # Relative
    my ($self, @args) = @_;
    $self->throw("Not enough arguments") unless @args >= 1;
    
    my @positions;
    my $rel;
    my $overlap = 0;
    if ($self eq "Bio::Map::PositionI") {
		$self = "Bio::Map::Position";
		$self->warn("calling static methods of an interface is deprecated; use $self instead");
	}
	if (ref $self) {
		push(@positions, $self);
	}
    if (ref $args[0] eq 'ARRAY') {
        push(@positions, @{shift(@args)});
    }
    else {
        push(@positions, shift(@args));
    }
    if ($args[0] && $args[0]->isa('Bio::Map::RelativeI')) {
        $rel = shift(@args);
        $overlap = shift(@args);
    }
    foreach my $arg (@args) {
        push(@positions, $arg) if $arg;
    }
    $self->throw("Need at least 2 Positions") unless @positions >= 2;
    
    my %known_maps;
    foreach my $pos (@positions) {
        $pos->isa('Bio::Map::PositionI') || $self->throw("Must supply only Bio::Map::PositionI objects, not [$pos]");
        my $map = $pos->map || next;
        $known_maps{$map->unique_id} = $map;
    }
    my %prior_positions;
    foreach my $map (values %known_maps) {
        foreach my $pos ($map->get_positions) {
            $prior_positions{$pos} = 1;
        }
    }
    
    my @outranges = ();
    foreach my $inrange (@positions) {
        my @outranges_new = ();
        my %overlapping_ranges = ();
        
        for (my $i=0; $i<@outranges; $i++) {
            my $outrange = $outranges[$i];
            if ($inrange->overlaps($outrange, undef, $rel, $overlap)) {
                my $union_able = $inrange->union($outrange, $rel); # using $inrange->union($outrange, $rel); gives >6x speedup,
                                                                   # but different answer, not necessarily incorrect...
                foreach my $pos ($union_able->get_positions) {
                    $overlapping_ranges{$pos->toString} = $pos; # we flatten down to a result on a single map
                                                                # to avoid creating 10s of thousands of positions during this process;
                                                                # we then apply the final answer to all maps at the very end
                    last;
                }
            }
            else {
                push(@outranges_new, $outrange);
            }
        }
        
        @outranges = @outranges_new;
        
        my @overlappers = values %overlapping_ranges;
        if (@overlappers) {
            if (@overlappers > 1) {
                my $merged_range_able = shift(@overlappers)->union(\@overlappers, $rel);
                push(@outranges, $merged_range_able->get_positions);
            }
            else {
                push(@outranges, @overlappers);
            }
        }
        else {
            push(@outranges, $self->new(-start => $inrange->start($rel), -end => $inrange->end($rel), -strand => $inrange->strand, -map => $inrange->map, -relative => $rel));
        }
    }
    
    # purge positions that were created whilst calculating the answer, but
    # aren't the final answer and weren't there previously
    my %answers = map { $_ => 1 } @outranges;
    foreach my $map (values %known_maps) {
        foreach my $pos ($map->get_positions) {
            if (! exists $prior_positions{$pos} && ! exists $answers{$pos}) {
                $map->purge_positions($pos);
            }
        }
    }
    
    my %post_positions;
    foreach my $map (values %known_maps) {
        foreach my $pos ($map->get_positions) {
            $post_positions{$pos} = 1;
        }
    }
    
    @outranges || return;
    
    # make an outrange on all known maps
    my @final_positions;
    foreach my $map (values %known_maps) {
        foreach my $pos (@outranges) {
            if ($pos->map eq $map) {
                push(@final_positions, $pos);
            }
            else {
                push(@final_positions, $pos->new(-start => $pos->start,
                                                 -end => $pos->end,
                                                 -relative => $pos->relative,
                                                 -map => $map));
            }
        }
    }
    
    # assign the positions to a result mappable
    my $result = Bio::Map::Mappable->new();
    $result->add_position(@final_positions); # sneaky, add_position can take a list of positions
    return $result;
}

# get start & end suitable for rangeI methods, taking relative into account
sub _pre_rangei {
    my ($self, $other, $rel) = @_;
    $self->throw("Must supply an object") unless $other;
    if ($rel) {
        $self->throw("Must supply an object for the Relative argument") unless ref($rel);
        $self->throw("This is [$rel], not a Bio::Map::RelativeI") unless $rel->isa('Bio::Map::RelativeI');
    }
    
    my ($other_start, $other_end);
    if (ref($other)) {
        if (ref($other) eq 'ARRAY') {
            $self->throw("_pre_rangei got an array");
        }
        $self->throw("This is [$other], not a Bio::RangeI object") unless defined $other && $other->isa('Bio::RangeI');
        
        if ($other->isa('Bio::Map::PositionI')) {
            # to get the desired start/end we need the position to be on a map;
            # if it isn't on one temporarily place it on self's map
            # - this lets us have 'generic' positions that aren't on any map
            # but have a relative defined and can thus be usefully compared to
            # positions that /are/ on maps
            my $other_map = $other->map;
            unless ($other_map) {
                my $self_map = $self->map || $self->throw("Trying to compare two positions but neither had been placed on a map");
                $other->map($self_map);
            }
            
            # want start and end positions relative to the supplied rel or map start
            $rel ||= $other->absolute_relative;
            $other_start = $other->start($rel);
            $other_end = $other->end($rel);
            
            unless ($other_map) {
                $self->map->purge_positions($other);
            }
        }
        else {
            $other_start = $other->start;
            $other_end = $other->end;
        }
    }
    else {
        $self->throw("not a number") unless looks_like_number($other);
        $other_start = $other_end = $other;
    }
    
	$other->throw("start is undefined") unless defined $other_start;
	$other->throw("end is undefined") unless defined $other_end;
    
    return ($other_start, $other_end);
}

1;
