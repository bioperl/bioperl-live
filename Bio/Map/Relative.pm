#
# BioPerl module for Bio::Map::Relative
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Sendu Bala
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::Relative - Represents what a Position's coordiantes are relative to.

=head1 SYNOPSIS

    # Get a Bio::Map::PositionI somehow
    my $pos = Bio::Map::Position->new(-value => 100);

    # its co-ordinates are implicitly relative to the start of its map
    my $implicit_relative = $pos->relative;
    my $type = $implicit_relative->type; # $type eq 'map'
    my $value = $implicit_relative->$type(); # $value == 0

    # make its co-ordinates relative to another Position
    my $pos_we_are_relative_to = Bio::Map::Position->new(-value => 200);
    my $relative = Bio::Map::Relative->new(-position => $pos_we_are_relative_to);
    $pos->relative($relative);

    # Get the start co-ordinate of $pos relative to $pos_we_are_relative_to
    my $start = $pos->start; # $start == 100

    # Get the start co-ordinate of $pos relative to the start of the map
    my $abs_start = $relative->absolute_conversion($pos); # $abs_start == 300
    # - or -
    $pos->absolute(1);
    my $abs_start = $pos->start; # $abs_start == 300
    $pos->absolute(0);

    # Get the start co-ordinate of $pos relative to a third Position
    my $pos_frame_of_reference = Bio::Map::Position->new(-value => 10);
    my $relative2 = Bio::Map::Relative->new(-position => $pos_frame_of_reference);
    my $start = $pos->start($relative2); # $start == 290

=head1 DESCRIPTION

A Relative object is used to describe what the co-ordinates (numerical(),
start(), end()) of a Position are relative to. By default they are
implicitly assumed to be relative to the start of the map the Position is on.
But setting the relative() of a Position to one of these objects lets us
define otherwise.

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

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::Relative;
use strict;
use Scalar::Util qw(looks_like_number);

use base qw(Bio::Root::Root Bio::Map::RelativeI);

=head2 new

 Title   : new
 Usage   : my $relative = Bio::Map::Relative->new();
 Function: Build a new Bio::Map::Relative object.
 Returns : Bio::Map::Relative object
 Args    : -map => int           : coordinates are relative to this point on the
                                   Position's map [default is map => 0, ie.
                                   relative to the start of the map],
           -element => Mappable  : or relative to this element's (a
                                   Bio::Map::MappableI) position in the map
                                   (only works if the given element has only one
                                   position in the map the Position belongs to),
           -position => Position : or relative to this other Position (a
                                   Bio::Map::PositionI, fails if the other
                                   Position is on a different map to this map)

           -description => string: Free text description of what this relative
                                   describes

           (To say a Position is relative to something and upstream of it,
            the Position's start() co-ordinate should be set negative)

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    my ($map, $element, $position, $desc) = 
	$self->_rearrange([qw( MAP ELEMENT POSITION DESCRIPTION )], @args);
    
    if (defined($map) + defined($element) + defined($position) > 1) {
        $self->throw("-map, -element and -position are mutually exclusive");
    }
    
    defined($map) && $self->map($map);
    $element && $self->element($element);
    $position && $self->position($position);
    $desc && $self->description($desc);
    
    return $self;
}

=head2 absolute_conversion

 Title   : absolute_conversion
 Usage   : my $absolute_coord = $relative->absolute_conversion($pos);
 Function: Convert the start co-ordinate of the supplied position into a number
           relative to the start of its map.
 Returns : scalar number
 Args    : Bio::Map::PositionI object

=cut

sub absolute_conversion {
    my ($self, $pos) = @_;
    $self->throw("Must supply an object") unless ref($pos);
    $self->throw("This is [$pos], not a Bio::Map::PositionI") unless $pos->isa('Bio::Map::PositionI');
    
    # get the raw start position of our position
    my $prior_abs = $pos->absolute;
    $pos->absolute(0) if $prior_abs;
    my $raw = $pos->start;
    $pos->absolute($prior_abs) if $prior_abs;
    $self->throw("Can't convert co-ordinates when start isn't set") unless defined($raw); #*** needed? return undef?
    
    # what are we relative to?
    my $type = $self->type;
    my $value = $self->$type;
    $self->throw("Details not yet set for this Relative, cannot convert") unless $type && defined($value);
    
    # get the absolute start of the thing we're relative to
    my $map = $pos->map;
    if ($type eq 'element') {
        $self->throw("Relative to a Mappable, but the Position has no map") unless $map;
        my @positions = $value->get_positions($map);
        $value = shift(@positions);
        $self->throw("Relative to a Mappable, but this Mappable has no positions on the supplied Position's map") unless $value;
    }
    if (ref($value)) {
        # psuedo-recurse
        my $rel = $value->relative;
        $value = $rel->absolute_conversion($value);
    }
    
    if (defined($value)) {
        return $value + $raw;
    }
    return;
}

=head2 type

 Title   : type
 Usage   : my $type = $relative->type();
 Function: Get the type of thing we are relative to. The types correspond
           to a method name, so the value of what we are relative to can
           subsequently be found by $value = $relative->$type;

           Note that type is set by the last method that was set, or during
           new().

 Returns : the string 'map', 'element' or 'position', or undef
 Args    : none

=cut

sub type {
    my $self = shift;
    return $self->{_use} || return;
}

=head2 map

 Title   : map
 Usage   : my $int = $relative->map();
           $relative->map($int);
 Function: Get/set the distance from the start of the map that the Position's
           co-ordiantes are relative to.
 Returns : int
 Args    : none to get, OR
           int to set; a value of 0 means relative to the start of the map.

=cut

sub map {
    my ($self, $num) = @_;
    if (defined($num)) {
        $self->throw("This is [$num], not a number") unless looks_like_number($num);
        $self->{_use} = 'map';
        $self->{_map} = $num;
    }
    return defined($self->{_map}) ? $self->{_map} : return;
}

=head2 element

 Title   : element
 Usage   : my $element = $relative->element();
           $relative->element($element);
 Function: Get/set the map element (Mappable) the Position is relative to. If
           the Mappable has more than one Position on the Position's map, we
           will be relative to the Mappable's first Position on the map.
 Returns : Bio::Map::MappableI
 Args    : none to get, OR
           Bio::Map::MappableI to set

=cut

sub element {
    my ($self, $element) = @_;
    if ($element) {
        $self->throw("Must supply an object") unless ref($element);
        $self->throw("This is [$element], not a Bio::Map::MappableI") unless $element->isa('Bio::Map::MappableI');
        $self->{_use} = 'element';
        $self->{_element} = $element;
    }
    return $self->{_element} || return;
}

=head2 position

 Title   : position
 Usage   : my $position = $relative->position();
           $relative->position($position);
 Function: Get/set the Position your Position is relative to. Your Position
           will be made relative to the start of this supplied Position. It
           makes no difference what maps the Positions are on.
 Returns : Bio::Map::PositionI
 Args    : none to get, OR
           Bio::Map::PositionI to set

=cut

sub position {
    my ($self, $pos) = @_;
    if ($pos) {
        $self->throw("Must supply an object") unless ref($pos);
        $self->throw("This is [$pos], not a Bio::Map::PositionI") unless $pos->isa('Bio::Map::PositionI');
        $self->{_use} = 'position';
        $self->{_position} = $pos;
    }
    return $self->{_position} || return;
}

=head2 description

 Title   : description
 Usage   : my $description = $relative->description();
           $relative->description($description);
 Function: Get/set a textual description of what this relative describes.
 Returns : string
 Args    : none to get, OR
           string to set

=cut

sub description {
    my $self = shift;
    if (@_) { $self->{desc} = shift }
    return $self->{desc} || '';
}

1;
