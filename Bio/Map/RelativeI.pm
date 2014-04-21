#
# BioPerl module for Bio::Map::RelativeI
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

Bio::Map::RelativeI - Interface for describing what a Position's coordiantes are
                      relative to.

=head1 SYNOPSIS

    # do not use this module directly
    # See Bio::Map::Relative for an example of
    # implementation.

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

package Bio::Map::RelativeI;
use strict;

use base qw(Bio::Root::RootI);

=head2 absolute_conversion

 Title   : absolute_conversion
 Usage   : my $absolute_coord = $relative->absolute_conversion($pos);
 Function: Convert the start co-ordinate of the supplied position into a number
           relative to the start of its map.
 Returns : scalar number
 Args    : Bio::Map::PositionI object

=cut

sub absolute_conversion {
    my $self = shift;
    $self->throw_not_implemented();
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
    $self->throw_not_implemented();
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
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 element

 Title   : element
 Usage   : my $element = $relative->element();
           $relative->element($element);
 Function: Get/set the map element (Mappable) the Position is relative to. If
           the Mappable has more than one Position on the Position's map, we
           will be relative to the Mappable's first Position on the map.
 Returns : Bio::Map::MappableI
 Args    : none got get, OR
           Bio::Map::MappableI to set

=cut

sub element {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 position

 Title   : position
 Usage   : my $position = $relative->position();
           $relative->position($position);
 Function: Get/set the Position your Position is relative to. Your Position
           will be made relative to the start of this supplied Position. It
           makes no difference what maps the Positions are on.
 Returns : Bio::Map::PositionI
 Args    : none got get, OR
           Bio::Map::PositionI to set

=cut

sub position {
    my $self = shift;
    $self->throw_not_implemented();
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
    $self->throw_not_implemented();
}

1;
