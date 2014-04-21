#
# BioPerl module for Bio::Map::MarkerI
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

Bio::Map::MarkerI - Interface for basic marker functionality

=head1 SYNOPSIS

    # do not use this module directly
    # See Bio::Map::Marker for an example of
    # implementation.

=head1 DESCRIPTION

A Marker is a Bio::Map::Mappable with some properties particular to markers.
It also offers a number of convienience methods to make dealing with map
elements easier.

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
Jason Stajich      jason@bioperl.org
Chad Matsalla      bioinformatics1@dieselwurks.com
Sendu Bala         bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Map::MarkerI;
use strict;

use base qw(Bio::Map::MappableI);

=head2 get_position_object

 Title   : get_position_class
 Usage   : my $position = $marker->get_position_object();
 Function: To get an object of the default Position class
           for this Marker. Subclasses should redefine this method.
           The Position returned needs to be a L<Bio::Map::PositionI> with
		   -element set to self.
 Returns : L<Bio::Map::PositionI>
 Args    : none for an 'empty' PositionI object, optionally
           Bio::Map::MapI and value string to set the Position's -map and -value
           attributes.

=cut

sub get_position_object {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 position

 Title   : position
 Usage   : my $position = $mappable->position();
		   $mappable->position($position);
 Function: Get/Set the Position of this Marker (where it is on which map),
           purging all other positions before setting.
 Returns : L<Bio::Map::PositionI>
 Args    : Bio::Map::PositionI
            OR
           Bio::Map::MapI AND
           scalar
            OR
           scalar, but only if the marker has a default map

=cut

sub position {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 positions

 Title   : positions
 Usage   : $marker->positions([$pos1, $pos2, $pos3]);
 Function: Add multiple Bio::Map::PositionI to this marker
 Returns : n/a
 Args    : array ref of $map/value tuples or array ref of Bio::Map::PositionI

=cut

sub positions {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 default_map

 Title   : default_map
 Usage   : my $map = $marker->default_map();
 Function: Get/Set the default map for the marker.
 Returns : L<Bio::Map::MapI>
 Args    : [optional] new L<Bio::Map::MapI>

=cut

sub default_map {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 in_map

 Title   : in_map
 Usage   : if ( $marker->in_map($map) ) {}
 Function: Tests if this marker is found on a specific map
 Returns : boolean
 Args    : a map unique id OR Bio::Map::MapI

=cut

1;
