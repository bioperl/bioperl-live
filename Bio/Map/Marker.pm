#
# BioPerl module for Bio::Map::Marker
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::Marker - An central map object representing a generic marker
that can have multiple location in several maps.

=head1 SYNOPSIS

  # get map objects somehow

  # a marker with complex localisation
  $o_usat = Bio::Map::Marker->new(-name=>'Chad Super Marker 2',
				  -positions => [ [$map1, $position1],
                                                  [$map1, $position2]
						] );

  # The markers deal with Bio::Map::Position objects which can also
  # be explicitly created and passed on to markers as an array ref:
  $o_usat2 = Bio::Map::Marker->new(-name=>'Chad Super Marker 3',
				  -positions => [ $pos1, 
                                                  $pos2
						] );

  # a marker with unique position in a map
  $marker1 = Bio::Map::Marker->new(-name=>'hypervariable1',
				   -map => $map1,
				   -position => 100
				   );

  # another way of creating a marker with unique position in a map:
  $marker2 = Bio::Map::Marker->new(-name=>'hypervariable2');
  $map1->add_element($marker2);
  $marker2->position(100);

  # position method is a short cut for get/setting unique positions
  # which overwrites previous values
  # to place a marker to other maps or to have multiple positions
  # for a map within the same map use add_position()

  $marker2->add_position(200);	# new position in the same map
  $marker2->add_position($map2,200); # new map

  # setting a map() in a marker or adding a marker into a map are
  # identical mathods. Both set the bidirectional connection which is
  # used by the marker to remember its latest, default map.

  # Regardes of how marker positions are created, they are stored and
  # returned as Bio::Map::PositionI objects:

  # unique position
  print $marker1->position->value, "\n";
  # several positions
  foreach $pos ($marker2->each_position($map1)) {
     print $pos->value, "\n";
  }

See L<Bio::Map::Position> and L<Bio::Map::PositionI> for more information.

=head1 DESCRIPTION

A Marker is a Bio::Map::Mappable with some properties particular to markers.
It also offers a number of convienience methods to make dealing with map
elements easier.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing list.  Your participation is much appreciated.

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

=head1 AUTHOR - Chad Matsalla

Email bioinformatics1@dieselwurks.com

=head1 CONTRIBUTORS

Heikki Lehvaslaiho heikki-at-bioperl-dot-org
Lincoln Stein      lstein@cshl.org
Jason Stajich      jason@bioperl.org
Sendu Bala         bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::Marker;
use strict;
use Bio::Map::Position;

use base qw(Bio::Map::Mappable Bio::Map::MarkerI);

=head2 new

 Title   : new
 Usage   : my $marker = Bio::Map::Marker->new( -name => 'Whizzy marker',
	                                          -position => $position);
 Function: Builds a new Bio::Map::Marker object
 Returns : Bio::Map::Marker
 Args    :
           -name    => name of this microsatellite
                       [optional], string,default 'Unknown'
           -default_map => the default map for this marker, a Bio::Map::MapI
           -position => map position for this marker, a Bio::Map::PositionI
           -positions => array ref of Bio::Map::PositionI objects

           position and positions can also take as values anything the
           corresponding methods can take

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
	bless($self, ref $class || $class);

    my ($name, $default_map, $map, $position, $positions) =
	$self->_rearrange([qw(NAME
				  DEFAULT_MAP
			      MAP
			      POSITION
			      POSITIONS
			      )], @args);

    if ($name) { $self->name($name); }
    else {$self->name('Unnamed marker'); }

    $map         && $self->default_map($map);
	$default_map && $self->default_map($default_map);
    $position    && $self->position($position);
    $positions   && $self->positions($positions);

    return $self;
}

=head2 default_map

 Title   : default_map
 Usage   : my $map = $marker->default_map();
 Function: Get/Set the default map for the marker.
 Returns : L<Bio::Map::MapI>
 Args    : [optional] new L<Bio::Map::MapI>

=cut

sub default_map {
	my ($self, $map) = @_;
	if (defined $map) {
		$self->thow("This is [$map], not Bio::Map::MapI object") unless $map->isa('Bio::Map::MapI');
		$self->{'_default_map'} = $map;
	}
	return $self->{'_default_map'} || return;
}

=head2 map

 Title   : map
 Function: This is a synonym of the default_map() method

		   *** does not actually add this marker to the map! ***

 Status  : deprecated, will be removed in next version

=cut

*map = \&default_map;

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
	my ($self, $map, $value) = @_;
	$map ||= $self->default_map;
	if ($value) {
		$self->throw("Value better be scalar, not [$value]") unless ref($value) eq '';
	}

	my $pos = Bio::Map::Position->new();
	$pos->map($map) if $map;
    $pos->value($value) if defined($value);
    $pos->element($self);
	return $pos;
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
    my ($self, $pos, $pos_actual) = @_;

	if ($pos) {
		$self->purge_positions;
		$self->add_position($pos, $pos_actual);
	}

    my @positions = $self->each_position;
    $self->warn('This marker has more than one Position, returning the most recently added') if scalar @positions > 1;
    return pop(@positions);
}

=head2 add_position

 Title   : add_position
 Usage   : $marker->add_position($position);
 Function: Add a Position to this marker
 Returns : n/a
 Args    : Bio::Map::PositionI
            OR
           Bio::Map::MapI AND
           scalar
            OR
           scalar, but only if the marker has a default map

=cut

sub add_position  {
    my ($self, $pos, $pos_actual) = @_;
    $self->throw("Must give a Position") unless defined $pos;

    my $map = $self->default_map;
	my $pos_map;
	if (ref($pos)) {
		if (ref($pos) eq 'ARRAY') {
			($pos, $pos_actual) = @{$pos};
			unless ($pos && $pos_actual && ref($pos)) {
				$self->throw("Supplied an array ref but did not contain two values, the first an object");
			}
		}

		if ($pos->isa('Bio::Map::PositionI')) {
			$pos_map = $pos->map;
			$self->default_map($pos_map) unless $map;
			$map = $pos_map if $pos_map;
		}
		elsif ($pos->isa('Bio::Map::MapI')) {
			$self->default_map($pos) unless $map;
			$map = $pos;
			$pos = $pos_actual;
		}
		else {
			$self->throw("This is [$pos], not a Bio::Map::PositionI or Bio::Map::MapI object");
		}
	}

	$self->throw("You need to give a marker a default map before you can set positions without explicit map!" ) unless $map;

	if (ref($pos) && $pos->isa('Bio::Map::PositionI')) {
		$pos->map($map) unless $pos_map;
		$self->SUPER::add_position($pos);
	}
	else {
		$self->get_position_object($map, $pos); # adds position to us
	}
}

=head2 positions

 Title   : positions
 Usage   : $marker->positions([$pos1, $pos2, $pos3]);
 Function: Add multiple Bio::Map::PositionI to this marker
 Returns : n/a
 Args    : array ref of $map/value tuples or array ref of Bio::Map::PositionI

=cut

sub positions {
    my ($self, $args_ref) = @_;

    foreach my $arg (@{$args_ref}) {
        if (ref($arg) eq 'ARRAY') {
            $self->add_position(@{$arg});
        }
        else {
            $self->add_position($arg);
        }
    }
}

=head2 in_map

 Title   : in_map
 Usage   : if ( $marker->in_map($map) ) {}
 Function: Tests if this marker is found on a specific map
 Returns : boolean
 Args    : a map unique id OR Bio::Map::MapI

=cut

sub in_map {
	my ($self, $query) = @_;
	$self->throw("Must supply an argument") unless defined($query);

	if (ref($query) eq '') {
		foreach my $map ($self->known_maps) {
			my $uid = $map->unique_id;
			if ($uid) {
				($uid eq $query) && return 1;
			}
		}
	}
    else {
		return $self->SUPER::in_map($query);
	}

    return 0;
}

1;
