### TO BE DELETED ###

# BioPerl module for Bio::Map::OrderedPositionWithDistance
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chad Matsalla <bioinformatics1@dieselwurks.com>
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::OrderedPositionWithDistance - Abstracts the notion of a member
	of an ordered list of markers. Each marker is a certain distance
	from the one in the ordered list before it.

=head1 SYNOPSIS

    use Bio::Map::OrderedPositionWithDistance;
	# the first marker in the sequence
    my $position = Bio::Map::OrderedPositionWithDistance->new(-positions => 1,
			-distance => 22.3 );
	# the second marker in the sequence, 15.6 units from the fist one
    my $position2 = Bio::Map::OrderedPositionWithDistance->new(-positions => 2,
			-distance => 15.6 );
	# the third marker in the sequence, coincidental with the second
	# marker
    my $position3 = Bio::Map::OrderedPositionWithDistance->new(-positions => 3,
			-distance => 0 );


=head1 DESCRIPTION

This object is an implementation of the PositionI interface and the
Position object handles the specific values of a position.
OrderedPositionWithDistance is intended to be slightly more specific
then Position but only specific enough for a parser from the MarkerIO
subsystem to create and then pass to a client application to bless into
the proper type. For an example of how this is intended to work, see the
Mapmaker.pm.

No units are assumed here - units are handled by context of which Map
a position is placed in.

Se Bio::Map::Position for additional information.

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

=head1 AUTHOR - Chad Matsalla

Email bioinformatics1@dieselwurks.com

=head1 CONTRIBUTORS

Lincoln Stein, lstein@cshl.org
Heikki Lehvaslaiho, heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Map::OrderedPositionWithDistance;
use strict;


use base qw(Bio::Map::Position);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Map::OrderedPositionWithDistance->new();
 Function: Builds a new Bio::Map::OrderedPositionWithDistance object 
 Returns : Bio::Map::OrderedPositionWithDistance
 Args    : -positions  - Should be a single value representing the order
	of this marker within the list of markers
	-distance - The distance this marker is from the marker before it.
		0 reflects coincidentality.

=cut

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  $self->{'_positions'} = [];
  my ($positions,$distance) = $self->_rearrange([qw(POSITIONS DISTANCE)], @args);
  if( ref($positions) =~ /array/i ) { 
      foreach my $p ( @$positions ) {
	  $self->add_position($p);
      }
  } else { 
      $self->add_position($positions);
  }
	$distance && $self->distance($distance);
	
  return $self;

}


=head2 distance($new_distance)

 Title   : distance($new_distance)
 Usage   : $position->distance(new_distance) _or_
        $position->distance()
 Function: get/set the distance of this position from the previous marker
 Returns : A scalar representing the current distance for this position.
 Args    : If $new_distance is provided the distance of this Position will
        be set to $new_distance

=cut

sub distance {
        my ($self,$distance) = @_;
        if ($distance) {
           $self->{'_distance'} = $distance;
        }
        return $self->{'_distance'};
}


1;
