# BioPerl module for Bio::Map::OrderedPosition
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

Bio::Map::OrderedPosition - Abstracts the notion of a member
	of an ordered list of markers. Each marker is a certain distance
	from the one in the ordered list before it.

=head1 SYNOPSIS

    use Bio::Map::OrderedPosition;
	# the first marker in the sequence
    my $position = Bio::Map::OrderedPosition->new(-order => 1,
			-positions => [ [ $map, 22.3] ] );
	# the second marker in the sequence, 15.6 units from the fist one
    my $position2 = Bio::Map::OrderedPosition->new(-order => 2,
			-positions => [ [ $map, 37.9] ] );
	# the third marker in the sequence, coincidental with the second
	# marker
    my $position3 = Bio::Map::OrderedPosition->new(-order => 3,
                        -posititions => [ [ $map, 37.9]] );

=head1 DESCRIPTION

This object is an implementation of the PositionI interface and the
Position object handles the specific values of a position.
OrderedPosition is intended to be slightly more specific then Position
but only specific enough for a parser from the MarkerIO subsystem to
create and then pass to a client application to bless into the proper
type. For an example of how this is intended to work, see the
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
Jason Stajich, jason@bioperl.org
Sendu Bala, bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::OrderedPosition;
use strict;


use base qw(Bio::Map::Position);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Map::OrderedPosition->new();
 Function: Builds a new Bio::Map::OrderedPosition object 
 Returns : Bio::Map::OrderedPosition
 Args    : -order : The order of this position

=cut

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    
    my ($order) = $self->_rearrange([qw(ORDER)], @args);
    $order && $self->order($order);
    
    return $self;
}

=head2 order

 Title   : order
 Usage   : $o_position->order($new_order);
           my $order = $o_position->order();
 Function: Get/set the order position of this position in a map.
 Returns : int, the order of this position
 Args    : none to get, OR int to set

=cut

sub order {
    my ($self, $order) = @_;
    if ($order) {
        $self->{'_order'} = $order;
    }
    return $self->{'_order'} || return;
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
    return $self->order;
}

=head2 equals

 Title   : equals
 Usage   : if ($mappable->equals($mapable2)) {...}
 Function: Test if a position is equal to another position.
 Returns : boolean
 Args    : Bio::Map::PositionI

=cut

sub equals {
   my ($self,$compare) = @_;
   return 0 if (! defined $compare || ! $compare->isa('Bio::Map::OrderedPosition'));
   return ($compare->order == $self->order);
}

# admittedly these aren't really the best comparisons in the world
# but it is a first pass we'll need to refine the algorithm or not 
# provide general comparisions and require these to be implemented
# by objects closer to the specific type of data

=head2 less_than

 Title   : less_than
 Usage   : if ($mappable->less_than($m2)) {...}
 Function: Tests if a position is less than another position
           It is assumed that 2 positions are in the same map.
 Returns : boolean
 Args    : Bio::Map::PositionI

=cut

sub less_than {
   my ($self,$compare) = @_;
   return 0 if (! defined $compare || ! $compare->isa('Bio::Map::OrderedPosition'));
   return ($compare->order < $self->order);
}

=head2 greater_than

 Title   : greater_than
 Usage   : if ($mappable->greater_than($m2)) {...}
 Function: Tests if position is greater than another position.
           It is assumed that 2 positions are in the same map.
 Returns : boolean
 Args    : Bio::Map::PositionI

=cut

sub greater_than {
   my ($self,$compare) = @_;
   return 0 if (! defined $compare || ! $compare->isa('Bio::Map::OrderedPosition'));
   return ($compare->order > $self->order);
}

1;