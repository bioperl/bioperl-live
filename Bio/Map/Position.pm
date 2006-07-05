# $Id$
#
# BioPerl module for Bio::Map::Position
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::Position - A single position of a Marker, or the range over which
                     that marker lies, in a Map

=head1 SYNOPSIS

    use Bio::Map::Position;
    my $position = new Bio::Map::Position(-map => $map, 
					  -marker => $marker,
					  -value => 100
					  );
					  
	my $position_with_range = new Bio::Map::Position(-map => $map, 
					  -marker => $marker,
					  -start => 100,
					  -length => 10
					  );

=head1 DESCRIPTION

This object is an implementation of the PositionI interface that
handles the specific values of a position.  This allows an element
(e.g. Marker) to have multiple positions within a map and still be
treated as a single entity.

This does not directly handle the concept of a relative map in which
no known exact positions exist but markers are just ordered relative
to one another - in that case arbitrary values must be assigned for
position values.

No units are assumed here - units are handled by context of which Map
a position is placed in or the subclass of this Position.

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
Chad Matsalla, bioinformatics1@dieselwurks.com
Sendu Bala, bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::Position;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::Map::PositionI;
use Scalar::Util qw(looks_like_number); # comes with perl 5.8, included in Bundle::Bioperl

@ISA = qw(Bio::Root::Root Bio::Map::PositionI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Map::Position();
 Function: Builds a new Bio::Map::Position object 
 Returns : Bio::Map::Position
 Args    : -map     a <Bio::Map::MapI> object
           -marker  a <Bio::Map::MarkerI> object
           * If this position has no range *
           -value => scalar             : something that describes the single
                                          point position of this Position, most
                                          likely an int
           
           * Or if this position has a range, at least two of *
           -start => int                : value of the start co-ordinate
           -end => int                  : value of the end co-ordinate
           -length => int               : length of the range

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
	
    my ($map, $marker, $value, $start, $end, $length) = 
	$self->_rearrange([qw( MAP 
			       MARKER 
			       VALUE
				   START
				   END
				   LENGTH
			       )], @args);
	
    my $do_range = defined($start) || defined($end);
    if ($value && $do_range) {
        $self->warn("-value and (-start|-end|-length) are mutually exclusive, ignoring value");
		$value = undef;
    }
	
    $map            && $self->map($map);
    $marker         && $self->marker($marker);
    defined($value) && $self->value($value);
	
    if ($do_range) {
		defined($start) && $self->start($start);
		defined($end)   && $self->end($end);
		if ($length) {
			if (defined($start) && ! defined($end)) {
				$self->end($start + $length - 1);
			}
			elsif (! defined($start)) {
				$self->start($end - $length + 1);
			}
		}
		defined($self->end) || $self->end($start);
    }
	
    return $self;
}

=head2 value

 Title   : value
 Usage   : my $pos = $position->value;
 Function: Get/Set the value for this postion
 Returns : scalar, value
 Args    : [optional] new value to set

=cut

sub value {
	my ($self, $value) = @_;
	if (defined $value) {
		$self->{'_value'} = $value;
		$self->start($self->numeric) unless defined($self->start);
		$self->end($self->numeric) unless defined($self->end);
	}
	return $self->{'_value'};
}

=head2 numeric

 Title   : numeric
 Usage   : my $num = $position->numeric;
 Function: Read-only method that is guaranteed to return a numeric 
           representation of the start of this position. 
 Returns : numeric (int or real) 
 Args    : none

=cut

sub numeric {
    my ($self) = @_;
    my $num = $self->{'_value'} || 0;
    
    $self->throw("This value [$num] is not numeric!") unless looks_like_number($num);
    
    return $num;
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
    return $self->numeric;
}

=head2 start

  Title   : start
  Usage   : $start = $position->start();
  Function: get/set the start of this range
  Returns : the start of this range
  Args    : optionally allows the start to be set
            using $range->start($start)

=cut

sub start {
	my ($self, $value) = @_;
    if (defined $value) {
		self->throw("This value [$value] is not numeric!") unless looks_like_number($value);
        $self->{'start'} = $value;
		$self->value($value) unless defined($self->value);
    }
    return $self->{'start'} || undef;
}

=head2 end

  Title   : end
  Usage   : $end = $position->end();
  Function: get/set the end of this range
  Returns : the end of this range
  Args    : optionally allows the end to be set
            using $range->end($end)

=cut

sub end {
	my ($self, $value) = @_;
    if (defined $value) {
		self->throw("This value [$value] is not numeric!") unless looks_like_number($value);
        $self->{'end'} = $value;
    }
    return $self->{'end'} || undef;
}

=head2 length

  Title   : length
  Usage   : $length = $position->length();
  Function: get the length of this range
  Returns : the length of this range
  Args    : none

=cut

sub length {
	my $self = shift;
	if (@_) {
		$self->warn(ref($self)."->length() is read-only");
	}
	
	if (defined($self->start) && defined($self->end)) {
		return $self->end - $self->start + 1;
	}
}

=head2 toString

  Title   : toString
  Usage   : print $position->toString(), "\n";
  Function: stringifies this range
  Returns : a string representation of the range of this Position

=cut

sub toString {
	my $self = shift;
	if (defined($self->start) && defined($self->end)) {
		return $self->start.'..'.$self->end;
	}
	return '';
}

1;
