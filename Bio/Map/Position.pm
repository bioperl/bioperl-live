#
# BioPerl module for Bio::Map::Position
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

Bio::Map::Position - A single position of a Marker, or the range over which
                     that marker lies, in a Map

=head1 SYNOPSIS

    use Bio::Map::Position;
    my $position = Bio::Map::Position->new(-map => $map, 
					  -element => $marker,
					  -value => 100
					  );

	my $position_with_range = Bio::Map::Position->new(-map => $map, 
					  -element => $marker,
					  -start => 100,
					  -length => 10
					  );

=head1 DESCRIPTION

This object is an implementation of the PositionI interface that
handles the specific values of a position. This allows a map element
(e.g. Marker) to have multiple positions within a map and still be
treated as a single entity.

This handles the concept of a relative map in which the order of
elements and the distance between them is known, but does not
directly handle the case when distances are unknown - in that case
arbitrary values must be assigned for position values.

No units are assumed here - units are handled by context of which Map
a position is placed in or the subclass of this Position.

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
use strict;

use Scalar::Util qw(looks_like_number);
use Bio::Map::Relative;

use base qw(Bio::Root::Root Bio::Map::PositionI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Map::Position->new();
 Function: Builds a new Bio::Map::Position object 
 Returns : Bio::Map::Position
 Args    : -map      => Bio::Map::MapI object
           -element  => Bio::Map::MappableI object
           -relative => Bio::Map::RelativeI object

           * If this position has no range, or if a single value can describe
             the range *
           -value => scalar             : something that describes the single
                                          point position or range of this
                                          Position, most likely an int

           * Or if this position has a range, at least two of *
           -start => int                : value of the start co-ordinate
           -end => int                  : value of the end co-ordinate
           -length => int               : length of the range

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
	
    my ($map, $marker, $element, $value, $start, $end, $length, $relative) = 
	$self->_rearrange([qw( MAP
			       MARKER
                   ELEMENT
			       VALUE
				   START
				   END
				   LENGTH
                   RELATIVE
			       )], @args);
	
    my $do_range = defined($start) || defined($end);
    if ($value && $do_range) {
        $self->warn("-value and (-start|-end|-length) are mutually exclusive, ignoring value");
		$value = undef;
    }
	
    $map            && $self->map($map);
    $marker         && $self->element($marker); # backwards compatibility
    $element        && $self->element($element);
    $relative       && $self->relative($relative);
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
    my ($self, $relative) = @_;
    if ($relative) {
        $self->throw("Must supply an object") unless ref($relative);
        $self->throw("This is [$relative], not a Bio::Map::RelativeI") unless $relative->isa('Bio::Map::RelativeI');
        $self->{_relative_not_implicit} = 1;
        $self->{_relative} = $relative;
    }
    return $self->{_relative} || $self->absolute_relative;
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
    if (@_) { $self->{_absolute} = shift }
    return $self->{_absolute} || 0;
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
 Returns : scalar numeric
 Args    : none to get the co-ordinate normally (see absolute() method), OR
           Bio::Map::RelativeI to get the co-ordinate converted to be
           relative to what this Relative describes.

=cut

sub numeric {
    my ($self, $value) = @_;
    my $num = $self->{'_value'};
    $self->throw("The value has not been set, can't convert to numeric") unless defined($num);
    $self->throw("This value [$num] is not numeric") unless looks_like_number($num);
    
    if (ref($value) && $value->isa('Bio::Map::RelativeI')) {
        # get the value after co-ordinate conversion
        my $raw = $num;
        my ($abs_start, $rel_start) = $self->_relative_handler($value);
        return $abs_start + $raw - $rel_start;
    }
    
    # get the value as per absolute
    if ($self->{_relative_not_implicit} && $self->absolute) {
        # this actually returns the start, but should be the same thing...
        return $self->relative->absolute_conversion($self);
    }
    
    return $num;
}

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

sub start {
	my ($self, $value) = @_;
    if (defined $value) {
        if (ref($value) && $value->isa('Bio::Map::RelativeI')) {
            # get the value after co-ordinate conversion
            my $raw = $self->{start};
            defined $raw || return;
            my ($abs_start, $rel_start) = $self->_relative_handler($value);
            return $abs_start + $raw - $rel_start;
        }
        else {
            # set the value
            $self->throw("This is [$value], not a number") unless looks_like_number($value);
            $self->{start} = $value;
            $self->value($value) unless defined($self->value);
        }
    }
    
    # get the value as per absolute
    if ($self->{_relative_not_implicit} && $self->absolute) {
        return $self->relative->absolute_conversion($self);
    }
    
    return defined($self->{start}) ? $self->{start} : return;
}

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

sub end {
	my ($self, $value) = @_;
    if (defined $value) {
        if (ref($value) && $value->isa('Bio::Map::RelativeI')) {
            # get the value after co-ordinate conversion
            my $raw = $self->{end};
            defined $raw || return;
            my ($abs_start, $rel_start) = $self->_relative_handler($value);
            return $abs_start + $raw - $rel_start;
        }
        else {
            # set the value
            $self->throw("This value [$value] is not numeric!") unless looks_like_number($value);
            $self->{end} = $value;
        }
    }
    
    # get the value as per absolute
    if ($self->{_relative_not_implicit} && $self->absolute) {
        my $raw = $self->{end} || return;
        my $abs_start = $self->relative->absolute_conversion($self) || return;
        return $abs_start + ($raw - $self->{start});
    }
    
    return defined($self->{end}) ? $self->{end} : return;
}

=head2 length

  Title   : length
  Usage   : $length = $position->length();
  Function: Get/set the length of this position's range, changing the end() if
            necessary. Getting and even setting the length will fail if both
            start() and end() are not already defined.
  Returns : the length of this range
  Args    : none to get, OR scalar numeric (>0) to set.

=cut

sub length {
	my ($self, $length) = @_;
	if ($length) {
        $length > 0 || return;
		my $existing_length = $self->length || return;
        return $length if $existing_length == $length;
        $self->end($self->{start} + $length - 1);
	}
	
	if (defined($self->start) && defined($self->end)) {
		return $self->end - $self->start + 1;
	}
    return;
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
    my ($self, $given_map) = @_;
    my $answer = $self->numeric($self->absolute_relative);
    return $answer;
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
	my ($self, $rel) = @_;
	if (defined($self->start) && defined($self->end)) {
		return $self->start($rel).'..'.$self->end($rel);
	}
	return '';
}

=head2 absolute_relative

 Title   : absolute_relative
 Usage   : my $rel = $position->absolute_relative();
 Function: Get a relative describing the start of the map. This is useful for
           supplying to the coordinate methods (start(), end() etc.) to get
           the temporary effect of having set absolute(1).
 Returns : Bio::Map::Relative
 Args    : none

=cut

sub absolute_relative {
    return Bio::Map::Relative->new(-map => 0, -description => 'start of map');
}

# get our own absolute start and that of the thing we want as a frame of
# reference
sub _relative_handler {
    my ($self, $value) = @_;
    
    my $own_relative = $self->relative;
    
    # if the requested relative position is the same as the actual
    # relative, the current co-ordinate values are correct so shortcut
    my ($own_type, $req_type) = ($own_relative->type, $value->type);
    if ($own_type && $req_type && $own_type eq $req_type && $own_relative->$own_type eq $value->$req_type) {
        return (0, 0);
    }
    
    my $abs_start = $own_relative->absolute_conversion($self);
    my $rel_start = $value->absolute_conversion($self);
    $self->throw("Unable to resolve co-ordinate because relative to something that ultimately isn't relative to the map start")
    unless defined($abs_start) && defined($rel_start);
    
    return ($abs_start, $rel_start);
}

1;
