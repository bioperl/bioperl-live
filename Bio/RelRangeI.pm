package Bio::RelRangeI;

# $Id $
# A Bio::RangeI with additional methods to support shifting between
# relative and absolute views.

=head1 NAME

Bio::RelRangeI -- A Bio::RangeI with additional methods to support
shifting between relative and absolute views.

=head1 SYNOPSIS

=head1 DESCRIPTION

A L<Bio::RangeI> is a range over a sequence, and may be defined
relative to another range.  This interface, L<Bio::RelRangeI>,
provides additional methods for accessing the range in relative and
absolute coordinate spaces.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR

Paul Edlefsen E<lt>paul@systemsbiology.orgE<gt>.
Copyright (c) 2003 Institute for Systems Biology

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=head1 CONTRIBUTORS

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...
use strict;
use vars qw( @ISA );

use Bio::RangeI;
@ISA = qw( Bio::RangeI );

use vars '$VERSION';
$VERSION = '1.00';

=head1 Bio::RelRangeI methods

These methods are unique to Bio::RelRangeI (that is, they are not
inherited from above).

=cut

=head2 absolute

  Title   : absolute
  Usage   : my $absolute_flag = $range->absolute( [$new_absolute_flag] );
  Function: Get/set the absolute flag.
  Returns : The current (or former, if used as a set method) value of the
            absolute flag.
  Args    : [optional] a new value for the absolute flag.

  If the absolute() flag is set then the start(), end(), and strand()
  methods will behave like the abs_start(), abs_end(), and abs_strand()
  methods, meaning that they will return values relative to abs_seq_id()
  rather than to seq_id().

=cut

sub absolute {
  shift->throw_not_implemented();
}

=head2 abs_seq_id

  Title   : abs_seq_id
  Usage   : my $abs_seq_id = $range->abs_seq_id();
  Function: Get the unique_id or primary_id of the L<Bio::PrimarySeqI>
            that this RangeI is defined over.
  Returns : The root seq_id, or undef if there is none.
  Args    : none

  Ranges may have no defined abs_seq_id, but this should be considered
  deprecated.  The concept of a 'range' requires that it is a range
  over some sequence; this method returns that sequence.  If the value
  of seq_id() is a string (the unique_id or primary_id of a
  L<Bio::PrimarySeqI>) then this method will be identical to seq_id().
  If the value of seq_id() is another L<Bio::RangeI>, then this method
  will return its seq_id() if that is a string, or keep searching up the
  tree until a string (or undef) is reached.

=cut

sub abs_seq_id {
  shift->throw_not_implemented();
}

=head2 abs_start

  Title   : abs_start
  Usage   : my $abs_start = $range->abs_start();
  Function: Get the absolute start position of this range.
  Returns : The current start position of this range, relative to the
            abs_seq_id.
  Args    : none

  Note the interdependence of abs_start() and start().  Changing start() will
  change abs_start().

  Note the interdependence of abs_start() and length().  Changing length() will
  change abs_start().

=cut

sub abs_start {
  shift->throw_not_implemented();
}

=head2 abs_end

  Title   : abs_end
  Usage   : my $abs_end = $range->abs_end();
  Function: Get the absolute end position of this range.
  Returns : The current absolute end position of this range, relative
            to the abs_seq_id.
  Args    : none

  Note the interdependence of abs_end() and end().  Changing end() will
  change abs_end().

  Note the interdependence of abs_end() and length().  Changing length() will
  change abs_end().

=cut

sub abs_end {
  shift->throw_not_implemented();
}

=head2 abs_strand

  Title   : abs_strand
  Usage   : my $abs_strand = $range->abs_strand();
  Function: Get the absolute strandedness (-1, 0, or 1) of this range.
  Returns : The current absolute strand value of this range.
  Args    : none

=cut

sub abs_strand {
  shift->throw_not_implemented();
}

=head2 abs_low

  Title   : abs_low
  Usage   : my $abs_low = $range->abs_low();
  Function: Get the least-valued absolute position of this range.
  Returns : The current lowest position of this range, relative to the
            abs_seq_id.
  Args    : none

  This will return either abs_start() or abs_end(), depending on which
  is lower.

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

=cut

sub abs_low {
  my $self = shift;
  my ( $a, $b ) = ( $self->abs_start(), $self->abs_end() );
  return ( ( $a < $b ) ? $a : $b );
} # abs_low()

=head2 abs_high

  Title   : abs_high
  Usage   : my $abs_high = $range->abs_high();
  Function: Get the greatest-valued absolute position of this range.
  Returns : The current highest position of this range, relative to the
            abs_seq_id.
  Args    : none

  This will return either abs_start() or abs_end(), depending on which
  is higher.

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

=cut

sub abs_high {
  my $self = shift;
  my ( $a, $b ) = ( $self->abs_start(), $self->abs_end() );
  return ( ( $a > $b ) ? $a : $b );
} # abs_high()

=head2 rel2abs

  Title   : rel2abs
  Usage   : my @abs_coords = $range->rel2abs( @rel_coords );
  Function: Convert relative coordinates into absolute coordinates
  Returns : a list of absolute coordinates
  Args    : a list of relative coordinates

  This function takes a list of positions in relative coordinates
  (relative to seq_id()), and converts them into absolute coordinates.

  Note that if absolute() is true this method still interprets
  incoming coordinates as if they were relative to what seq_id() would
  be if absolute() were false.

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.  Note that this implementation
  uses abs_start() and abs_strand(), so these methods should not be
  defined in terms of rel2abs(), lest a vicious cycle occur.

=cut

sub rel2abs {
  my $self = shift;

  my @result;

  my $abs_start = $self->abs_start();
  my $abs_strand = $self->abs_strand();
  @result = (
             ( $abs_strand < 0 ) ?
             map { $abs_start - $_ + 1 } @_ :
             map { $_ + $abs_start - 1 } @_
            );

  # if called with a single argument, caller will expect a single scalar reply
  # not the size of the returned array!
  return $result[ 0 ] if ( ( @result == 1 ) and !wantarray );
  return @result;
} # rel2abs(..)

=head2 abs2rel

  Title   : abs2rel
  Usage   : my @rel_coords = $range->abs2rel( @abs_coords )
  Function: Convert absolute coordinates into relative coordinates
  Returns : a list of relative coordinates
  Args    : a list of absolute coordinates

  This function takes a list of positions in absolute coordinates
  and converts them into relative coordinates (relative to seq_id()).

  Note that if absolute() is true this method still produces
  coordinates relative to what seq_id() would be if absolute() were
  false.

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.  Note that this implementation
  uses abs_start() and abs_strand(), so these methods should not be
  defined in terms of abs2rel(), lest a vicious cycle occur.

=cut

sub abs2rel {
  my $self = shift;
  my @result;

  my $abs_start = $self->abs_start();
  my $abs_strand = $self->abs_strand();
  @result = (
             ( $abs_strand < 0 ) ?
             map { $abs_start - $_ + 1 } @_ :
             map { $_ - $abs_start + 1 } @_
            );

  # if called with a single argument, caller will expect a single scalar reply
  # not the size of the returned array!
  return $result[ 0 ] if ( ( @result == 1 ) and !wantarray );
  return @result;
} # abs2rel(..)

=head1 Bio::RangeI methods

These methods are inherited from L<Bio::RangeI>.  Changes between this
interface and that one are noted in the pod, but to be sure you might
want to check L<Bio::RangeI> in case they have gotten out of sync.

=cut

#                   --Coders beware!--
# Changes to the Bio::RangeI pod need to be copied to here.

=head2 seq_id

  Title   : seq_id
  Usage   : my $seq_id = $range->seq_id( [new_seq_id] );
  Function: Get/Set a unique_id or primary_id of a L<Bio::PrimarySeqI>
            or another L<Bio::RangeI> that this RangeI is defined
            over or relative to.  If absolute() is true, this will be
            identical to abs_seq_id().
  Returns : The current (or former, if used as a set method) value of
            the seq_id.
  Args    : [optional] A new (string or L<Bio::RangeI> seq_id value

  Ranges may have no defined seq_id, but this should be considered
  deprecated.  The concept of a 'range' requires that it is a range
  over some sequence; this method returns (and optionally sets) that
  sequence.  It is also possible to specify another range, to support
  relative ranges.  If the value of seq_id is another L<Bio::RangeI>,
  then this RelRangeI's positions are relative to that RangeI's
  positions (unless absolute() is true, in which case they are
  relative to the root seq_id).  If seq_id is the id of a sequence then
  it should provide enough information for a user of a RelRangeI to
  retrieve that sequence; ideally it should be a
  L<Bio::GloballyIdentifiableI> unique_id.

  You may not set the seq_id when absolute() is true.

  The behavior of this method differs from its behavior in
  L<Bio::RangeI> when absolute() is true.

=cut

sub seq_id {
  shift->throw_not_implemented();
}

=head2 start

  Title   : start
  Usage   : my $start = $range->start( [$new_start] );
  Function: Get/set the start of this range.
  Returns : The current (or former, if used as a set method) start position
            of this range.  If absolute() is true then this value will
            be relative to the abs_seq_id; otherwise it will be
            relative to the seq_id.
  Args    : [optional] a new start position

  Note the interdependence of start() and abs_start().  Changing start() will
  change abs_start().

  You may not set start() when absolute() is true.

  Note the interdependence of start() and length().  Changing start() will
  change length().

  The behavior of this method differs from its behavior in
  L<Bio::RangeI> when absolute() is true.  Also, the RangeI interface
  does not specify whether the return value should be the new value or
  the old value.  This interface specifies that it should be the old
  value.

=cut

sub start {
  shift->throw_not_implemented();
}

=head2 end

  Title   : end
  Usage   : my $end = $range->end( [$new_end] );
  Function: Get/set the end of this range.
  Returns : The current (or former, if used as a set method) end position
            of this range.  If absolute() is true then this value will
            be relative to the abs_seq_id; otherwise it will be
            relative to the seq_id.
  Args    : [optional] a new end position

  Note the interdependence of end() and abs_end().  Changing end() will
  change abs_end().

  You may not set end() when absolute() is true.

  Note the interdependence of end() and length().  Changing one will
  change the other.

  The behavior of this method differs from its behavior in
  L<Bio::RangeI> when absolute() is true.  Also, the RangeI interface
  does not specify whether the return value should be the new value or
  the old value.  This interface specifies that it should be the old
  value.

=cut

sub end {
  shift->throw_not_implemented();
}

=head2 strand

  Title   : strand
  Usage   : my $strand = $range->strand( [$new_strand] );
  Function: Get/set the strandedness (-1, 0, or 1) of this range.
  Returns : The current (or former, if used as a set method) strand value
            of this range.  If absolute() is true then this value will
            be absolute.  Otherwise it will be relative to the
            strandedness (if any) of seq_id.
  Args    : [optional] a new strand value.

  You may not set strand() when absolute() is true.

  The behavior of this method differs from its behavior in
  L<Bio::RangeI> when absolute() is true.  Also, the RangeI interface
  does not specify whether the return value should be the new value or
  the old value.  This interface specifies that it should be the old
  value.

=cut

sub strand {
  shift->throw_not_implemented();
}

=head2 length

  Title   : length
  Usage   : my $length = $range->length( [$new_length] );
  Function: Get/set the length of this range.
  Returns : The current (or former, if used as a set method) length
            of this range.
  Args    : [optional] a new length

  length = ( ( end - start ) + 1 ) = ( ( abs_high - abs_low ) + 1 ).

  Note the interdependence of start()|end()|abs_start()|abs_end() and
  length().  Changing start() or end() will change the length.
  Changing the length will change the end() (and consequently abs_end()).

  You may not set the length when absolute() is true.

  This method differs from L<Bio::RangeI> in that it here accepts an
  argument to modify the length.  If the length is modified then the
  return value will be the former length.

=cut

sub length {
  shift->throw_not_implemented();
}

=head2 overlaps

  Title   : overlaps
  Usage   : if( $r1->overlaps( $r2 ) ) { do stuff }
  Function: tests if $r2 overlaps $r1
  Args    : arg #1 = a L<Bio::RangeI> to compare this one to (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') (optional)
  Returns : true if the ranges overlap, false otherwise

  The second argument's values may be:
   "strong"        ranges must have the same strand
   "weak"          ranges must have the same strand or no strand
   "ignore"        ignore strand information (default)

  The behavior of this method differs from its behavior in
  L<Bio::RangeI> when the other range is also a RelRangeI.  When both
  are RelRangeI objects, the abs_seq_ids must be the same and
  absolute coordinates are always used for the test.  If either range
  has no defined abs_seq_id then abs_seq_id will be ignored in the test.

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

=cut

sub overlaps {
  my $self = shift;
  my ( $other, $strand_option ) = @_;
  if( defined( $other ) && $other->isa( 'Bio::RelRangeI' ) ) {
    return (
            $self->_testStrand( $other, $strand_option ) and
            ( ( defined( $self->abs_seq_id() ) &&
                defined( $other->abs_seq_id() ) ) ?
              ( $self->abs_seq_id() eq $other->abs_seq_id() ) : 1 ) and
            not (
                 ( $self->abs_low() > $other->abs_high() or
                   $self->abs_high() < $other->abs_low() ) )
           );
  } else {
    return $self->SUPER::overlaps(@_);
  }
} # overlaps(..)

=head2 contains

  Title   : contains
  Usage   : if( $r1->contains( $r2 ) ) { do stuff }
  Function: tests if $r2 is totally contained within $r1
  Args    : arg #1 = a L<Bio::RangeI> to compare this one to,
                     or an integer position (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') (optional)
  Returns : true iff this range wholly contains the given range

  The second argument's values may be:
   "strong"        ranges must have the same strand
   "weak"          ranges must have the same strand or no strand
   "ignore"        ignore strand information (default)

  The behavior of this method differs from its behavior in
  L<Bio::RangeI> when the other range is also a RelRangeI.  When both
  are RelRangeI objects, the abs_seq_ids must be the same and
  absolute coordinates are always used for the test.  If either range
  has no defined abs_seq_id then abs_seq_id will be ignored in the test.

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

=cut

sub contains {
  my $self = shift;
  my ( $other, $strand_option ) = @_;
  if( defined( $other ) && $other->isa( 'Bio::RelRangeI' ) ) {
    return (
            $self->_testStrand( $other, $strand_option ) and
            ( ( defined( $self->abs_seq_id() ) &&
                defined( $other->abs_seq_id() ) ) ?
              ( $self->abs_seq_id() eq $other->abs_seq_id() ) : 1 ) and
            ( $self->abs_low() >= $other->abs_low() ) and
            ( $self->abs_high() <= $other->abs_high() )
           );
  } else {
    return $self->SUPER::contains(@_);
  }
} # contains(..)

=head2 equals

  Title   : equals
  Usage   : if( $r1->equals( $r2 ) ) { do something }
  Function: Test whether $r1 has the same abs_start, abs_end, length,
            and abs_seq_id as $r2.
  Args    : arg #1 = a L<Bio::RangeI> to compare this one to (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') (optional)
  Returns : true iff they are describing the same range

  If either range has no defined (abs) seq_id then (abs) seq_id will
  be ignored in the test.

  The behavior of this method differs from its behavior in
  L<Bio::RangeI> when the other range is a RelRangeI.  When both
  are RelRangeI objects, the abs_seq_ids must be the same (instead of
  the seq_ids) and absolute coordinates are always used for the test.

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

=cut

sub equals {
  my $self = shift;
  my ( $other, $strand_option ) = @_;
  if( defined( $other ) && $other->isa( 'Bio::RelRangeI' ) ) {
    return (
            $self->_testStrand( $other, $strand_option ) and
            ( ( defined( $self->abs_seq_id() ) &&
                defined( $other->abs_seq_id() ) ) ?
              ( $self->abs_seq_id() eq $other->abs_seq_id() ) : 1 ) and
            ( $self->abs_low() == $other->abs_low() ) and
            ( $self->abs_high() == $other->abs_high() )
           );
  } else {
    return $self->SUPER::equals( @_ );
  }
} # equals(..)

=head2 intersection

  Title   : intersection
  Usage   : my $intersection_range = $r1->intersection( $r2 ) (scalar context)
            OR
            my ( $start, $end, $strand ) = $r1->intersection( $r2 )
             (list context)
  Function: gives the range that is contained by both ranges
  Args    : arg #1 = a range to compare this one to (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') (optional)
  Returns : undef if they do not overlap,
            or new range object containing the overlap
            or (in list context) the start, end, and strand of that range.

  The behavior of this method differs from its behavior in
  L<Bio::RangeI> when the other range is also a RelRangeI.  When both
  are RelRangeI objects abs_seq_ids must not be different and
  absolute coordinates are always used.  The returned object will have
  as its seq_id the abs_seq_id of this range (or, if that is undef, the
  abs_seq_id of the other range).

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

=cut

sub intersection {
  my $self = shift;
  my ( $other, $strand_option ) = @_;

  if( $other->isa( 'Bio::RelRangeI' ) {
    unless( $self->_testStrand( $other, $strand_option ) &&
            ( ( defined( $self->abs_seq_id() ) &&
                defined( $other->abs_seq_id() ) ) ?
              ( $self->abs_seq_id() eq $other->abs_seq_id() ) : 1 ) ) {
      return undef;
    }

    my @low = sort { $a <=> $b } ( $self->abs_low(), $other->abs_low() );
    my @high   = sort { $a <=> $b } ( $self->abs_high(), $other->abs_high() );

    my $low = pop @low;
    my $high = shift @high;
    if( $low > $high ) {
      return undef;
    }

    if( wantarray ) {
      return ( $low, $high, ( ( $self->strand() == $other->strand() ) ?
                              $self->strand() : 0 ) );
    }
    return $self->new( '-seq_id' =>
                         ( defined( $self->abs_seq_id() ) ?
                           $self->abs_seq_id() :
                           $other->abs_seq_id() ),
                       '-start' => $low,
                       '-end' => $high,
                       '-strand' =>
                         ( ( $self->strand() == $other->strand() ) ?
                           $self->strand() :
                           0 )
                     );
  } else {
    return $self->SUPER::intersection( @_ );
  }
} # intersection(..)

=head2 union

  Title   : union
  Usage   : my $union_range = $r1->union( @other_ranges ); (scalar context)
            OR
            my ( $start, $end, $strand ) = $r1->union( @other_ranges );
              (list context)
            OR
            my $union_range = Bio::RelRangeI->union( @ranges );
              (scalar context)
            OR
            my ( $start, $end, $strand ) = Bio::RelRangeI->union( @ranges );
              (list context)
  Function: finds the minimal range that contains all of the ranges
  Args    : a range or list of ranges to find the union of
  Returns : a new range object that contains all of the given ranges, or
            (in list context) the start, end, and strand of that range object.

  The behavior of this method differs from its behavior in
  L<Bio::RangeI> when another range is also a RelRangeI.  When both
  are RelRangeI objects abs_seq_ids must not be different and
  absolute coordinates are always used.  RangeIs may be given, mixed
  with RelRangeIs, but they will be treated as if they were RelRangeIs
  with absolute() set to true.  The returned object will have
  as its seq_id the abs_seq_id of this range (or, if that is undef, the
  abs_seq_id of the first other range with one defined).

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

=cut

sub union {
  my $self = shift;
  my @ranges = @_;
  if( ref $self ) {
    unshift @ranges, $self;
  }
  my $abs_seq_id = $self->abs_seq_id();
  my $union_strand = $self->strand();

  my ( $low, $high );
  foreach my $other ( @ranges ) {
    if( defined( $union_strand ) ) {
      if( $union_strand != $other->strand() ) {
        $union_strand = 0;
      }
    } else {
      $union_strand = $other->strand();
    }
    if( $other->isa( 'Bio::RelRangeI' ) ) {
      if( defined( $abs_seq_id ) &&
          defined( $other->abs_seq_id() ) &&
          ( $abs_seq_id ne $other->abs_seq_id() )
        ) {
        $self->throw( "At least one of the given RelRangeI objects has an incompatible abs_seq_id() value." );
      }
      unless( defined( $abs_seq_id ) ) {
        $abs_seq_id = $other->abs_seq_id();
      }
      if( !defined( $low ) or ( $low > $other->abs_low() ) ) {
        $low = $other->abs_low();
      }
      if( !defined( $high ) or ( $high < $other->abs_high() ) ) {
        $high = $other->abs_high();
      }
    } else { # It's not a RelRangeI; must be a RangeI.  Assume absoluteness.
      if( defined( $abs_seq_id ) &&
          defined( $other->seq_id() ) &&
          ( $abs_seq_id ne $other->seq_id() )
        ) {
        $self->throw( "At least one of the given RangeI objects has an incompatible seq_id() value.  It must be the same as the abs_seq_id() value of this RelRangeI object." );
      }
      unless( defined( $abs_seq_id ) ) {
        $abs_seq_id = $other->seq_id();
      }
      if( !defined( $low ) or ( $low > $other->start() ) ) {
        $low = $other->start();
      }
      if( !defined( $high ) or ( $high < $other->end() ) ) {
        $high = $other->end();
      }
    }
  }
  if( wantarray ) {
    return ( $low, $high, $union_strand );
  }
  return $self->new( -seq_id => $abs_seq_id,
	             -start  => $low,
	             -stop   => $high,
                     -strand => $union_strand
                   );
} # union(..)

=head2 overlap_extent

 Title   : overlap_extent
 Usage   : my ( $a_unique, $common, $b_unique ) = $a->overlap_extent( $b );
 Function: Provides actual amount of overlap between two different ranges.
 Returns : 3-tuple consisting of:
           - the number of positions unique to a
           - the number of positions common to both
           - the number of positions unique to b
 Args    : a L<Bio::RangeI> object

  The behavior of this method differs from its behavior in
  L<Bio::RangeI> when the other range is also a RelRangeI.  When both
  are RelRangeI objects abs_seq_ids must not be different and
  absolute coordinates are always used.

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

=cut

sub overlap_extent {
  my $self = shift;
  my ( $other ) = @_;

  if( $other->isa( 'Bio::RelRangeI' ) ) {
    unless( $self->overlaps( $other ) ) {
      return ( $self->length(), 0 , $other->length() );
    }

    my ( $self_unique, $other_unique );
    if( $self->abs_low() < $other->abs_low() ) {
      $self_unique = $other->abs_low() - $self->abs_low();
    } else {
      $other_unique = $self->abs_low() - $other->abs_low();
    }

    if( $self->abs_high() > $other->abs_high() ) {
      $self_unique += $self->abs_high() - $other->abs_high();
    } else {
      $other_unique += $other->abs_high() - $self->abs_high();
    }
    my $intersection = $self->intersection( $other );
    
    return (
            $self_unique,
            ( 1 + $intersection->end() - $intersection->start() ),
            $other_unique
           );
  } else {
    return $self->SUPER::overlap_extent( @_ );
  }
} # overlap_extent(..)

1;

__END__
