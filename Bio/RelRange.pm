package Bio::RelRange;

# $Id $
# A simple implementation of Bio::RelRangeI, a Bio::RangeI with
# additional methods to support shifting between relative and absolute
# views.

=head1 NAME

Bio::RelRange -- A simple, final implementation of the L<Bio::RelRangeI>
interface, which is a L<Bio::RangeI> with additional methods to
support shifting between relative and absolute views.

=head1 SYNOPSIS

=head1 DESCRIPTION

A L<Bio::RangeI> is a range over a sequence, and may be defined
relative to another range.  This class implements the
L<Bio::RelRangeI> interface, which provides additional methods for
accessing the range in relative and absolute coordinate spaces.

WARNING: This is a final implementation.  That means that its methods
are defined using each other, so if you override a method you may
break everything.  You may feel free to subclass it, so long as you
don't override the methods in it.

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

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...
use strict;
use vars qw( @ISA );

use Bio::Root::Root;
use Bio::RelRangeI;
use Exporter;
@ISA = qw( Bio::Root::Root Bio::RelRangeI Exporter );
@EXPORT_OK = qw( absSeqId );

use vars '$VERSION';
$VERSION = '1.00';

=head1 Exported functions

=head2 absSeqId

  Title   : absSeqId
  Usage   : use Bio::RelRange qw( absSeqId );
            my $abs_seq_id = absSeqId( $range );
  Function: Get the unique_id or primary_id of the L<Bio::PrimarySeqI>
            that the given L<Bio::RangeI> object is defined over.
  Returns : The root seq_id, or undef if there is none.
  Args    : a L<Bio::RangeI> object.

  This is a utility function provided by the RelRange package to find
  the root seq_id of any L<Bio::RangeI> objects (even those that are not
  L<Bio::RelRangeI> implementers).

=cut

sub absSeqId {
  my $range = shift;
  unless( ref( $range ) && $range->isa( 'Bio::RangeI' ) ) {
    return undef;
  }
  if( $range->isa( 'Bio::RelRangeI' ) ) {
    return $range->abs_seq_id();
  }

  ## Okay so it's a RangeI but not a RelRangeI.
  my $seq_id = $range->seq_id();
  while( defined( $seq_id ) &&
         ref( $seq_id ) ) {
    ## Assertion: ( $seq_id isa 'Bio::RangeI' )
    $seq_id = $seq_id->seq_id();
  }
  return $seq_id;
} # absSeqId()

=head1 Constructors

=head2 new

  Title   : new
  Usage   : $range = Bio::RelRange->new(
                       -seq_id => $another_range_or_a_sequence_id,
                       -start => 100,
                       -end => 200,
                       -strand => +1,
                       -absolute => 1
                     );
  Function: generates a new Bio::RelRange object
  Returns : a new Bio::RelRange object
  Args    : two of (-start, -end, -length) - the third is calculated
          : -strand (defaults to 0)
          : -seq_id (not required but highly recommended)
          : -absolute (defaults to 0)

=cut

sub new {
  my ( $caller, @args ) = @_;
  my $self = $caller->SUPER::new( @args );
  my ( $seq_id, $strand, $start, $end, $length, $absolute ) = 
    $self->_rearrange( [ qw( SEQ_ID
                             STRAND
                             START
                             END
                             LENGTH
                             ABSOLUTE
                           ) ], @args );
  $self->seq_id( $seq_id );
  $self->absolute( $absolute || 0);
  $self->strand($strand || 0);

  if( defined $start ) {
    $self->start( $start );
    if( defined $end ) {
      $self->end( $end );
    } elsif( defined $length ) {
      $self->length( $length );
    }
  } elsif( defined $end && defined $length ) {
    $self->end( $end );
    $self->start( $self->end() - $length + 1 );
  }
  return $self;
} # new(..)

=head1 Bio::RelRangeI methods

These methods implement the L<Bio::RelRangeI> interface (and
consequently the L<Bio::RangeI> interface).

=cut

#                   --Coders beware!--
# Changes to the Bio::RangeI pod need to be copied to here.
#
# This is a final implementation.  That means that its methods are
# defined using each other, so if you override a method you may break
# everything.  You have been warned.

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
  my $self = shift;
  my ( $new_val ) = @_;
  my $old_val = $self->{ '_absolute' };
  if( defined( $new_val ) ) {
    $self->{ '_absolute' } = $new_val;
  }
  return $old_val;
} # absolute(..)

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
  my $self = shift;
  my $ancestor_seq_id = $self->{ '_seq_id' };

  while( defined( $ancestor_seq_id ) &&
         ref( $ancestor_seq_id ) ) {
    ## Assertion: ( $ancestor_seq_id isa 'Bio::RangeI' )
    $ancestor_seq_id = $ancestor_seq_id->seq_id();
  }
  return $ancestor_seq_id;
} # abs_seq_id()

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
  my $self = shift;

  my $ancestor_seq_id = $self->{ '_seq_id' };
  my $abs_strand = $self->{ '_strand' };
  my $abs_start = $self->{ '_start' };

  while( defined( $ancestor_seq_id ) &&
         ref( $ancestor_seq_id ) ) {
    ## Assertion: ( $ancestor_seq_id isa 'Bio::RangeI' )
    $ancestor_seq_id = $ancestor_seq_id->seq_id();
    $abs_strand *= $ancestor_seq_id->strand();
    if( $abs_strand < 0 ) {
      $abs_start = ( $ancestor_seq_id->start() - $abs_start + 1 );
    } else {
      $abs_start += ( $ancestor_seq_id->start() - 1 );
    }
  }
  return $abs_start;
} # abs_start()

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
  my $self = shift;

  my $ancestor_seq_id = $self->{ '_seq_id' };
  my $abs_strand = $self->{ '_strand' };
  my $abs_end = $self->{ '_end' };

  while( defined( $ancestor_seq_id ) &&
         ref( $ancestor_seq_id ) ) {
    ## Assertion: ( $ancestor_seq_id isa 'Bio::RangeI' )
    $ancestor_seq_id = $ancestor_seq_id->seq_id();
    $abs_strand *= $ancestor_seq_id->strand();
    if( $abs_strand < 0 ) {
      $abs_end = ( $ancestor_seq_id->end() - $abs_end + 1 );
    } else {
      $abs_end += ( $ancestor_seq_id->end() - 1 );
    }
  }
  return $abs_end;
} # abs_end()

=head2 abs_strand

  Title   : abs_strand
  Usage   : my $abs_strand = $range->abs_strand();
  Function: Get the absolute strandedness (-1, 0, or 1) of this range.
  Returns : The current absolute strand value of this range.
  Args    : none

=cut

sub abs_strand {
  my $self = shift;

  my $ancestor_seq_id = $self->{ '_seq_id' };
  my $abs_strand = $self->{ '_strand' };

  while( defined( $ancestor_seq_id ) &&
         ref( $ancestor_seq_id ) ) {
    unless( $abs_strand ) {
      return $abs_strand;
    }
    ## Assertion: ( $ancestor_seq_id isa 'Bio::RangeI' )
    $ancestor_seq_id = $ancestor_seq_id->seq_id();
    $abs_strand *= $ancestor_seq_id->strand();
  }
  return $abs_strand;
} # abs_strand()

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

=head2 rel2abs

  Title   : rel2abs
  Usage   : my @abs_coords = $range->rel2abs( @rel_coords );
  Function: Convert relative coordinates into absolute coordinates
  Returns : a list of absolute coordinates
  Args    : a list of relative coordinates

  This function takes a list of positions in relative coordinates
  (relative to seq_id()), and converts them into absolute coordinates.
  Note that if absolute() is true this method does nothing.

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

=head2 abs2rel

  Title   : abs2rel
  Usage   : my @rel_coords = $range->abs2rel( @abs_coords )
  Function: Convert absolute coordinates into relative coordinates
  Returns : a list of relative coordinates
  Args    : a list of absolute coordinates

  This function takes a list of positions in absolute coordinates
  and converts them into relative coordinates (relative to seq_id()).
  Note that if absolute() is true this method does nothing.

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

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

=cut

sub seq_id {
  my $self = shift;
  my ( $new_val ) = @_;
  my $old_val =
    ( $self->absolute() ?
      $self->abs_seq_id() :
      $self->{ '_seq_id' } );
  if( defined( $new_val ) ) {
    if( $self->absolute() ) {
      $self->throw( "Unable to set the absolute seq_id.  If you would like to set the relative seq_id, first turn the absolute() flag off, then try again." );
    }
    if( ( ( ref $new_val ) && !$new_val->isa( 'Bio::RangeI' ) ) ||
        ( ref \$new_val ne 'STRING' ) ) {
      $self->throw( "The given value is neither a string nor a Bio::RangeI object." );
    }
    $self->{ '_seq_id' } = $new_val;
  }
  return $old_val;
} # seq_id(..)

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

=cut

sub start {
  my $self = shift;
  my ( $new_val ) = @_;
  my $old_val =
    ( $self->absolute() ?
      $self->abs_start() :
      $self->{ '_start' } );
  if( defined( $new_val ) ) {
    if( $self->absolute() ) {
      $self->throw( "Unable to set the absolute start value.  If you would like to set the relative start, first turn the absolute() flag off, then try again." );
    }
    unless( $new_val =~ /^[-+]?\d+$/ ) {
      $self->throw( "'$new_val' is not an integer.\n" );
    }
    $self->{ '_start' } = $new_val;
  }
  return $old_val;
} # start(..)

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

=cut

sub end {
  my $self = shift;
  my ( $new_val ) = @_;
  my $old_val =
    ( $self->absolute() ?
      $self->abs_end() :
      $self->{ '_end' } );
  if( defined( $new_val ) ) {
    if( $self->absolute() ) {
      $self->throw( "Unable to set the absolute end value.  If you would like to set the relative end, first turn the absolute() flag off, then try again." );
    }
    unless( $new_val =~ /^[-+]?\d+$/ ) {
      $self->throw( "'$new_val' is not an integer.\n" );
    }
    $self->{ '_end' } = $new_val;
  }
  return $old_val;
} # end(..)

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

=cut

sub strand {
  my $self = shift;
  my ( $new_val ) = @_;
  my $old_val =
    ( $self->absolute() ?
      $self->abs_strand() :
      $self->{ '_strand' } );
  if( defined( $new_val ) ) {
    if( $self->absolute() ) {
      $self->throw( "Unable to set the absolute strand value.  If you would like to set the relative strand, first turn the absolute() flag off, then try again." );
    }
    if( ( $new_val eq '+1' ) || ( $new_val eq '+' ) ) {
      $new_val = 1;
    } elsif( $new_val eq '-' ) {
      $new_val = -1;
    } elsif( $new_val eq '.' ) {
      $new_val = 0;
    }
    unless( ( $new_val == -1 ) || ( $new_val == 0 ) || ( $new_val == 1 ) ) {
      $self->throw( "'$new_val' is not in the set { '-1', '0', '1' }." );
    }
    $self->{ '_strand' } = $new_val;
  }
  return $old_val;
} # strand(..)

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

=cut

sub length {
  my $self = shift;
  # Note that it shouldn't matter whether we use absolute or relative
  # positions, but we differentiate anyway just to be on the safe side.
  my $old_val =
    ( $self->absolute() ?
      ( ( $self->end() - $self->start() ) + 1 ) :
      ( ( $self->abs_high() - $self->abs_low() ) + 1 ) );
  my ( $new_val ) = @_;

  if( defined( $new_val ) ) {
    if( $self->absolute() ) {
      $self->throw( "Unable to set the length value on a range that has its absolute() flag set.  If you would like to set the length, first turn the absolute() flag off, then try again." );
    }
    unless( $new_val ) {
      $self->throw( "Ranges must have a non-zero length." );
    }
    $self->end( $self->start() + ( $new_val - 1 ) );
  }
  return $old_val;
} # length(..)

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

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

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

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

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

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

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

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

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

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

=head2 overlap_extent

 Title   : overlap_extent
 Usage   : my ( $a_unique, $common, $b_unique ) = $a->overlap_extent( $b );
 Function: Provides actual amount of overlap between two different ranges.
 Returns : 3-tuple consisting of:
           - the number of positions unique to a
           - the number of positions common to both
           - the number of positions unique to b
 Args    : a L<Bio::RangeI> object

  This method is implemented in the interface, and need not be
  overridden in concrete subclasses.

=cut

1;

__END__
