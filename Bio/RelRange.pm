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
use vars qw( @ISA @EXPORT_OK );

use Bio::Root::Root;
use Bio::RelRangeI;
require Exporter;
@ISA = qw( Bio::Root::Root Bio::RelRangeI Exporter );
@EXPORT_OK = qw( absRange absSeqId absStart absEnd absStrand );

use vars '$VERSION';
$VERSION = '1.00';

use Bio::DB::GFF::Util::Rearrange; # for &rearrange.

=head1 Exported functions

=head2 absRange

  Title   : absRange
  Usage   : use Bio::RelRange qw( absRange );
            my $abs_range = absRange( $range );
  Function: Get the range of the L<Bio::PrimarySeqI> (or whatever)
            that the given L<Bio::RangeI> object is defined over.
  Returns : The root range, or undef if there is none.
  Args    : a L<Bio::RangeI> object.

  This is a utility function provided by the RelRange package to find
  the root range of any L<Bio::RangeI> objects (even those that are not
  L<Bio::RelRangeI> implementers).

=cut

sub absRange {
  my $range = shift;
  if( $range->isa( 'Bio::RelRangeI' ) ) {
    return $range->abs_range();
  }

  ## Okay so it's a RangeI but not a RelRangeI.
  my $seq_id = $range->seq_id();
  my $next_seq_id;
  while( defined( $seq_id ) &&
         ref( $seq_id ) &&
         $seq_id->isa( 'Bio::RangeI' ) ) {
    $next_seq_id = $seq_id->seq_id();
    # If next_seq_id isn't a range then $seq_id *is* the
    # root sequence's range, and that's what we're looking for.
    last unless( defined( $next_seq_id ) &&
                 ref( $next_seq_id ) &&
                 $next_seq_id->isa( 'Bio::RangeI' ) );
    $seq_id = $next_seq_id;
  }
  return $seq_id;
} # absRange()

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
  L<Bio::RelRangeI> implementers).  For convenience this method will
  return the given value if the given value is not a L<Bio::RangeI>
  (this is convenient when the given value is *already* the id of some
  sequence).

=cut

sub absSeqId {
  my $range = shift;
  unless( $range && ref( $range ) && $range->isa( 'Bio::RangeI' ) ) {
    return $range;
  }
  if( $range->isa( 'Bio::RelRangeI' ) ) {
    return $range->abs_seq_id();
  }

  ## Okay so it's a RangeI but not a RelRangeI.
  my $seq_id = $range->seq_id();
  while( defined( $seq_id ) &&
         ref( $seq_id ) &&
         $seq_id->isa( 'Bio::RelRangeI' ) ) {
    $seq_id = $seq_id->seq_id();
  }
  return $seq_id;
} # absSeqId()

=head2 absLow

  Title   : absLow
  Usage   : use Bio::RelRange qw( absLow );
            my $abs_low = absLow( $range );
  Function: Get the absolute start or end position of the given range,
            whichever is lesser-valued.
  Returns : The current absolute start or end position of the given range
            (whichever is lower), relative to absSeqId( $range ).
  Args    : a L<Bio::RangeI> object.

  This is a utility function provided by the RelRange package to convert
  the coordinate of any L<Bio::RangeI> objects (even those that are not
  L<Bio::RelRangeI> implementers) into its absolute equivalent.

=cut

sub absLow {
  my $range = shift;
  unless( $range && ref( $range ) && $range->isa( 'Bio::RangeI' ) ) {
    return undef;
  }
  if( $range->isa( 'Bio::RelRangeI' ) ) {
    return $range->abs_low();
  }

  ## Okay so it's a RangeI but not a RelRangeI.
  my $seq_id = $range->seq_id();
  my $abs_strand = 1;
  my ( $abs_low, $abs_high );
  if( $range->strand() < 0 ) {# then abs_low starts high (confusing, ain't it!)
    $abs_low = $range->_end();
    $abs_high = $range->_start();
    if( $abs_high > $abs_low ) {
      ( $abs_low, $abs_high ) = ( $abs_high, $abs_low );
    }
  } else {
    $abs_low = $range->_start();
    $abs_high = $range->_end();
    if( $abs_high > $abs_low ) {
      ( $abs_low, $abs_high ) = ( $abs_high, $abs_low );
    }
  }

  my ( $next_seq_id, $low, $high );
  while( defined( $seq_id ) &&
         ref( $seq_id ) &&
         $seq_id->isa( 'Bio::RangeI' ) ) {
    $next_seq_id = $seq_id->seq_id();
    $abs_strand *= $seq_id->strand();
    ( $low, $high ) = ( $seq_id->start(), $seq_id->end() );
    if( $low > $high ) {
      ( $low, $high ) = ( $high, $low );
    }

    # If next_seq_id isn't a range then $seq_id *is* the
    # root sequence's range (so it really should be 1..root sequence length).
    if( ( defined( $next_seq_id ) && ref( $next_seq_id ) &&
          $next_seq_id->isa( 'Bio::RangeI' ) ) ||
        ( $low != 1 ) ) {
      if( $abs_strand < 0 ) {
        $abs_low = ( $high - $abs_low + 1 );
        $abs_high = ( $high - $abs_high + 1 );
      } else {
        $abs_low += ( $low - 1 );
        $abs_high += ( $low - 1 );
      }
    }
    $seq_id = $next_seq_id;
  }
  if( $abs_low < $abs_high ) {
    return $abs_low;
  } else {
    return $abs_high;
  }
} # absLow()

=head2 absStart

  Title   : absStart
  Usage   : use Bio::RelRange qw( absStart );
            my $abs_start = absStart( $range );
  Function: Get the absolute start position of the given range.
  Returns : The current start position of the given range, relative to
            absSeqId( $range ).
  Args    : a L<Bio::RangeI> object.

  This is a utility function provided by the RelRange package to
  convert the start coordinate of any L<Bio::RangeI> objects (even
  those that are not L<Bio::RelRangeI> implementers) into their
  absolute equivalents.  If the object is a L<Bio::RelRangeI> then
  this method will delegate to its abs_start method, which uses its
  orientation_policy to determine whether the returned value should be
  the greater or lesser absolute position.  If it is not a RelRange
  then an orientation policy of 'dependent' will be used, so the
  returned value will be the greater value iff the absStrand of the
  given range is negative.

=cut

sub absStart {
  my $range = shift;
  unless( $range && ref( $range ) && $range->isa( 'Bio::RangeI' ) ) {
    return undef;
  }
  if( $range->isa( 'Bio::RelRangeI' ) ) {
    return $range->abs_start();
  }

  ## Okay so it's a RangeI but not a RelRangeI.
  my $abs_strand = absStrand( $range );
  if( $abs_strand < 0 ) {
    return absHigh( $range );
  } else {
    return absLow( $range );
  }
} # absStart()

=head2 absHigh

  Title   : absHigh
  Usage   : use Bio::RelRange qw( absHigh );
            my $abs_high = absHigh( $range );
  Function: Get the absolute high position of the given range.
  Returns : The current high position of the given range, relative to
            absSeqId( $range ).
  Args    : a L<Bio::RangeI> object.

  This is a utility function provided by the RelRange package to convert
  the high coordinate of any L<Bio::RangeI> objects (even those that are not
  L<Bio::RelRangeI> implementers) into their absolute equivalents.

=cut

sub absHigh {
  my $range = shift;
  unless( $range && ref( $range ) && $range->isa( 'Bio::RangeI' ) ) {
    return undef;
  }
  if( $range->isa( 'Bio::RelRangeI' ) ) {
    return $range->abs_low();
  }

  ## Okay so it's a RangeI but not a RelRangeI.
  my $seq_id = $range->seq_id();
  my $abs_strand = 1;
  my ( $abs_low, $abs_high );
  if( $range->strand() < 0 ) {# then abs_high starts low (confusing, ain't it!)
    $abs_low = $range->_end();
    $abs_high = $range->_start();
    if( $abs_high > $abs_low ) {
      ( $abs_low, $abs_high ) = ( $abs_high, $abs_low );
    }
  } else {
    $abs_low = $range->_start();
    $abs_high = $range->_end();
    if( $abs_high > $abs_low ) {
      ( $abs_low, $abs_high ) = ( $abs_high, $abs_low );
    }
  }

  my ( $next_seq_id, $low, $high );
  while( defined( $seq_id ) &&
         ref( $seq_id ) &&
         $seq_id->isa( 'Bio::RangeI' ) ) {
    $next_seq_id = $seq_id->seq_id();
    $abs_strand *= $seq_id->strand();
    ( $low, $high ) = ( $seq_id->start(), $seq_id->end() );
    if( $low > $high ) {
      ( $low, $high ) = ( $high, $low );
    }

    # If next_seq_id isn't a range then $seq_id *is* the
    # root sequence's range (so it really should be 1..root sequence length).
    if( ( defined( $next_seq_id ) && ref( $next_seq_id ) &&
          $next_seq_id->isa( 'Bio::RangeI' ) ) ||
        ( $low != 1 ) ) {
      if( $abs_strand < 0 ) {
        $abs_low = ( $high - $abs_low + 1 );
        $abs_high = ( $high - $abs_high + 1 );
      } else {
        $abs_low += ( $low - 1 );
        $abs_high += ( $low - 1 );
      }
    }
    $seq_id = $next_seq_id;
  }
  if( $abs_low < $abs_high ) {
    return $abs_high;
  } else {
    return $abs_low;
  }
} # absHigh()

=head2 absEnd

  Title   : absEnd
  Usage   : use Bio::RelRange qw( absEnd );
            my $abs_end = absEnd( $range );
  Function: Get the absolute end position of the given range.
  Returns : The current end position of the given range, relative to
            absSeqId( $range ).
  Args    : a L<Bio::RangeI> object.

  This is a utility function provided by the RelRange package to
  convert the end coordinate of any L<Bio::RangeI> objects (even
  those that are not L<Bio::RelRangeI> implementers) into their
  absolute equivalents.  If the object is a L<Bio::RelRangeI> then
  this method will delegate to its abs_end method, which uses its
  orientation_policy to determine whether the returned value should be
  the greater or lesser absolute position.  If it is not a RelRange
  then an orientation policy of 'dependent' will be used, so the
  returned value will be the lesser value iff the absStrand of the
  given range is negative.

=cut

sub absEnd {
  my $range = shift;
  unless( $range && ref( $range ) && $range->isa( 'Bio::RangeI' ) ) {
    return undef;
  }
  if( $range->isa( 'Bio::RelRangeI' ) ) {
    return $range->abs_end();
  }

  ## Okay so it's a RangeI but not a RelRangeI.
  my $abs_strand = absStrand( $range );
  if( $abs_strand < 0 ) {
    return absLow( $range );
  } else {
    return absHigh( $range );
  }
} # absEnd()

=head2 absStrand

  Title   : absStrand
  Usage   : use Bio::RelRange qw( absStrand );
            my $abs_strand = absStrand( $range );
  Function: Get the absolute strandedness (-1, 0, or 1) of the given range.
  Returns : The current absolute strand value of the given range.
  Args    : a L<Bio::RangeI> object.

  This is a utility function provided by the RelRange package to convert
  the strand value of any L<Bio::RangeI> objects (even those that are not
  L<Bio::RelRangeI> implementers) into their absolute equivalents.

=cut

sub absStrand {
  my $range = shift;
  unless( $range && ref( $range ) && $range->isa( 'Bio::RangeI' ) ) {
    return undef;
  }
  if( $range->isa( 'Bio::RelRangeI' ) ) {
    return $range->abs_strand();
  }

  ## Okay so it's a RangeI but not a RelRangeI.
  my $seq_id = $range->seq_id();
  my $abs_strand = $range->strand();

  while( defined( $seq_id ) &&
         ref( $seq_id ) &&
         $seq_id->isa( 'Bio::RangeI' ) ) {
    unless( $abs_strand ) {
      return $abs_strand;
    }
    $abs_strand *= $seq_id->strand();
    $seq_id = $seq_id->seq_id();
  }
  return $abs_strand;
} # absStrand()

=head1 Constructors

=head2 new

  Title   : new
  Usage   : $range = Bio::RelRange->new( $another_range_to_copy_from );
            OR
            $range = Bio::RelRange->new(
                       -seq_id => $another_range_or_a_sequence_id,
                       -start => 100,
                       -end => 200,
                       -strand => +1,
                       -absolute => 1
                     );
  Function: generates a new Bio::RelRange object
  Returns : a new Bio::RelRange object
  Args    : a L<Bio::RangeI> to copy from
            OR
              two of (-start, -end, -length) - the third is calculated
              -strand (defaults to 0)
              -seq_id (not required but highly recommended)
              -absolute (defaults to 0)
              -orientation_policy (defaults to 'ignorant')

    Note that if you pass a RelRangeI that is in absolute() mode, the
    copied values will be absolute, and the relative values will not
    be preserved.  To work around this, set absolute() to false before
    passing it in.  The new RelRangeI will never be in absolute() mode
    when the single-argument form of this method is used.

    If the -absolute argument is true, the other values will be
    interpreted relative to the given -seq_id, but then absolute()
    will be set to true.

    If the given start value is greater than the given end value then
    the strand will be forced to negative.

=cut

sub new {
  my $caller = shift;
  my $self = $caller->SUPER::new( @_ );

  # Hack to deal with the fact that RelRange calls start(), end(), & strand()
  # which will lead to an error in Bio::Search::HSP::BlastHSP 
  # because parsing hasn't yet occurred.
  # TODO: Remove this when BlastHSP doesn't do lazy parsing.
  $self->{ '_initializing' } = 1;

  my ( $seq_id, $strand, $start, $end, $length, $absolute, $policy );

  if( scalar( @_ ) && ( $_[ 0 ] =~ /^-/ ) ) {
    ( $seq_id, $strand, $start, $end, $length, $absolute, $policy ) = 
      rearrange( [ [ qw( SEQ_ID SEQID SEQ ID ) ],
                   qw( STRAND ),
                   [ qw( START BEGIN LOW ) ],
                   [ qw( STOP END HIGH ) ],
                   [ qw( LENGTH LEN ) ],
                   [ qw( ABSOLUTE ABS ) ],
                   [ qw( ORIENTATION_POLICY ORIENTATIONPOLICY
                         ORIENTATION POLICY ) ]
                 ], @_ );
  } else {
    my $copy_from = shift;
    if( $copy_from && ref( $copy_from ) && $copy_from->isa( 'Bio::RangeI' ) ) {
      $seq_id = $copy_from->seq_id();
      $strand = $copy_from->strand();
      $start = $copy_from->start();
      $end = $copy_from->end();
      if( $copy_from->isa( 'Bio::RelRangeI' ) ) {
        $policy = $copy_from->orientation_policy();
        ## Notice that we purposefully don't copy the absoluteness.
      }
    }
  }

  $self->seq_id( $seq_id );

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
  } elsif( defined $end ) {
    $self->end( $end );
    $self->start( 1 );
  } elsif( defined $length ) {
    $self->start( 1 );
    $self->length( $length );
  } else {
    # If there is no start/end/length info, try to inherit the full
    # parent range.
    if( ref( $seq_id ) && $seq_id->isa( 'Bio::RangeI' ) ) {
      $self->start( 1 );
      $self->length( $seq_id->length() );
    }
  }

  ## Force the strand to negative if the coordinates are given in reverse.
  if( $self->start() > $self->end() ) {
    $strand = -1;
  }

  $self->orientation_policy( $policy || 'ignorant' );
  $self->ensure_orientation();

  $self->strand( $strand || 0 );

  $self->absolute( $absolute || 0 );
  return $self;
} # new(..)

=head2 new_from_relrange

 Title   : new_from_relrange
 Usage   : my $new_relrange =
             Bio::RelRange->new_from_relrange( $copy_from );
 Function: Create a new Bio::RelRange object by copying values from
           another RelRange object.
 Returns : A new L<Bio::RelRange> object
 Args    : Another L<Bio::RelRange> object
 Status  : Protected

  This is a special copy constructor.  It forces the new range into
  the L<Bio::RelRange> package, regardless of the package that it is
  called from.  This causes subclass-specfic information to be dropped.

  As a special bonus you may also pass an existing hash and it will be the
  blessed an anointed object that is returned, like so:
    $new_relrange =
      Bio::RelRange->new_from_relrange(
        $copy_from,
        $new_relrange
      );

=cut

sub new_from_relrange {
  my $pack = shift; # ignored
  my $relrange = shift || $pack;
  my $new_relrange = shift || Bio::RelRange->new();
  @{ $new_relrange }{ qw( _absolute _seq_id _start _end _strand ) } =
    @{ $relrange }{ qw( _absolute _seq_id _start _end _strand ) };
  return bless $new_relrange, __PACKAGE__;
} # new_from_relrange(..)

=head1 Bio::RelRangeI methods

These methods implement the L<Bio::RelRangeI> interface (and
consequently the L<Bio::RangeI> interface).

=cut

=head2 orientation_policy

  Title   : orientation_policy
  Usage   : my $orientation_policy =
              $range->orientation_policy( [new_policy] );
  Function: Get/Set the oriention policy that this RelRangeI uses.
  Returns : The current (or former, if used as a set method) value of
            the orientation policy.
  Args    : [optional] A new (string) orientation_policy value.

  The BioPerl community has various opinions about how coordinates
  should be returned when the strand is negative.  Some folks like the
  start to be the lesser-valued position in all circumstances
  ('independent' of the strand value).  Others like the start to be
  the lesser-valued position when the strand is 0 or 1 and the
  greater-valued position when the strand is -1 ('dependent' on the
  strand value).  Others expect that the start and end values are
  whatever they were set to ('ignorant' of the strand value).

  Legal values of orientation_policy are:
      Value          Assertion
   ------------   -------------------
   'independent'  ( start() <= end() )
   'dependent'    (( strand() < 0 )?( end() <= start() ):( start() <= end() ))
   'ignorant'     # No assertion.

  See also ensure_orientation().  Note that changing the
  orientation_policy will not automatically ensure that the
  orientation policy assertion holds, so you should call
  ensure_orientation() also.

=cut

sub orientation_policy {
  my $self = shift;
  my $new_value = shift;
  my $old_value = $self->{ '_orientation_policy' };
  if( defined( $new_value ) ) {
    unless( ( $new_value eq 'independent' ) ||
            ( $new_value eq 'dependent' ) ||
            ( $new_value eq 'ignorant' ) ) {
      $self->throw( "Illegal orientation_policy value: '$new_value'." );
    }
    $self->{ '_orientation_policy' } = $new_value;
  }
  return $old_value;
} # orientation_policy(..)

=head2 ensure_orientation

  Title   : ensure_orientation
  Usage   : $range->ensure_orientation();
  Function: After calling this method, the orientation_policy assertion
            will be true.
  Returns : nothing
  Args    : none

  The orientation_policy is an assertion about the relative values of
  the start() and end() positions.  This assertion might fail when the
  start and end positions change.  This method reorients the values in
  case the assertion fails.  After calling this method the assertion
  will be true.

=cut

sub ensure_orientation {
  my $self = shift;

  # If they care at all about such things, we make sure that the
  # actual stored value of _start is less than that of _end, and then
  # in retrieval we can deal with the 'dependent' case by flipping
  # when the strand is negative.
  if( ( $self->{ '_orientation_policy' } eq 'independent' ) ||
      ( $self->{ '_orientation_policy' } eq 'dependent' ) ) {
    unless( $self->{ '_start' } <= $self->{ '_end' } ) {
      my $tmp = $self->{ '_start' };
      $self->{ '_start' } = $self->{ '_end' };
      $self->{ '_end' } = $tmp;
    }
  }
} # ensure_orientation(..)

#                   --Coders beware!--
# Changes to the Bio::RelRangeI pod need to be copied to here.
#
# This is a final implementation.  That means that its methods are
# defined using each other, so if you override a method you may break
# everything.  You have been warned.

=head2 rel2abs_strand

  Title   : rel2abs_strand
  Usage   : my $abs_strand = $range->rel2abs_strand( $rel_strand );
  Function: Convert a strand that is relative to seq_id() into one that
            is relative to abs_seq_id().
  Returns : a strand value (-1, 0, or 1).
  Args    : a strand value (-1, 0, or 1).

  This function takes a strand value that is relative to seq_id()
  and converts it so that it is absolute (ie. relative to abs_seq_id()).

  Note that if absolute() is true this method still interprets
  the argument strand as it were relative to what seq_id() would
  be if absolute() were false.

=cut

sub rel2abs_strand {
  my $self = shift;
  my $rel_strand = shift;

  if( ( $rel_strand eq '+1' ) || ( $rel_strand eq '+' ) ) {
    $rel_strand = 1;
  } elsif( $rel_strand eq '-' ) {
    $rel_strand = -1;
  } elsif( $rel_strand eq '.' ) {
    $rel_strand = 0;
  }
  unless( ( $rel_strand == -1 ) ||
          ( $rel_strand == 0 ) ||
          ( $rel_strand == 1 ) ) {
    $self->throw( "'$rel_strand' is not in the set { '-1', '0', '1' }." );
  }

  my $ancestor_seq_id = $self->_seq_id();
  my $abs_strand = $rel_strand;

  while( defined( $ancestor_seq_id ) &&
         ref( $ancestor_seq_id ) &&
         $ancestor_seq_id->isa( 'Bio::RangeI' ) ) {
    unless( $abs_strand ) {
      return $abs_strand;
    }
    $abs_strand *= $ancestor_seq_id->strand();
    $ancestor_seq_id = $ancestor_seq_id->seq_id();
  }
  return $abs_strand;
} # rel2abs_strand(..)

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

=head2 abs_range

  Title   : abs_range
  Usage   : my $abs_range = $range->abs_range();
  Function: Get the range of the abs_seq_id that this RangeI is defined over.
  Returns : The root range, or undef if there is none.
  Args    : none

  Ranges may have no defined abs_range, but this should be considered
  deprecated.  The concept of a 'range' requires that it is a range
  over some sequence; this method returns the range of that sequence.
  If the value of seq_id() is a string (the unique_id or primary_id of
  a L<Bio::PrimarySeqI>) then this method will return a range that is
  equal to this one (to $self).  If the value of seq_id() is another
  L<Bio::RangeI>, then this method will return it if its seq_id() if
  is a string, or keep searching up the tree until a range with a
  seq_id that is a string (or undef) is reached, and return that
  range.

=cut

sub abs_range {
  my $self = shift;
  my $ancestor_seq_id = $self->_seq_id();

  unless( defined( $ancestor_seq_id ) &&
          ref( $ancestor_seq_id ) &&
          $ancestor_seq_id->isa( 'Bio::RangeI' ) ) {
    return $self;
  }

  my $next_ancestor_seq_id;
  while( defined( $ancestor_seq_id ) &&
         ref( $ancestor_seq_id ) &&
         $ancestor_seq_id->isa( 'Bio::RangeI' ) ) {
    $next_ancestor_seq_id = $ancestor_seq_id->seq_id();
    # If next_ancestor_seq_id isn't a range then $ancestor_seq_id *is* the
    # root sequence's range, and that's what we're looking for.
    last unless( defined( $next_ancestor_seq_id ) &&
                 ref( $next_ancestor_seq_id ) &&
                 $next_ancestor_seq_id->isa( 'Bio::RangeI' ) );
    $ancestor_seq_id = $next_ancestor_seq_id;
  }
  return $ancestor_seq_id;
} # abs_range()

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
  my $ancestor_seq_id = $self->_seq_id();

  while( defined( $ancestor_seq_id ) &&
         ref( $ancestor_seq_id ) &&
         $ancestor_seq_id->isa( 'Bio::RangeI' ) ) {
    $ancestor_seq_id = $ancestor_seq_id->seq_id();
  }
  return $ancestor_seq_id;
} # abs_seq_id()

=head2 abs_low

  Title   : abs_low
  Usage   : my $abs_low = $range->abs_low();
  Function: Get the least-valued absolute position of this range.
  Returns : The current lowest position of this range, relative to the
            abs_seq_id.
  Args    : none

=cut

sub abs_low {
  my $self = shift;

  my $ancestor_seq_id = $self->_seq_id();
  my $abs_strand = 1;
  my ( $abs_low, $abs_high );
  if( $self->strand() < 0 ) { # then abs_low starts high (confusing, ain't it!)
    $abs_low = $self->_end();
    $abs_high = $self->_start();
    if( $abs_high > $abs_low ) {
      ( $abs_low, $abs_high ) = ( $abs_high, $abs_low );
    }
  } else {
    $abs_low = $self->_start();
    $abs_high = $self->_end();
    if( $abs_high < $abs_low ) {
      ( $abs_low, $abs_high ) = ( $abs_high, $abs_low );
    }
  }

  ## TODO: REMOVE
  #print STDERR "BEGIN abs_low()\n";

  my ( $next_ancestor_seq_id, $low, $high );
  while( defined( $ancestor_seq_id ) &&
         ref( $ancestor_seq_id ) &&
         $ancestor_seq_id->isa( 'Bio::RangeI' )
       ) {
    $next_ancestor_seq_id = $ancestor_seq_id->seq_id();
    $abs_strand *= $ancestor_seq_id->strand();
    ( $low, $high ) = ( $ancestor_seq_id->start(), $ancestor_seq_id->end() );
    if( $low > $high ) {
      ( $low, $high ) = ( $high, $low );
    }

    ## TODO: REMOVE
    #print STDERR "   abs_low is $abs_low\n";
    #print STDERR "   abs_high is $abs_high\n";
    #print STDERR "   low is $low, high is $high\n";
    #print STDERR "   asid is $ancestor_seq_id\n";
    #print STDERR "   abs_strand is $abs_strand.\n";
    #print STDERR "   next_asid is $next_ancestor_seq_id\n";

    # If next_ancestor_seq_id isn't a range then $ancestor_seq_id *is* the
    # root sequence's range (so it really should be 1..root sequence length).
    if( ( defined( $next_ancestor_seq_id ) && ref( $next_ancestor_seq_id ) &&
          $next_ancestor_seq_id->isa( 'Bio::RangeI' ) ) ||
        ( $low != 1 ) ) {
      if( $abs_strand < 0 ) {
        $abs_low = ( $high - $abs_low + 1 );
        $abs_high = ( $high - $abs_high + 1 );
      } else {
        $abs_low += ( $low - 1 );
        $abs_high += ( $low - 1 );
      }
    }
    $ancestor_seq_id = $next_ancestor_seq_id;
  }
  if( $abs_low < $abs_high ) {
    ## TODO: REMOVE
    #print STDERR "END $abs_low\n";
    return $abs_low;
  } else {
    ## TODO: REMOVE
    #print STDERR "END $abs_high\n";
    return $abs_high;
  }
} # abs_low()

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

  This method is defined in terms of abs_high(), abs_low(),
  abs_strand(), and orientation_policy().

=cut

sub abs_start {
  my $self = shift;

  if( $self->orientation_policy() eq 'dependent' ) {
    if( $self->abs_strand() < 0 ) {
      return $self->abs_high();
    } else {
      return $self->abs_low();
    }
  } else {
    return $self->abs_low();
  }
} # abs_start()

=head2 abs_high

  Title   : abs_high
  Usage   : my $abs_high = $range->abs_high();
  Function: Get the greatest-valued absolute position of this range.
  Returns : The current highest position of this range, relative to the
            abs_seq_id.
  Args    : none

=cut

sub abs_high {
  my $self = shift;

  my $ancestor_seq_id = $self->_seq_id();
  my $abs_strand = 1;
  my ( $abs_low, $abs_high );
  if( $self->strand() < 0 ) { # then abs_high starts low (confusing, ain't it!)
    $abs_low = $self->_end();
    $abs_high = $self->_start();
    if( $abs_high > $abs_low ) {
      ( $abs_low, $abs_high ) = ( $abs_high, $abs_low );
    }
  } else {
    $abs_low = $self->_start();
    $abs_high = $self->_end();
    if( $abs_high < $abs_low ) {
      ( $abs_low, $abs_high ) = ( $abs_high, $abs_low );
    }
  }

  ## TODO: REMOVE
  #print STDERR "BEGIN abs_high()\n";

  my ( $next_ancestor_seq_id, $low, $high );
  while( defined( $ancestor_seq_id ) &&
         ref( $ancestor_seq_id ) &&
         $ancestor_seq_id->isa( 'Bio::RangeI' )
       ) {
    $next_ancestor_seq_id = $ancestor_seq_id->seq_id();
    $abs_strand *= $ancestor_seq_id->strand();
    ( $low, $high ) = ( $ancestor_seq_id->start(), $ancestor_seq_id->end() );
    if( $low > $high ) {
      ( $low, $high ) = ( $high, $low );
    }

    ## TODO: REMOVE
    #print STDERR "   abs_high is $abs_high\n";
    #print STDERR "   abs_low is $abs_low\n";
    #print STDERR "   high is $high, low is $low\n";
    #print STDERR "   asid is $ancestor_seq_id\n";
    #print STDERR "   abs_strand is $abs_strand.\n";
    #print STDERR "   next_asid is $next_ancestor_seq_id\n";

    # If next_ancestor_seq_id isn't a range then $ancestor_seq_id *is* the
    # root sequence's range (so it really should be 1..root sequence length).
    if( ( defined( $next_ancestor_seq_id ) && ref( $next_ancestor_seq_id ) &&
          $next_ancestor_seq_id->isa( 'Bio::RangeI' ) ) ||
        ( $low != 1 ) ) {
      if( $abs_strand < 0 ) {
        $abs_low = ( $high - $abs_low + 1 );
        $abs_high = ( $high - $abs_high + 1 );
      } else {
        $abs_low += ( $low - 1 );
        $abs_high += ( $low - 1 );
      }
    }
    $ancestor_seq_id = $next_ancestor_seq_id;
  }
  if( $abs_low > $abs_high ) {
    ## TODO: REMOVE
    #print STDERR "END $abs_low\n";
    return $abs_low;
  } else {
    ## TODO: REMOVE
    #print STDERR "END $abs_high\n";
    return $abs_high;
  }
} # abs_high()

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

  This method is defined in terms of abs_high(), abs_low(),
  abs_strand(), and orientation_policy().

=cut

sub abs_end {
  my $self = shift;

  if( $self->orientation_policy() eq 'dependent' ) {
    if( $self->abs_strand() < 0 ) {
      return $self->abs_low();
    } else {
      return $self->abs_high();
    }
  } else {
    return $self->abs_high();
  }
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

  my $ancestor_seq_id = $self->_seq_id();
  my $abs_strand = $self->_strand();

  while( defined( $ancestor_seq_id ) &&
         ref( $ancestor_seq_id ) &&
         $ancestor_seq_id->isa( 'Bio::RangeI' )
       ) {
    unless( $abs_strand ) {
      return $abs_strand;
    }
    $abs_strand *= $ancestor_seq_id->strand();
    $ancestor_seq_id = $ancestor_seq_id->seq_id();
  }
  return $abs_strand;
} # abs_strand()

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

=head2 abs2rel_strand

  Title   : abs2rel_strand
  Usage   : my $rel_strand = $range->abs2rel_strand( $abs_strand )
  Function: Convert a strand that is relative to abs_seq_id() into one that
            is relative to seq_id().
  Returns : a strand value (-1, 0, or 1).
  Args    : a strand value (-1, 0, or 1).

  This function takes a strand value that is absolute (ie. relative to
  abs_seq_id()) and converts it so that it is relative to seq_id().

  Note that if absolute() is true this method still returns the strand
  relative to what seq_id() would be if absolute() were false.

  This method turns out to be identical to rel2abs_strand, so it is
  implemented in the interface as an inheritable alias for
  rel2abs_strand.

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

## Note: We pass through any additional args to _seq_id(..), but only
## when setting the new value.
sub seq_id {
  my $self = shift;
  my ( $new_val, @other_args ) = @_;
  my $old_val =
    ( $self->absolute() ?
      $self->abs_seq_id() :
      $self->_seq_id() );
  if( defined( $new_val ) ) {
    if( $self->absolute() ) {
      $self->throw( "Unable to set the absolute seq_id.  If you would like to set the relative seq_id, first turn the absolute() flag off, then try again." );
    }
    $self->_seq_id( $new_val, @other_args );
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
      $self->_start() );
  if( defined( $new_val ) ) {
    if( $self->absolute() ) {
      $self->throw( "Unable to set the absolute start value.  If you would like to set the relative start, first turn the absolute() flag off, then try again." );
    }
    unless( $new_val =~ /^[-+]?\d+$/ ) {
      $self->throw( "'$new_val' is not an integer.\n" );
    }
    $self->_start( $new_val );
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
      $self->_end() );
  if( defined( $new_val ) ) {
    if( $self->absolute() ) {
      $self->throw( "Unable to set the absolute end value.  If you would like to set the relative end, first turn the absolute() flag off, then try again." );
    }
    unless( $new_val =~ /^[-+]?\d+$/ ) {
      $self->throw( "'$new_val' is not an integer.\n" );
    }
    $self->_end( $new_val );
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
      $self->_strand() );
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
    $self->_strand( $new_val );
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

  length = ( ( abs_high - abs_low ) + 1 ).

  Note the interdependence of start()|end()|abs_start()|abs_end() and
  length().  Changing start() or end() will change the length.
  Changing the length will change the end() (and consequently abs_end()).

  You may not set the length when absolute() is true.

=cut

sub length {
  my $self = shift;
  my $old_val = ( ( $self->high() - $self->low() ) + 1 );
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

# Internal overridable getter/setter for the actual stored value of seq_id.
sub _seq_id {
  my $self = shift;
  my ( $new_val ) = @_;
  my $old_val = $self->{ '_seq_id' };
  if( defined $new_val ) {
    $self->{ '_seq_id' } = $new_val;
  }
  return $old_val;
} # _seq_id(..)

# Internal overridable getter/setter for the actual stored value of start.
# Uses the '_end' value if $self->_strand is negative, just to be confusing.
sub _start {
  my $self = shift;
  my ( $new_val ) = @_;
  my $old_val =
    ( ( ( $self->{ '_orientation_policy' } eq 'dependent' ) &&
        ( $self->{ '_strand' } < 0 ) ) ?
      $self->{ '_end' } :
      $self->{ '_start' } );
  if( defined $new_val ) {
    ( ( ( $self->{ '_orientation_policy' } eq 'dependent' ) &&
        ( $self->{ '_strand' } < 0 ) ) ?
      $self->{ '_end' } :
      $self->{ '_start' } ) = $new_val;
  }
  return $old_val;
} # _start(..)

# Internal overridable getter/setter for the actual stored value of end.
# Uses the '_start' value if $self->_strand is negative, just to be confusing.
sub _end {
  my $self = shift;
  my ( $new_val ) = @_;
  my $old_val =
    ( ( ( $self->{ '_orientation_policy' } eq 'dependent' ) &&
        ( $self->{ '_strand' } < 0 ) ) ?
      $self->{ '_start' } :
      $self->{ '_end' } );
  if( defined $new_val ) {
    ( ( ( $self->{ '_orientation_policy' } eq 'dependent' ) &&
        ( $self->{ '_strand' } < 0 ) ) ?
      $self->{ '_start' } :
      $self->{ '_end' } ) = $new_val;
  }
  return $old_val;
} # _end(..)

# Internal overridable getter/setter for the actual stored value of strand.
sub _strand {
  my $self = shift;
  my ( $new_val ) = @_;
  my $old_val = $self->{ '_strand' };
  if( defined $new_val ) {
    $self->{ '_strand' } = $new_val;
  }
  return $old_val;
} # _strand(..)

1;

__END__
