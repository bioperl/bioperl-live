# $Id$
#
# BioPerl module for Bio::RangeI
#
# Cared for by Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Matthew Pocock
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::RangeI - Range interface

=head1 SYNOPSIS

  #Do not run this module directly

=head1 DESCRIPTION

This provides a standard BioPerl range interface that should be
implemented by any object that wants to be treated as a range. This
serves purely as an abstract base class for implementers and can not
be instantiated.

Ranges are modeled as having (start, end, length, strand). They use
Bio-coordinates - all points E<gt>= start and E<lt>= end are within the
range. End is always greater-than or equal-to start, and length is
greather than or equal to 1. The behaviour of a range is undefined if
ranges with negative numbers or zero are used.

So, in summary:

  length = end - start + 1
  end >= start
  strand = (-1 | 0 | +1)

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                         - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki@ebi.ac.uk

=head1 CONTRIBUTORS

Juha Muilu (muilu@ebi.ac.uk)

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::RangeI;

use strict;
use Carp;
use Bio::Root::RootI;
use vars qw(@ISA);
use integer;
use vars qw( @ISA %STRAND_OPTIONS );

@ISA = qw( Bio::Root::RootI );

BEGIN {
# STRAND_OPTIONS contains the legal values for the strand options
    %STRAND_OPTIONS = map { $_, '_'.$_ }
    (
     'strong', # ranges must have the same strand
     'weak',   # ranges must have the same strand or no strand
     'ignore', # ignore strand information
     );
}

# utility methods
#

# returns true if strands are equal and non-zero
sub _strong {
    my ($r1, $r2) = @_;
    my ($s1, $s2) = ($r1->strand(), $r2->strand());
    
    return 1 if $s1 != 0 && $s1 == $s2;
}

# returns true if strands are equal or either is zero
sub _weak {
    my ($r1, $r2) = @_;
    my ($s1, $s2) = ($r1->strand(), $r2->strand());
    return 1 if $s1 == 0 || $s2 == 0 || $s1 == $s2;
}

# returns true for any strandedness
sub _ignore {
    return 1;
}

# works out what test to use for the strictness and returns true/false
# e.g. $r1->_testStrand($r2, 'strong')
sub _testStrand() {
    my ($r1, $r2, $comp) = @_;
    return 1 unless $comp;
    my $func = $STRAND_OPTIONS{$comp};
    return $r1->$func($r2);
}

=head1 Abstract methods

These methods must be implemented in all subclasses.

=head2 seq_id

  Title   : seq_id
  Usage   : my $seq_id = $range->seq_id( [new_seq_id] );
  Function: Get/Set a unique_id or primary_id of a L<Bio::PrimarySeqI>
            or another L<Bio::RangeI> that this RangeI is defined
            over or relative to.
  Returns : The current (or former, if used as a set method) value of
            the seq_id.
  Args    : [optional] A new (string or L<Bio::RangeI> seq_id value

  Ranges may have no defined seq_id, but this should be considered
  deprecated.  The concept of a 'range' requires that it is a range
  over some sequence; this method returns (and optionally sets) that
  sequence.  It is also possible to specify another range, to support
  relative ranges.  If the value of seq_id is another L<Bio::RangeI>,
  then this RangeI's positions are relative to that RangeI's
  positions.  If seq_id is the id of a sequence then it should provide
  enough information for a user of a RangeI to retrieve that sequence;
  ideally it should be a L<Bio::GloballyIdentifiableI> unique_id.

=cut

sub seq_id {
    shift->throw_not_implemented();
}

=head2 start

  Title   : start
  Usage   : $start = $range->start();
  Function: get/set the start of this range
  Returns : the start of this range
  Args    : optionaly allows the start to be set
            using $range->start($start)

=cut

sub start {
    shift->throw_not_implemented();
}

=head2 end

  Title   : end
  Usage   : $end = $range->end();
  Function: get/set the end of this range
  Returns : the end of this range
  Args    : optionaly allows the end to be set
            using $range->end($end)

=cut

sub end {
    shift->throw_not_implemented();
}

=head2 length

  Title   : length
  Usage   : $length = $range->length();
  Function: get/set the length of this range
  Returns : the length of this range
  Args    : optionaly allows the length to be set
             using $range->length($length)

=cut

sub length {
    shift->throw_not_implemented();
}

=head2 strand

  Title   : strand
  Usage   : $strand = $range->strand();
  Function: get/set the strand of this range
  Returns : the strandidness (-1, 0, +1)
  Args    : optionaly allows the strand to be set
            using $range->strand($strand)

=cut

sub strand {
    shift->throw_not_implemented();
}

=head1 Boolean Methods

These methods return true or false. They throw an error if start and
end are not defined.

  $range->overlaps($otherRange) && print "Ranges overlap\n";

=head2 overlaps

  Title   : overlaps
  Usage   : if($r1->overlaps($r2)) { do stuff }
  Function: tests if $r2 overlaps $r1
  Args    : arg #1 = a range to compare this one to (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') (optional)
  Returns : true if the ranges overlap, false otherwise

=cut

sub overlaps {
    my ($self, $other, $so) = @_;
    
    $self->throw("start is undefined") unless defined $self->start;
    $self->throw("end is undefined") unless defined $self->end;
    $self->throw("not a Bio::RangeI object") unless defined $other && 
	$other->isa('Bio::RangeI');
    $other->throw("start is undefined") unless defined $other->start;
    $other->throw("end is undefined") unless defined $other->end;
    
    return
	($self->_testStrand($other, $so) 
	 and not (
		  ($self->start() > $other->end() or
		   $self->end() < $other->start()   )
		  ));
}

=head2 contains

  Title   : contains
  Usage   : if($r1->contains($r2) { do stuff }
  Function: tests whether $r1 totally contains $r2 
  Args    : arg #1 = a range to compare this one to (mandatory)
	             alternatively, integer scalar to test
            arg #2 = strand option ('strong', 'weak', 'ignore') (optional)
  Returns : true if the argument is totaly contained within this range

=cut

sub contains {
  my ($self, $other, $so) = @_;
  $self->throw("start is undefined") unless defined $self->start;
  $self->throw("end is undefined") unless defined $self->end;

  if(defined $other && ref $other) { # a range object?
      $other->throw("Not a Bio::RangeI object") unless  $other->isa('Bio::RangeI');
      $other->throw("start is undefined") unless defined $other->start;
      $other->throw("end is undefined") unless defined $other->end;

      return ($self->_testStrand($other, $so)      and
	      $other->start() >= $self->start() and
	      $other->end() <= $self->end());
  } else { # a scalar?
      $self->throw("'$other' is not an integer.\n") unless $other =~ /^[-+]?\d+$/;
      return ($other >= $self->start() and $other <= $self->end());
  }
}

=head2 equals

  Title   : equals
  Usage   : if( $r1->equals( $r2 ) ) { do something }
  Function: Test whether $r1 has the same start, end, length, and seq_id as $r2
  Args    : a range to test for equality
  Returns : true if they are describing the same range

  If either range has no defined seq_id then seq_id will be ignored in
  the test.

=cut

sub equals {
    my ($self, $other, $so) = @_;

    $self->throw("start is undefined") unless defined $self->start;
    $self->throw("end is undefined") unless defined $self->end;
    $other->throw("Not a Bio::RangeI object") unless  $other->isa('Bio::RangeI');
    $other->throw("start is undefined") unless defined $other->start;
    $other->throw("end is undefined") unless defined $other->end;

    return ($self->_testStrand($other, $so)   and
            ( ( defined( $self->seq_id() ) && defined( $other->seq_id() ) ) ?
              ( $self->seq_id() eq $other->seq_id() ) : 1 ) and
	    $self->start() == $other->start() and
	    $self->end()   == $other->end() );
}

=head1 Geometrical methods

These methods do things to the geometry of ranges, and return
Bio::RangeI compliant objects or triplets (start, stop, strand) from
which new ranges could be built.


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

=cut

sub intersection {
    my ($self, $other, $so) = @_;
    return unless $self->_testStrand($other, $so);

    $self->throw("start is undefined") unless defined $self->start;
    $self->throw("end is undefined") unless defined $self->end;
    $other->throw("Not a Bio::RangeI object") unless  $other->isa('Bio::RangeI');
    $other->throw("start is undefined") unless defined $other->start;
    $other->throw("end is undefined") unless defined $other->end;

    my @start = sort {$a<=>$b}
    ($self->start(), $other->start());
    my @end   = sort {$a<=>$b}
    ($self->end(),   $other->end());

    my $start = pop @start;
    my $end = shift @end;

    my $union_strand;  # Strand for the union range object.

    if($self->strand == $other->strand) {
	$union_strand = $other->strand;
    } else {
	$union_strand = 0;
    }

    if($start > $end) {
	return undef;
    } else {
      if( wantarray ) {
        return ($start, $end, $union_strand);
      }
      return $self->new('-start' => $start,
                        '-end' => $end,
                        '-strand' => $union_strand
                       );
    }
}

=head2 union

  Title   : union
  Usage   : my $union_range = $r1->union( @other_ranges ); (scalar context)
            OR
            my ( $start, $end, $strand ) = $r1->union( @other_ranges );
              (list context)
            OR
            my $union_range = Bio::RangeI->union( @ranges );
              (scalar context)
            OR
            my ( $start, $end, $strand ) = Bio::RangeI->union( @ranges );
              (list context)
  Function: finds the minimal range that contains all of the ranges
  Args    : a range or list of ranges to find the union of
  Returns : a new range object that contains all of the given ranges, or
            (in list context) the start, end, and strand of that range object.

=cut

sub union {
    my $self = shift;
    my @ranges = @_;
    if(ref $self) {
	unshift @ranges, $self;
    }

    my @start = sort {$a<=>$b}
    map( { $_->start() } @ranges);
    my @end   = sort {$a<=>$b}
    map( { $_->end()   } @ranges);

    my $start = shift @start;
    while( !defined $start ) {
	$start = shift @start;
    }

    my $end = pop @end;

    my $union_strand;  # Strand for the union range object.

    foreach(@ranges) {
	if(! defined $union_strand) {
	    $union_strand = $_->strand;
	    next;
	} else {
	    if($union_strand ne $_->strand) {
		$union_strand = 0;
		last;
	    }
	}
    }
    return undef unless $start or $end;
    if( wantarray() ) {
	return ( $start,$end,$union_strand);
    } else { 
	return $self->new('-start' => $start,
			  '-end' => $end,
			  '-strand' => $union_strand
			  );
    }
}

=head2 overlap_extent

 Title   : overlap_extent
 Usage   : ($a_unique,$common,$b_unique) = $a->overlap_extent($b)
 Function: Provides actual amount of overlap between two different
           ranges.
 Example :
 Returns : array of values for 
           - the amount unique to a
           - the amount common to both
           - the amount unique to b
 Args    : a range


=cut

sub overlap_extent{
    my ($a,$b) = @_;

    $a->throw("start is undefined") unless defined $a->start;
    $a->throw("end is undefined") unless defined $a->end;
    $b->throw("Not a Bio::RangeI object") unless  $b->isa('Bio::RangeI');
    $b->throw("start is undefined") unless defined $b->start;
    $b->throw("end is undefined") unless defined $b->end;
    
    my ($au,$bu,$is,$ie);
    if( ! $a->overlaps($b) ) {
	return ($a->length,0,$b->length);
    }

    if( $a->start < $b->start ) {
	$au = $b->start - $a->start;
    } else {
	$bu = $a->start - $b->start;
    }

    if( $a->end > $b->end ) {
	$au += $a->end - $b->end;
    } else {
	$bu += $b->end - $a->end;
    }
    my $intersect = $a->intersection($b);
    $ie = $intersect->end;
    $is = $intersect->start;

    return ($au,$ie-$is+1,$bu);
}

1;

__END__