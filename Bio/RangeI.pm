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
greater than or equal to 1. The behaviour of a range is undefined if
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

=head2 start

  Title   : start
  Usage   : $start = $range->start();
  Function: get/set the start of this range
  Returns : the start of this range
  Args    : optionally allows the start to be set
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
  Args    : optionally allows the end to be set
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
  Args    : optionally allows the length to be set
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
  Args    : optionally allows the strand to be set
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
  Returns : true if the argument is totally contained within this range

=cut

sub contains {
  my ($self, $other, $so) = @_;
  $self->throw("start is undefined") unless defined $self->start;
  $self->throw("end is undefined") unless defined $self->end;

  if(defined $other && ref $other) { # a range object?
      $other->throw("Not a Bio::RangeI object: $other") unless  $other->isa('Bio::RangeI');
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
  Usage   : if($r1->equals($r2))
  Function: test whether $r1 has the same start, end, length as $r2
  Args    : a range to test for equality
  Returns : true if they are describing the same range

=cut

sub equals {
    my ($self, $other, $so) = @_;

    $self->throw("start is undefined") unless defined $self->start;
    $self->throw("end is undefined") unless defined $self->end;
    $other->throw("Not a Bio::RangeI object") unless  $other->isa('Bio::RangeI');
    $other->throw("start is undefined") unless defined $other->start;
    $other->throw("end is undefined") unless defined $other->end;

    return ($self->_testStrand($other, $so)   and
	    $self->start() == $other->start() and
	    $self->end()   == $other->end()       );
}

=head1 Geometrical methods

These methods do things to the geometry of ranges, and return
Bio::RangeI compliant objects or triplets (start, stop, strand) from
which new ranges could be built.


=head2 intersection

  Title   : intersection
  Usage   : ($start, $stop, $strand) = $r1->intersection($r2)
            my $containing_range = $r1->intersection($r2)
  Function: gives the range that is contained by both ranges
  Args    : arg #1 = a range to compare this one to (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') (optional)
  Returns : undef if they do not overlap, 
            or a range if they do overlap (in the form of an object 
            like the calling one, OR a three element array)

=cut

sub intersection {
    my ($self, $other, $so) = @_;
    $self->throw("missing arg: you need to pass in another feature")
      unless $other;
    return unless $self->_testStrand($other, $so);

    $self->throw("start is undefined") unless defined $self->start;
    $self->throw("end is undefined") unless defined $self->end;
    $self->throw("Not a Bio::RangeI object: $other") unless ref($other);
    $other->throw("Not a Bio::RangeI object: $other") unless $other->isa('Bio::RangeI');
    $other->throw("start is undefined") unless defined $other->start;
    $other->throw("end is undefined") unless defined $other->end;

    my @start = sort {$a<=>$b}
    ($self->start(), $other->start());
    my @end   = sort {$a<=>$b}
    ($self->end(),   $other->end());

    my $start = pop @start;
    my $end = shift @end;

    my $union_strand;  # Strand for the union range object.

    if(defined($self->strand) &&
       defined($other->strand) &&
       $self->strand == $other->strand) {
	$union_strand = $other->strand;
    } else {
	$union_strand = 0;
    }

    if($start > $end) {
	return undef;
    } else {
	if( wantarray() ) {
	    return ($start, $end, $union_strand);
	}
	else {
	    return $self->new('-start' => $start,
			      '-end' => $end,
			      '-strand' => $union_strand
			     );
	}

    }
}

=head2 union

   Title   : union
    Usage   : ($start, $stop, $strand) = $r1->union($r2);
            : ($start, $stop, $strand) = Bio::Range->union(@ranges);
              my $newrange = Bio::Range->union(@ranges);
    Function: finds the minimal range that contains all of the ranges
    Args    : a range or list of ranges
    Returns : the range object containing all of the ranges
              (in the form of an object like the calling one, OR 
              a three element array)

=cut

sub union {
    my $self = shift;
    my @ranges = @_;
    if ($self eq "Bio::RangeI") {
	$self = "Bio::Range";
	$self->warn("calling static methods of an interface is deprecated; use $self instead");
    }
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
	    if(not defined $_->strand or $union_strand ne $_->strand) {
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

    if( ! $a->overlaps($b) ) {
	return ($a->length,0,$b->length);
    }

    my ($au,$bu) = (0, 0);
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
    my $ie = $intersect->end;
    my $is = $intersect->start;

    return ($au,$ie-$is+1,$bu);
}

=head2 disconnected_ranges

    Title   : disconnected_ranges
    Usage   : my @disc_ranges = Bio::Range->disconnected_ranges(@ranges);
    Function: finds the minimal set of ranges such that each input range
              is fully contained by at least one output range, and none of
              the output ranges overlap
    Args    : a list of ranges
    Returns : a list of objects of the same type as the input 
              (conforms to RangeI)

=cut

sub disconnected_ranges {
    my $self = shift;
    if ($self eq "Bio::RangeI") {
	$self = "Bio::Range";
	$self->warn("calling static methods of an interface is deprecated; use $self instead");
    }
    my @inranges = @_;
    if(ref $self) {
	unshift @inranges, $self;
    }

    my @outranges = (); # disconnected ranges

    # iterate through all input ranges $inrange,
    # adding each input range to the set of output ranges @outranges,
    # provided $inrange does not overlap ANY range in @outranges
    # - if it does overlap an outrange, then merge it
    foreach my $inrange (@inranges) {
	my $intersects = 0;
	my @outranges_new = ();
	my @intersecting_ranges = ();

        # iterate through all @outranges, testing if it intersects
        # current $inrange; if it does, merge and add to list
        # of @intersecting_ranges, otherwise add $outrange to
        # the new list of outranges that do NOT intersect
	for (my $i=0; $i<@outranges; $i++) {
	    my $outrange = $outranges[$i];
	    my $intersection = $inrange->intersection($outrange);
	    if ($intersection) {
		$intersects = 1;
		my $union = $inrange->union($outrange);
		push(@intersecting_ranges, $union);
	    }
	    else {
		push(@outranges_new, $outrange);
	    }
	}
	@outranges = @outranges_new;
        # @outranges now contains a list of non-overlapping ranges
        # that do not intersect the current $inrange

	if (@intersecting_ranges) {
	    if (@intersecting_ranges > 1) {
		# this sf intersected > 1 range, which means that
		# all the ranges it intersects should be joined
		# together in a new range
                my $merged_range =
                  $self->union(@intersecting_ranges);
		push(@outranges, $merged_range);
		     
	    }
	    else {
		# exactly 1 intersecting range
		push(@outranges, @intersecting_ranges);
	    }
	}
	else {
	    # no intersections found - new range
	    push(@outranges, 
		 $self->new('-start'=>$inrange->start,
			    '-end'=>$inrange->end,
			    '-strand'=>$inrange->strand,
			   ));
	}
    }
    return @outranges;
}


1;
