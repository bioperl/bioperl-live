# $Id$
#
# BioPerl module for Bio::RangeI
#
# Cared for by Matthew Pocock <mrp@sanger.ac.uk>
#
# Copywright Matthew Pocock
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::RangeI - Range interface

=head1 SYNOPSIS

None.

=head1 DESCRIPTION

This provides a standard BioPerl range interface that should be
implemented by any object that wants to be treated as a range. This
serves purely as an abstract base class for implementers and can
not be instantiated.

Ranges are modeled as having (start, end, length, strand). They
use Bio-coordinates - all points >= start and <= end are within
the range. End is always greater-than or equal-to start, and
length is greather than or equal to 1. The behaviour of a range
is undefined if ranges with negative numbers or zero are used.

So, in summary:

  length = end - start + 1
  end >= start
  strand = (-1 | 0 | +1)

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Matthew Pocock

Email mrp@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal metho ds are usually preceded with a _

=cut

package Bio::RangeI;

use strict;
use Carp;

use vars '%STRAND_OPTIONS';

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
# Prints out a method like:
# Abstract method stop defined in interface Bio::RangeI not implemented by package You::BadRange
sub _abstractDeath {
  my $self = shift;
  my $package = ref $self || $self;
  my $caller = (caller)[1];
  
  confess "Abstract method '$caller' defined in interface Bio::RangeI not implemented by pacakge $package";
}

# returns true if strands are equal and non-zero
sub _strong {
  my ($r1, $r2) = @_;
  my ($s1, $s2) = ($r1->strand(), $r2->strand());
  
  return $s1 != 0 and $s2 != 0 and
         $s1 == $s2;
}

# returns true if strands are equal or either is zero
sub _weak {
  my ($r1, $r2) = @_;
  my ($s1, $s2) = ($r1->strand(), $r2->strand());

  return $s1 == 0 or $s2 == 0 or
         $s1 == $s2;
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

=head2

  Title   : new
  Function: confesses if you try to instantiate a RangeI
          : RangeI is an interface, so RangeI->new should never be called
	  : To make a range, instantiate one of the implementing classes. e.g.
	  : $range = Bio::Range->new(-start=>20, -end=>2000, -strand=>1)

=cut

sub new {
  shift->_abstractDeath();
}

=head2

  Title   : start
  Usage   : $start = $range->start();
  Function: get/set the start of this range
  Returns : the start of this range
  Args    : optionaly allows the start to be set
          : using $range->start($start)

=cut

sub start {
  shift->_abstractDeath();
}

=head2

  Title   : end
  Usage   : $end = $range->end();
  Function: get/set the end of this range
  Returns : the end of this range
  Args    : optionaly allows the end to be set
          : using $range->end($start)

=cut

sub end {
  shift->_abstractDeath();
}

=head2

  Title   : length
  Usage   : $length = $range->length();
  Function: get/set the length of this range
  Returns : the length of this range
  Args    : optionaly allows the length to be set
          : using $range->length($length)

=cut

sub length {
  shift->_abstractDeath();
}

=head2

  Title   : strand
  Usage   : $strand = $range->strand();
  Function: get/set the strand of this range
  Returns : the strandidness (-1, 0, +1)
  Args    : optionaly allows the strand to be set
          : using $range->strand($strand)

=cut

sub strand {
  shift->_abstractDeath();
}

=head1 Boolean Methods

These methods return true or false.

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
  return
    ($self->_testStrand($other, $so) and 
    not (
      ($self->start() > $other->end() or
       $self->end() < $other->start()   )
        ));
}

=head2 contains

  Title   : contains
  Usage   : if($r1->contains($r2) { do stuff }
  Function: tests whether $r1 totally contains $r2
  Args    : arg #1 = a range to compare this one to (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') (optional)
  Returns : true if the argument is totaly contained within this range

=cut

sub contains {
  my ($self, $other, $so) = @_;
  if(ref $other) { # a range object?
    return ($self->_testStrand($other, $so)      and
      $other->start() >= $self->start() and
      $other->end() <= $self->end());
  } else { # a scalar?
    return ($other >= $self->start() and $other <= $self->end());
  }
}

=head2

  Title   : equals
  Usage   : if($r1->equals($r2))
  Function: test whether $r1 has the same start, end, length as $r2
  Args    : a range to test for equality
  Returns : true if they are describing the same range

=cut

sub equals {
  my ($self, $other, $so) = @_;
  return ($self->_testStrand($other, $so)   and
          $self->start() == $other->start() and
          $self->end()   == $other->end()       );
}

=head1 Geometrical methods

These methods do things to the geometry of ranges, and return
triplets (start, stop, strand) from which new ranges could be built.

=head2

  Title   : intersection
  Usage   : ($start, $stop, $strand) = $r1->intersection($r2)
  Function: gives the range that is contained by both ranges
  Args    : arg #1 = a range to compare this one to (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') (optional)
  Returns : nothing if they do not overlap, or the range that they do overlap

=cut

sub intersection {
  my ($self, $other, $so) = @_;
  return unless $self->_testStrand($other, $so);

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
    return;     
  } else {
    return ($start, $end, $union_strand);
  }
}

=head2

  Title   : union
  Usage   : ($start, $stop, $strand) = $r1->union($r2);
          : ($start, $stop, $strand) = Bio::RangeI->union(@ranges);
  Function: finds the minimal range that contains all of the ranges
  Args    : a range or list of ranges to find the union of
  Returns : the range containing all of the ranges

=cut

sub union {
  my $thingy = shift;
  my @ranges = @_;
  if(ref $thingy) {
    unshift @ranges, $thingy;
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
  
  return ($start, $end, $union_strand);
}

=head2 overlap_extent

 Title   : overlap_extent
 Usage   : ($a_unique,$common,$b_unique) = $a->overlap_extent($b)
 Function: Provides actual amount of overlap between two different
           ranges, being the amount unique to a, the amount unique to
           b and the amount common to both.
 Example :
 Returns : 
 Args    :


=cut

sub overlap_extent{
   my ($a,$b) = @_;

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

   ($is,$ie) = $a->intersection($b);

  # print STDERR "Returning $au,$ie,$is,$bu\n";

   return ($au,$ie-$is+1,$bu);
}

1;
