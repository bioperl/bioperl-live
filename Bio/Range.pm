#
# BioPerl module for Bio::Range
#
# Cared for by Matthew Pocock <mrp@sanger.ac.uk>
#
# Copywright Matthew Pocock
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
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

Bio::Range - Pure perl RangeI implementation

=head1 DESCRIPTION

This provides a pure perl implementation of the  BioPerl range interface.

Ranges are modeled as having (start, end, length, strand). They
use Bio-coordinates - all points >= start and <= end are within
the range. End is always greater-than or equal-to start, and
length is greather than or equal to 1. The behaviour of a range
is undefined if ranges with negative numbers or zero are used.

So, in summary:

  length = end - start + 1
  end >= start
  strand = (-1 | 0 | +1)

=head1 SYNOPSIS

  $range = new Bio::Range(-start=>10, -end=>30, -strand=>+1);
  $r2 = new Bio::Range(-start=>15, -end=>200, -strand=>+1);

  print join(', ', $range->union($r2), "\n";
  print join(', ', $range->intersection($r2), "\n";
  print $range->overlaps($r2), "\n";
  print $range->contains($r2), "\n";

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

The rest of the documentation details each of the object methods. Internal metho
ds are usually preceded with a _

=cut

package Bio::Range;

use strict;
use Carp;

use Bio::RangeI;

use vars qw(@ISA);

@ISA = qw(Bio::RangeI);

=head1 Constructors

=head2 new

  Title   : new
  Usage   : $range = Bio::Range->new(-start => 100, -end=> 200, -strand = +1);
  Function: generates a new Bio::Range
  Returns : a new range
  Args    : two of (-start, -end, '-length') - the third is calculated
          : -strand (defaults to 0)

=cut

sub new {
  my $thingy = shift;
  my $package = ref($thingy) || $thingy;
  my $self = bless {}, $package;
  my $usageMessage = "Specify exactly two of -start, -end, '-length'";
  my %args = @_;
  $self->strand($args{-strand} || 0);
  
  if(exists($args{'-start'}) && exists($args{'-end'}) &&
     exists($args{'-length'})) {
      confess $usageMessage;
  }
  
  if(exists($args{'-start'})) {
      $self->start($args{'-start'});
      if(exists($args{'-end'})) {
	  $self->end($args{'-end'});
      } elsif($args{'-length'}) {
	  $self->end($self->start()+$args{'-length'}-1);
      } else {
	  confess $usageMessage;
      }
  } elsif(exists($args{'-end'}) && exists($args{'-length'})) {
      $self->end($args{'-end'});
      $self->start($self->end() - $args{'-length'} + 1);
  } else {
      confess $usageMessage;
  }
  return $self;
}

=head1 Member variable access

These methods let you get at and set the member variables

=head2 start

  Title    : start
  Function : return or set the start co-ordinate
  Example  : $s = $range->start(); $range->start(7);
  Returns  : the value of the start co-ordinate
  Args     : optionaly, the new start co-ordinate
  Overrides: Bio::RangeI::start

=cut

sub start {
  my $self = shift;
  @_ ? $self->{start} = shift
     : $self->{start};
}

=head2 end

  Title    : end
  Function : return or set the end co-ordinate
  Example  : $e = $range->end(); $range->end(2000);
  Returns  : the value of the end co-ordinate
  Args     : optionaly, the new end co-ordinate
  Overrides: Bio::RangeI::end

=cut

sub end {
  my $self = shift;
  @_ ? $self->{end} = shift
     : $self->{end};
}

=head2 strand

  Title    : strand
  Function : return or set the strandidness
  Example  : $st = $range->strand(); $range->strand(-1);
  Returns  : the value of the strandedness (-1, 0 or 1)
  Args     : optionaly, the new strand - (-1, 0, 1) or (-, ., +).
  Overrides: Bio::RangeI::Strand

=cut

sub strand {
  my $self = shift;
  if(@_) {
    my $val = shift;
    $val =~ tr/+/1/;
    $val =~ tr/-/-1/;
    $val =~ tr/./0/;
    if($val == -1 || $val == 0 || $val == 1 ) {
      $self->{strand} = $val;
    }
  }  
  return $self->{strand};
}

=head2 length

  Title    : length
  Function : returns the length of this range
  Example  : $length = $range->length();
  Returns  : the length of this range, equal to end - start + 1
  Args     : if you attempt to set the length, and exeption will be thrown
  Overrides: Bio::RangeI::Length

=cut

sub length {
  my $self = shift;
  if(@_) {
    confess ref($self), "->length() is read-only";
  }
  return $self->end() - $self->start() + 1;
}

=head2 toString

  Title   : toString
  Function: stringifies this range
  Example : print $range->toString(), "\n";
  Returns : a string representation of this range

=cut

sub toString {
  my $self = shift;
  return  "(${\$self->start}, ${\$self->end}) strand=${\$self->strand}";
}

=head1 Boolean Methods

These methods return true or false.

 $range->overlaps($otherRange) && print "Ranges overlap\n";

=head2 overlaps

  Title    : overlaps
  Usage    : if($r1->overlaps($r2)) { do stuff }
  Function : tests if $r2 overlaps $r1
  Args     : a range to test for overlap with
  Returns  : true if the ranges overlap, false otherwise
  Inherited: Bio::RangeI

=head2 contains

  Title    : contains
  Usage    : if($r1->contains($r2) { do stuff }
  Function : tests wether $r1 totaly contains $r2
  Args     : a range to test for being contained
  Returns  : true if the argument is totaly contained within this range
  Inherited: Bio::RangeI

=head2

  Title    : equals
  Usage    : if($r1->equals($r2))
  Function : test whether $r1 has the same start, end, length as $r2
  Args     : a range to test for equality
  Returns  : true if they are describing the same range
  Inherited: Bio::RangeI

=head1 Geometrical methods

These methods do things to the geometry of ranges, and return
triplets (start, end, strand) from which new ranges could be built.

=head2

  Title    : intersection
  Usage    : ($start, $stop, $strand) = $r1->intersection($r2)
  Function : gives the range that is contained by both ranges
  Args     : a range to compare this one to
  Returns  : nothing if they don't overlap, or the range that they do overlap
  Inherited: Bio::RangeI

=head2

  Title    : union
  Usage    : ($start, $stop, $strand) = $r1->union($r2);
           : ($start, $stop, $strand) = Bio::RangeI->union(@ranges);
  Function : finds the minimal range that contains all of the ranges
  Args     : a range or list of ranges to find the union of
  Returns  : the range containing all of the ranges
  Inherited: Bio::RangeI

1;
