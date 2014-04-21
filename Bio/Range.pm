#
# BioPerl module for Bio::Range
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
#
# Copywright Matthew Pocock
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
#

=head1 NAME

Bio::Range - Pure perl RangeI implementation

=head1 SYNOPSIS

  $range = Bio::Range->new(-start=>10, -end=>30, -strand=>+1);
  $r2 = Bio::Range->new(-start=>15, -end=>200, -strand=>+1);

  print join(', ', $range->union($r2)), "\n";
  print join(', ', $range->intersection($r2)), "\n";

  print $range->overlaps($r2), "\n";
  print $range->contains($r2), "\n";

=head1 DESCRIPTION

This provides a pure perl implementation of the BioPerl range
interface.

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
the bugs and their resolution.  Bug reports can be submitted via  the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Heikki Lehvaslaiho

Email heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Range;

use strict;
use Carp;
use integer;


use base qw(Bio::Root::Root Bio::RangeI);

=head1 Constructors

=head2 new

  Title   : new
  Usage   : $range = Bio::Range->new(-start => 100, -end=> 200, -strand = +1);
  Function: generates a new Bio::Range
  Returns : a new range
  Args    : -strand (defaults to 0) and any two of (-start, -end, -length),
            the third will be calculated

=cut

sub new {
  my ($caller, @args) = @_;
  my $self = $caller->SUPER::new(@args);
  my ($strand, $start, $end, $length) = 
      $self->_rearrange([qw(STRAND 
			    START
			    END 
			    LENGTH
			    )],@args);
  $self->strand($strand || 0);

  if(defined $start ) {
      $self->start($start);
      if(defined $end) {
	  $self->end($end);
      } elsif(defined $length) {
	  $self->end($self->start()+ $length - 1);
      }
  } elsif(defined $end && defined $length ) {
      $self->end($end);
      $self->start($self->end() - $length + 1);
  }
  return $self;
}

=head2 unions

 Title   : unions
 Usage   : @unions = Bio::Range->unions(@ranges);
 Function: generate a list of non-intersecting Bio::Range objects
           from a list of Bio::Range objects which may intersect
 Returns : a list of Bio::Range objects
 Args    : a list of Bio::Range objects


=cut

sub unions {
  my ($class,@i) = @_;

  my $i = 0;
  my %i = map { $i++ => $_ } @i;

  my $lastsize = scalar(keys %i);

  do {

    foreach my $j (sort { $i{$a}->start <=> $i{$b}->start } keys %i){
      foreach my $k (sort { $i{$a}->start <=> $i{$b}->start } keys %i){

        #it may have been replaced by a union under the key of
        #the overlapping range, we are altering the hash in-place
        next unless $i{$j};

        next if $i{$k}->end   < $i{$j}->start;
        last if $i{$k}->start > $i{$j}->end;

        if($i{$j}->overlaps($i{$k})){
          my($start,$end,$strand) = $i{$j}->union($i{$k});
          delete($i{$k});
          $i{$j} = Bio::Range->new( -start => $start , -end => $end , -strand => $strand );
        }
      }
    }

    goto DONE if scalar(keys %i) == $lastsize;
    $lastsize = scalar(keys %i);

    #warn $lastsize;

  } while(1);

  DONE:

  return values %i;
}


=head1 Member variable access

These methods let you get at and set the member variables

=head2 start

  Title    : start
  Function : return or set the start co-ordinate
  Example  : $s = $range->start(); $range->start(7);
  Returns  : the value of the start co-ordinate
  Args     : optionally, the new start co-ordinate
  Overrides: Bio::RangeI::start

=cut

sub start {
    my ($self,$value) = @_;
    if( defined $value) {
	$self->throw("'$value' is not an integer.\n") 
	    unless $value =~ /^[-+]?\d+$/;
        $self->{'start'} = $value;
    }
    return $self->{'start'};
}

=head2 end

  Title    : end
  Function : return or set the end co-ordinate
  Example  : $e = $range->end(); $range->end(2000);
  Returns  : the value of the end co-ordinate
  Args     : optionally, the new end co-ordinate
  Overrides: Bio::RangeI::end

=cut

sub end {

    my ($self,$value) = @_;
    if( defined $value) {
	$self->throw("'$value' is not an integer.\n") 
	    unless $value =~ /^[-+]?\d+$/;
        $self->{'end'} = $value;
    }
    return $self->{'end'};
}

=head2 strand

  Title    : strand
  Function : return or set the strandedness
  Example  : $st = $range->strand(); $range->strand(-1);
  Returns  : the value of the strandedness (-1, 0 or 1)
  Args     : optionally, the new strand - (-1, 0, 1) or (-, ., +).
  Overrides: Bio::RangeI::strand

=cut

{

my %VALID_STRAND = (
    -1      => -1,
    0       => 0,
    1       => 1,
    '+'     => 1,
    '-'     => -1,
    '.'     => 0
);

sub strand {
  my $self = shift;
  if(@_) {
    my $val = shift;
    if (exists $VALID_STRAND{$val}) {
        $self->{'strand'} = $VALID_STRAND{$val};
    } else {
        $self->throw("Invalid strand: $val");
    }
  }
  return $self->{'strand'};
}

}

=head2 length

  Title    : length
  Function : returns the length of this range
  Example  : $length = $range->length();
  Returns  : the length of this range, equal to end - start + 1
  Args     : if you attempt to set the length an exception will be thrown
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
  Function : tests whether $r1 totally contains $r2
  Args     : a range to test for being contained
  Returns  : true if the argument is totally contained within this range
  Inherited: Bio::RangeI

=head2 equals

  Title    : equals
  Usage    : if($r1->equals($r2))
  Function : test whether $r1 has the same start, end, length as $r2
  Args     : a range to test for equality
  Returns  : true if they are describing the same range
  Inherited: Bio::RangeI

=head1 Geometrical methods

These methods do things to the geometry of ranges, and return
triplets (start, end, strand) from which new ranges could be built.

=head2 intersection

  Title    : intersection
  Usage    : ($start, $stop, $strand) = $r1->intersection($r2)
  Function : gives the range that is contained by both ranges
  Args     : a range to compare this one to
  Returns  : nothing if they do not overlap, or the range that they do overlap
  Inherited: Bio::RangeI::intersection

=cut

=head2 union

  Title    : union
  Usage    : ($start, $stop, $strand) = $r1->union($r2);
           : ($start, $stop, $strand) = Bio::Range->union(@ranges);
  Function : finds the minimal range that contains all of the ranges
  Args     : a range or list of ranges
  Returns  : the range containing all of the ranges
  Inherited: Bio::RangeI::union

=cut

1;
