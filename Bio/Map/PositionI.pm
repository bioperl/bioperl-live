# $Id$
#
# BioPerl module for Bio::Map::PositionI
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::PositionI - Abstracts the notion of a position having a value in the context of a marker and a Map

=head1 SYNOPSIS

    # do not use directly
    my $marker_position = $position->value;
    my $map             = $position->map;
    my $marker          = $position->marker;

=head1 DESCRIPTION

This object stores one of the postions that a mappable object
(e.g. Marker) may have in a map.

The method numeric() returns the position in a form that can be
compared between other positions of the same type.

A 'position', in addition to being a single point, can also be an area and so
can be imagined as a range and compared with other positions on the basis of
overlap, intersection etc.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Lincoln Stein, lstein-at-cshl.org
Heikki Lehvaslaiho, heikki-at-bioperl-dot-org
Sendu Bala, bix-at-sendu-dot-me-dot-uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Map::PositionI;
use vars qw(@ISA);
use strict;
use Bio::Root::RootI;
use Bio::RangeI;
use Carp;

@ISA = qw(Bio::Root::RootI Bio::RangeI);


=head2 map

 Title   : map
 Usage   : my $id = map->map;
 Function: Get/Set the map the position is in.
 Returns : L<Bio::Map::MapI>
 Args    : [optional] new L<Bio::Map::MapI>

=cut

sub map {
   my ($self, $value) = @_;
   $self->throw_not_implemented();
}

=head2 marker

 Title   : marker
 Usage   : my $id = marker->marker;
 Function: Get/Set the marker the position is in.
 Returns : L<Bio::Map::MarkerI>
 Args    : [optional] new L<Bio::Map::MarkerI>

=cut

sub marker {
   my ($self, $value) = @_;
   $self->throw_not_implemented();
}


=head2 value

 Title   : value
 Usage   : my $pos = $position->value;
 Function: Get/Set the value for this position
 Returns : scalar, value
 Args    : [optional] new value to set

=cut

sub value {
   my ($self, $value) = @_;
   $self->throw_not_implemented();
}

=head2 numeric

 Title   : numeric
 Usage   : my $num = $position->numeric;
 Function: Read-only method that is guarantied to return 
           representation for this position that can be compared with others
 Returns : numeric (int, real or range)
 Args    : none

=cut

sub numeric {
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 start

  Title   : start
  Usage   : $start = $range->start();
  Function: get/set the start of this range
  Returns : the start of this range
  Args    : optionally allows the start to be set
            using $range->start($start)

=cut

=head2 end

  Title   : end
  Usage   : $end = $range->end();
  Function: get/set the end of this range
  Returns : the end of this range
  Args    : optionally allows the end to be set
            using $range->end($end)

=cut

=head2 length

  Title   : length
  Usage   : $length = $range->length();
  Function: get/set the length of this range
  Returns : the length of this range
  Args    : optionally allows the length to be set
             using $range->length($length)

=cut

=head2 strand

  Title   : strand
  Usage   : $strand = $range->strand();
  Function: get the strand of this position; it is always 1 (the forward strand)
  Returns : 1, the strandedness
  Args    : none

=cut

sub strand {
   return 1;
}

=head1 Boolean Methods

These range-related methods return true or false. They throw an error if start and
end are not defined.

  $range->overlaps($otherRange) && print "Ranges overlap\n";

=head2 overlaps

  Title   : overlaps
  Usage   : if($r1->overlaps($r2)) { do stuff }
  Function: tests if $r2 overlaps $r1
  Args    : arg #1 = a range to compare this one to (mandatory)
            arg #2 = optional strand-testing arg ('strong', 'weak', 'ignore')
  Returns : true if the ranges overlap, false otherwise

=cut

=head2 contains

  Title   : contains
  Usage   : if($r1->contains($r2) { do stuff }
  Function: tests whether $r1 totally contains $r2 
  Args    : arg #1 = a range to compare this one to (mandatory)
	             alternatively, integer scalar to test
            arg #2 = optional strand-testing arg ('strong', 'weak', 'ignore')
  Returns : true if the argument is totally contained within this range

=cut

=head1 Geometrical methods

These methods do things to the geometry of ranges, and return
Bio::Map::PositionI compliant objects or triplets (start, stop, strand) from
which new positions could be built.

=head2 intersection

 Title   : intersection
 Usage   : ($start, $stop, $strand) = $r1->intersection($r2)
           ($start, $stop, $strand) = Bio::Map::PositionI->intersection(\@positions);
           my $containing_range = $r1->intersection($r2);
           my $containing_range = Bio::Map::PositionI->intersection(\@positions);
 Function: gives the range that is contained by all ranges
 Returns : undef if they do not overlap, or
           the range that they do overlap in an object like the calling one
           (with map transfered across from self or the first range in the array
           ref of ranges supplied, as long as the map is shared by all inputs),
           OR a three element array
 Args    : arg #1 = [REQUIRED] a position/range to compare this one to,
                    or an array ref of positions/ranges
           arg #2 = strand option ('strong', 'weak', 'ignore')

=cut

sub intersection {
    my $self = shift;
    
    @_ > 0 && $_[0] or return;
    
    unless (ref($_[0]) eq 'ARRAY') {
        if (wantarray()) {
            return $self->SUPER::intersection(@_);
        }
        my $output = $self->SUPER::intersection(@_);
        $output->map($self->map);
        return $output;
    }
    
    my @input_ranges = @{shift(@_)};
    @input_ranges > 1 or return;
    
    my ($output, $map);
    while (@input_ranges > 1) {
        unless ($output) {
            $output = shift(@input_ranges);
            $map = $output->map;
        }
        
        my $input_range = shift(@input_ranges);
        if ($map && (my $compare_map = $input_range->map)) {
            ($compare_map eq $map) or ($map = undef);
        }
        $output = $output->SUPER::intersection($input_range, @_);
    }
    
    $output->map($map);
    return $output->intersection(shift(@input_ranges), @_);
}

=head2 union

 Title   : union
 Usage   : ($start, $stop, $strand) = $r1->union($r2);
           ($start, $stop, $strand) = Bio::Map::PositionI->union(@positions);
           my $newpos = Bio::Map::PositionI->union(@positions);
 Function: finds the minimal position/range that contains all of the positions
 Returns : the position object containing all of the ranges in an object like
           the calling one (with map transfered across from self or the first
           range in the array ref of ranges supplied, as long as the map is
           shared by all inputs), OR a three element array
 Args    : a range or list of positions/ranges to find the union of

=cut
sub union {
    my $self = shift;
    
    if (wantarray()) {
        return $self->SUPER::union(@_);
    }
    
    my $output = $self->SUPER::union(@_);
    
    my $map;
    if (@_ > 1) {
        my @inputs = @_;
        $map = shift(@inputs)->map;
        $map or return $output;
        
        foreach my $input (@inputs) {
            if (my $compare_map = $input->map) {
                $compare_map ne $map and return $output;
            }
        }
        
        $output->map($map);
    }
    else {
        $output->map($self->map);
    }
    
    return $output;
}

=head2 overlap_extent

 Title   : overlap_extent
 Usage   : ($a_unique,$common,$b_unique) = $a->overlap_extent($b)
 Function: Provides actual amount of overlap between two different
           positions
 Example :
 Returns : array of values containing the length unique to the calling 
           position, the length common to both, and the length unique to 
           the argument position
 Args    : a position

=cut

#*** should this be overridden from RangeI?

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

#*** this should be overridden from RangeI

1;
