# $Id$
#
# BioPerl module for Bio::Map::CytoMarker
#
# Cared for by Heikki Lehvaslaiho heikki@ebi.ac.uk
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::CytoMarker - An object representing a marker.

=head1 SYNOPSIS

  $o_usat = new Bio::Map::CytoMarker(-name=>'Chad Super Marker 2',
				 -position => $pos);

=head1 DESCRIPTION

This object handles markers with a positon in a cytogenetic map known.
This marker will have a name and a position.


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Heikki Lehvaslaiho 

Email heikki@ebi.ac.uk

=head1 CONTRIBUTORS

Chad Matsalla      bioinformatics1@dieselwurks.com
Lincoln Stein      lstein@cshl.org
Jason Stajich      jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::CytoMarker;
use vars qw(@ISA);
use strict;
use Bio::Map::Marker;
use Bio::Map::CytoPosition;
use Bio::RangeI;

@ISA = qw(Bio::Map::Marker Bio::RangeI  );


=head2 Bio::Map::MarkerI methods

=cut

=head2 get_position_object

 Title   : get_position_class
 Usage   : my $pos = $marker->get_position_object();
 Function: To get an object of the default Position class
           for this Marker. Subclasses should redefine this method.
           The Position needs to be L<Bio::Map::PositionI>.
 Returns : L<Bio::Map::CytoPosition>
 Args    : none

=cut

sub get_position_object {
   my ($self) = @_;
   return new Bio::Map::CytoPosition();
}


=head2 Comparison methods

The numeric values for cutogeneic loctions go from the p tip of
chromosome 1, down to the q tip and similarly throgh consecutive
chromosomes, through X and end the the q tip of X. See
L<Bio::Map::CytoPosition::cytorange> for more details.

The numeric values for cytogenetic positions are ranges of type
L<Bio::Range>, so MarkerI type of operators (equals, less_than,
greater_than) are not very meaningful, but they might be of some use
combined with L<Bio::RangeI> methods (overlaps, contains, equals,
intersection, union). equals(), present in both interfaces is treated
as a more precice RangeI method.

CytoMarker has a method L<get_chr> which might turn out to be useful
in this context.

The L<less_than> and L<greater_than> methods are implemented by
comparing the end values of the range, so you better first check that
markers do not overlap, or you get an opposite result than expected.
The numerical values are not metric, so avarages are not meaningful.

Note: These methods always return a value. A false value (0) might
mean that you have not set the position! Check those warnings.

=cut

=head2 Bio::Map::MarkerI comparison methods

=cut

=head2 tuple

 Title   : tuple
 Usage   : ($me, $you) = $self->_tuple($compare)
 Function: Utility method to extract numbers and test for missing values.
 Returns : two ranges or tuple of -1
 Args    : Bio::Map::MappableI or Bio::Map::PositionI

=cut

=head2 less_than

 Title   : less_than
 Usage   : if( $mappable->less_than($m2) ) ...
 Function: Tests if a position is less than another position
 Returns : boolean
 Args    : Bio::Map::MappableI  or Bio::Map::PositionI

=cut


sub less_than {
    my ($self,$compare) = @_;

    my ($me, $you) = $self->tuple($compare);
    return 0 if $me == -1 or $you == -1 ;

    $me  = $me->end;
    $you  = $you->start;

    print STDERR "me=$me, you=$you\n" if $self->verbose > 0;
    return $me < $you;
}

=head2 greater_than

 Title   : greater_than
 Usage   : if( $mappable->greater_than($m2) ) ...
 Function: Tests if position is greater than another position
 Returns : boolean
 Args    : Bio::Map::MappableI or Bio::Map::PositionI

=cut


sub greater_than {
    my ($self,$compare) = @_;

    my ($me, $you) = $self->tuple($compare);
    return 0 if $me == -1 or $you == -1 ;

    $me  = $me->start;
    $you  = $you->end;
    print STDERR "me=$me, you=$you\n" if $self->verbose > 0;
    return $me > $you;
}

=head2 RangeI methods

=cut


=head2 equals

 Title   : equals
 Usage   : if( $mappable->equals($mapable2)) ...
 Function: Test if a position is equal to another position
 Returns : boolean
 Args    : Bio::Map::MappableI or Bio::Map::PositionI

=cut

sub equals {
    my ($self,$compare) = @_;

    my ($me, $you) = $self->tuple($compare);
    return 0 unless $me->isa('Bio::RangeI') and $you->isa('Bio::RangeI');

    return $me->equals($you);
}

=head2 overlaps

  Title    : overlaps
  Usage    : if($r1->overlaps($r2)) { do stuff }
  Function : tests if $r2 overlaps $r1
  Args     : a range to test for overlap with
  Returns  : true if the ranges overlap, false otherwise
  Inherited: Bio::RangeI

=cut

sub overlaps {
    my ($self,$compare) = @_;

    my ($me, $you) = $self->tuple($compare);
    return 0 unless $me->isa('Bio::RangeI') and $you->isa('Bio::RangeI');

    return $me->overlaps($you);
}

=head2 contains

  Title    : contains
  Usage    : if($r1->contains($r2) { do stuff }
  Function : tests wether $r1 totaly contains $r2
  Args     : a range to test for being contained
  Returns  : true if the argument is totaly contained within this range
  Inherited: Bio::RangeI

=cut

sub contains {
    my ($self,$compare) = @_;

    my ($me, $you) = $self->tuple($compare);
    return 0 unless $me->isa('Bio::RangeI') and $you->isa('Bio::RangeI');
    print STDERR "me=", $me->start. "-", $me->end, " ",
    "you=", $you->start. "-", $you->end, "\n"
	if $self->verbose > 0;

    return $me->contains($you);
}

=head2 intersection

  Title    : intersection
  Usage    : ($start, $stop, $strand) = $r1->intersection($r2)
  Function : gives the range that is contained by both ranges
  Args     : a range to compare this one to
  Returns  : nothing if they do not overlap, or the range that they do overlap
  Inherited: Bio::RangeI::intersection

=cut

sub intersection {
    my ($self,$compare) = @_;

    my ($me, $you) = $self->tuple($compare);
    return 0 unless $me->isa('Bio::RangeI') and $you->isa('Bio::RangeI');

    return $me->intersection($you);
}

=head2 union

  Title    : union
  Usage    : ($start, $stop, $strand) = $r1->union($r2);
           : ($start, $stop, $strand) = Bio::Range->union(@ranges);
  Function : finds the minimal range that contains all of the ranges
  Args     : a range or list of ranges to find the union of
  Returns  : the range containing all of the ranges
  Inherited: Bio::RangeI::union

=cut

sub union {
    my ($self,$compare) = @_;

    my ($me, $you) = $self->tuple($compare);
    return 0 unless $me->isa('Bio::RangeI') and $you->isa('Bio::RangeI');

    return $me->union($you);
}


=head2 New methods

=cut


=head2 get_chr

 Title   : get_chr
 Usage   : my $mychr = $marker->get_chr();
 Function: Read only method for the  chromosome string of the location.
           A shotrcut to $marker->position->chr().
 Returns : chromosome value
 Args    : [optional] new chromosome value

=cut


sub get_chr {
    my ($self) = @_;
    return undef unless $self->position;
    return $self->position->chr;
}

1;

