#
# BioPerl module for Bio::Map::CytoMarker
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::CytoMarker - An object representing a marker.

=head1 SYNOPSIS

  $o_usat = Bio::Map::CytoMarker->new(-name=>'Chad Super Marker 2',
				 -position => $pos);

=head1 DESCRIPTION

This object handles markers with a positon in a cytogenetic map known.
This marker will have a name and a position.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Heikki Lehvaslaiho 

Email heikki-at-bioperl-dot-org

=head1 CONTRIBUTORS

Chad Matsalla      bioinformatics1@dieselwurks.com
Lincoln Stein      lstein@cshl.org
Jason Stajich      jason@bioperl.org
Sendu Bala         bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::CytoMarker;
use strict;
use Bio::Map::CytoPosition;

use base qw(Bio::Map::Marker);


=head2 Bio::Map::MarkerI methods

=cut

=head2 get_position_object

 Title   : get_position_class
 Usage   : my $position = $marker->get_position_object();
 Function: To get an object of the default Position class
           for this Marker. Subclasses should redefine this method.
           The Position returned needs to be a L<Bio::Map::PositionI> with
		   -element set to self.
 Returns : L<Bio::Map::PositionI>
 Args    : none for an 'empty' PositionI object, optionally
           Bio::Map::MapI and value string to set the Position's -map and -value
           attributes.

=cut

sub get_position_object {
   my ($self, $map, $value) = @_;
   $map ||= $self->default_map;
   if ($value) {
	  $self->throw("Value better be scalar, not [$value]") unless ref($value) eq '';
   }
   
   my $pos = Bio::Map::CytoPosition->new();
   $pos->map($map) if $map;
   $pos->value($value) if $value;
   $pos->element($self);
   return $pos;
}


=head2 Comparison methods

The numeric values for cutogeneic loctions go from the p tip of
chromosome 1, down to the q tip and similarly throgh consecutive
chromosomes, through X and end the the q tip of X. See
L<Bio::Map::CytoPosition::cytorange> for more details.

=cut

=head2 New methods

=cut

=head2 get_chr

 Title   : get_chr
 Usage   : my $mychr = $marker->get_chr();
 Function: Read only method for the  chromosome string of the location.
           A shortcut to $marker->position->chr().
 Returns : chromosome value
 Args    : [optional] new chromosome value

=cut

sub get_chr {
    my ($self) = @_;
    return unless $self->position;
    return $self->position->chr;
}

1;

