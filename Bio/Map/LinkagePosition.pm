# BioPerl module for Bio::Map::LinkagePosition
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::LinkagePosition - Create a Position for a Marker that will be placed
	                        on a Bio::Map::LinkageMap

=head1 SYNOPSIS

    use Bio::Map::Position;
    my $position = Bio::Map::LinkagePosition->new(-positions => 1,
						 -distance => 22.1 );

	# can get listing of positions
    my @positions = $position->each_position;


=head1 DESCRIPTION

Position for a Bio::Map::MarkerI compliant object that will be
placed on a Bio::Map::LinkageMap. See L<Bio::Map::MarkerI> and
L<Bio::Map::LinkageMap> for details

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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

=head1 AUTHOR - Chad Matsalla

Email bioinformatics1@dieselwurks.com

=head1 CONTRIBUTORS

Lincoln Stein, lstein@cshl.org
Heikki Lehvaslaiho, heikki-at-bioperl-dot-org
Jason Stajich jason@bioperl.org
Sendu Bala bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::LinkagePosition;
use strict;


use base qw(Bio::Map::OrderedPosition);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Map::LinkagePosition->new(-positions => $position,
				                                   -distance => $distance);
 Function: Builds a new Bio::Map::LinkagePosition object
 Returns : Bio::Map::LinkagePosition
 Args    : -order => the relative order of this marker on a linkage map
 	       -positions => positions on a map

=cut

=head2 Bio::Map::PositionI methods

=cut

=head2 order

 Title   : order
 Usage   : $o_position->order($order)
           my $order = $o_position->order()
 Function: get/set the order position of this position in a map
 Returns : int
 Args    : none to get, int to set

=cut


1;
