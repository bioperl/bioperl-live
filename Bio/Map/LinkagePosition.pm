# BioPerl module for Bio::Map::LinkagePosition
#
# Cared for by Chad Matsalla <bioinformatics1@dieselwurks.com>
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
    my $position = new Bio::Map::LinkagePosition(-positions => 1,
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

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Chad Matsalla

Email bioinformatics1@dieselwurks.com

=head1 CONTRIBUTORS

Lincoln Stein, lstein@cshl.org
Heikki Lehvaslaiho, heikki@ebi.ac.uk
Jason Stajich jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Map::LinkagePosition;
use vars qw(@ISA);
use strict;
require 'dumpvar.pl';

use Bio::Map::OrderedPosition;

@ISA = qw(Bio::Map::OrderedPosition);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Map::LinkagePosition(-positions => $position,
				-distance => $distance );
 Function: Builds a new Bio::Map::LinkagePosition object
 Returns : Bio::Map::LinkagePosition
 Args    : -order => the relative order of this marker on a linkage map
 	   -positions => positions on a map
=cut

=head2 Bio::Map::PositionI methods

=cut

=head2 order

 Title   : order
 Usage   : $o_position->order($new_position) _or_
           $o_position->order()
 Function: get/set the order position of this position in a map
 Returns :
 Args    : If $new_position is provided, the current position of this Position
           will be set to $new_position.

=cut


1;
