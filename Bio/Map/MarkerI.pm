# $Id$
#
# BioPerl module for Bio::Map::MarkerI
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::MarkerI - Interface for basic marker functionality

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the interface here

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
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Heikki Lehvaslaiho heikki@ebi.ac.uk
Lincoln Stein      lstein@cshl.org
Jason Stajich      jason@bioperl.org
Chad Matsalla      bioinformatics1@dieselwurks.com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Map::MarkerI;
use vars qw(@ISA);
use strict;
use Bio::Map::MappableI;

@ISA = qw(Bio::Map::MappableI);

=head2 name($new_name)

 Title   : name($new_name)
 Usage   : my $name = $o_usat->name($new_name) _or_
	   my $name = $o_usat->name()
 Function: Get/Set the name for this Microsatellite
 Returns : A scalar representing the current name of this Microsatellite
 Args    : If provided, the current name of this Microsatellite
	   will be set to $new_name.

=cut

sub name {
    my ($self) = @_;
    $self->_abstractDeath('name');
}

=head2 Bio::Map::MappableI methods

=head2 map_position

 Title   : position
 Usage   : my $position = $mappable->map_position(); 
 Function: Get/Set the Bio::Map::PositionI for a mappable element
 Returns : Bio::Map::PositionI
 Args    : (optional) Bio::Map::PositionI

=head2 equals

 Title   : equals
 Usage   : if( $mappable->equals($mapable2)) ...
 Function: Test if a position is equal to another position
 Returns : boolean
 Args    : Bio::Map::MappableI

=head2 less_than

 Title   : less_than
 Usage   : if( $mappable->less_than($m2) ) ...
 Function: Tests if a position is less than another position
 Returns : boolean
 Args    : Bio::Map::MappableI

=head2 greater_than

 Title   : greater_than
 Usage   : if( $mappable->greater_than($m2) ) ...
 Function: Tests if position is greater than another position
 Returns : boolean
 Args    : Bio::Map::MappableI

=cut

1;
