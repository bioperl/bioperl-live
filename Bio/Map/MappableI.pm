# $Id$
#
# BioPerl module for Bio::Map::MappableI
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::MappableI - An object that can be placed in a map

=head1 SYNOPSIS

    # get a Bio::Map::MappableI somehow
    my $position = $element->map_position();
    # these methods will be important for building sorted lists
    if( $position->equals($p2) ) {
	# do something
    } elsif( $position->less_tha($p2) ) {} 
      elsif( $position->greater_than($p2) ) { }    


=head1 DESCRIPTION

This object handles the generic notion of an element placed on a
(linear) Map.  Elements can have multiple positions in a map such as
is the case of Restriction enzyme cut sites on a sequence map.  For
exact information about an element's position in a map one must query
the associate PositionI object which is accessible through the
position() method.

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

=head1 CONTRIBUTORS

Heikki Lehvaslaiho heikki@ebi.ac.uk
Lincoln Stein      lstein@cshl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# 'Let the code begin...


package Bio::Map::MappableI;
use vars qw(@ISA);
use strict;
use Bio::Root::RootI;
use Carp;

@ISA = qw(Bio::Root::RootI);

=head2 position

 Title   : position
 Usage   : my $position = $mappable->position(); 
 Function: Get/Set the Bio::Map::PositionI for a mappable element
 Returns : Bio::Map::PositionI
 Args    : (optional) Bio::Map::PositionI

=cut

sub position{
   my ($self,@args) = @_;
   $self->_abstractDeath('position');
}

=head2 equals

 Title   : equals
 Usage   : if( $mappable->equals($mapable2)) ...
 Function: Test if a position is equal to another position
 Returns : boolean
 Args    : Bio::Map::MappableI or Bio::Map::PositionI

=cut

sub equals{
   my ($self,@args) = @_;
   $self->_abstractDeath('equals');
}

=head2 less_than

 Title   : less_than
 Usage   : if( $mappable->less_than($m2) ) ...
 Function: Tests if a position is less than another position
 Returns : boolean
 Args    : Bio::Map::MappableI or Bio::Map::PositionI

=cut

sub less_than{
   my ($self,@args) = @_;
   $self->_abstractDeath('less_than');
}

=head2 greater_than

 Title   : greater_than
 Usage   : if( $mappable->greater_than($m2) ) ...
 Function: Tests if position is greater than another position
 Returns : boolean
 Args    : Bio::Map::MappableI or Bio::Map::PositionI

=cut

sub greater_than{
   my ($self,@args) = @_;
   $self->_abstractDeath('greater_than');
}

1;
