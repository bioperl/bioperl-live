# MapI.pm,v 1.5 2002/04/18 12:52:48 jason Exp
#
# BioPerl module for Bio::Map::MapI
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::MapI - Interface for describing Map objects in bioperl 

=head1 SYNOPSIS

    # get a MapI somehowe
    my $name   = $map->name();     # string
    my $length = $map->length();   # integer
    my $species= $map->species;    # Bio::Species
    my $type   = $map->type();     # genetic/sts/rh/

=head1 DESCRIPTION

This object describes the basic functionality of a Map in bioperl.
Maps are anything from Genetic Map to Sequence Map to and Assembly Map
to Restriction Enzyme to FPC.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org            - General discussion
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

Lincoln Stein, lstein@cshl.org
Heikki Lehvaslaiho, heikki@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::MapI;
use vars qw(@ISA);
use strict;
use Bio::Root::RootI;
use Carp;

@ISA = qw(Bio::Root::RootI);

=head2 species

 Title   : species
 Usage   : my $species = $map->species;
 Function: Get/Set Species for a map
 Returns : L<Bio::Species> object
 Args    : (optional) Bio::Species

=cut

sub species{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 units

 Title   : units
 Usage   : $map->units('cM');
 Function: Get/Set units for a map
 Returns : units for a map
 Args    : units for a map (string)

=cut

sub units{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 type

 Title   : type
 Usage   : my $type = $map->type
 Function: Get/Set Map type
 Returns : String coding map type
 Args    : (optional) string

=cut

sub type {
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 name

 Title   : name
 Usage   : my $name = $map->name
 Function: Get/Set Map name
 Returns : Map name
 Args    : (optional) string

=cut

sub name {
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 length

 Title   : length
 Usage   : my $length = $map->length();
 Function: Retrieves the length of the map, 
           It is possible for the length to be unknown
           for maps such as Restriction Enzyme, will return undef
           in that case
 Returns : integer representing length of map in current units
           will return undef if length is not calculateable
 Args    : none

=cut

sub length{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 unique_id

 Title   : unique_id
 Usage   : my $id = $map->unique_id;
 Function: Get/Set the unique ID for this map
 Returns : a unique identifier
 Args    : [optional] new identifier to set 

=cut

sub unique_id{
   my ($self,$id) = @_;
   $self->throw_not_implemented();
}

=head2 add_element

 Title   : add_element
 Usage   : $map->add_element($marker)
 Function: Add a Bio::Map::MappableI object to the Map
 Returns : none
 Args    : Bio::Map::MappableI object

=cut

sub add_element{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 each_element

 Title   : each_element
 Usage   : my @elements = $map->each_element;
 Function: Retrieves all the elements in a map
           unordered
 Returns : Array of Map elements (L<Bio::Map::MarkerI>)
 Args    :


=cut

sub each_element{
   my ($self) = @_;
   $self->throw_not_implemented();
}

1;
