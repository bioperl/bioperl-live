# PositionI.pm,v 1.7 2002/04/18 12:52:48 jason Exp
#
# BioPerl module for Bio::Map::PositionI
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::PositionI - Abstracts the notion of a position  having
  a value in the context of a marker and a Map

=head1 SYNOPSIS

    # do not use directly

=head1 DESCRIPTION

This object stores one of the postions a that a mappable object
(e.g. Marker) may have in a map (e.g. restriction enzymes or a SNP
mapped to several chromosomes).

The method numeric() returns the position in a form that can be
compared between other positions of the same type.

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

Lincoln Stein, lstein@cshl.org
Heikki Lehvaslaiho, heikki@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Map::PositionI;
use vars qw(@ISA);
use strict;
use Bio::Root::RootI;
use Carp;

@ISA = qw(Bio::Root::RootI);


=head2 map

 Title   : map
 Usage   : my $id = map->$map;
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
 Usage   : my $id = marker->$marker;
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

1;
