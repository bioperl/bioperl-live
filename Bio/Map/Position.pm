# $Id$
#
# BioPerl module for Bio::Map::Position
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::Position - A single position of a Marker in a Map

=head1 SYNOPSIS

    use Bio::Map::Position;
    my $position = new Bio::Map::Position(-map => $map, 
					  -marker => $marker
					  -value => 100
					  );

=head1 DESCRIPTION

This object is an implementation of the PositionI interface that
handles the specific values of a position.  This allows an element
(e.g. Marker) to have multiple positions within a map and still be
treated as a single entity.  

This does not directly handle the concept of a relative map in which
no known exact positions exist but markers are just ordered relative
to one another - in that case arbitrary values must be assigned for
position values.

No units are assumed here - units are handled by context of which Map
a position is placed in or the subclass of this Position.

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
Chad Matsalla, bioinformatics1@dieselwurks.com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Map::Position;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::Map::PositionI;

@ISA = qw(Bio::Root::Root Bio::Map::PositionI );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Map::Position();
 Function: Builds a new Bio::Map::Position object 
 Returns : Bio::Map::Position
 Args    : -map     a <Bio::Map::MapI> object
           -marker  a <Bio::Map::MarkerI> object
           -value   string or number

=cut

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
  
    my ($map, $marker, $value) = 
	$self->_rearrange([qw( MAP 
			       MARKER 
			       VALUE
			       )], @args);

    $map     && $self->map($map);
    $marker  && $self->marker($marker);
    $value   && $self->value($value);

    return $self;
}

=head2 map

 Title   : map
 Usage   : my $id = map->$map;
 Function: Get/Set the map the position is in.
 Returns : L<Bio::Map::MapI>
 Args    : [optional] new L<Bio::Map::MapI>

=cut

sub map {
   my ($self,$map) = @_;
   if( defined $map ) {
       $self->throw("This is [$map], not a Bio::Map::MapI object")
	   unless $map->isa('Bio::Map::MapI');
       $self->{'_map'} = $map;
   }
   return $self->{'_map'};
}

=head2 marker

 Title   : marker
 Usage   : my $id = marker->$marker;
 Function: Get/Set the marker the position is in.
 Returns : L<Bio::Map::MarkerI>
 Args    : [optional] new L<Bio::Map::MarkerI>

=cut

sub marker {
   my ($self,$marker) = @_;
   if( defined $marker ) {
       $self->thow("This is [$marker], not a Bio::Map::MarkerI object")
	   unless $marker->isa('Bio::Map::MarkerI');
       $self->{'_marker'} = $marker;
   }
   return $self->{'_marker'};
}

=head2 value

 Title   : value
 Usage   : my $pos = $position->value;
 Function: Get/Set the value for this postion
 Returns : scalar, value
 Args    : [optional] new value to set

=cut

sub value {
   my ($self,$value) = @_;
   if( defined $value ) {
       $self->{'_value'} = $value;
   }
   return $self->{'_value'};
}

=head2 numeric

 Title   : numeric
 Usage   : my $num = $position->numeric;
 Function: Read-only method that is guarantied to return a numeric 
           representation for this position. 
 Returns : numeric (int or real) 
 Args    : none

=cut

sub numeric {
   my ($self) = @_;
   my $num = $self->{'_value'} || 0;

   # expand this to cover scientific notation, too!
   $self->throw("This value [$num] is not numeric!")
       unless $num && $num =~ /^[+-]?[\d.]+$/;
   return $num;
}

1;
