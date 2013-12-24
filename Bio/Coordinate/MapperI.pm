package Bio::Coordinate::MapperI;
use utf8;
use strict;
use warnings;
use parent qw(Bio::Root::RootI);

# ABSTRACT: Interface describing coordinate mappers.
# AUTHOR:   Heikki Lehvaslaiho <heikki@bioperl.org>
# OWNER:    Heikki Lehvaslaiho
# LICENSE:  Perl_5

=head1 SYNOPSIS

  # not to be used directly

=head1 DESCRIPTION

MapperI defines methods for classes capable for mapping locations
between coordinate systems.

=cut

=head2 in

 Title   : in
 Usage   : $obj->in('peptide');
 Function: Set and read the input coordinate system.
 Example :
 Returns : value of input system
 Args    : new value (optional), Bio::LocationI

=cut

sub in {
   my ($self,$value) = @_;

   $self->throw_not_implemented();

}

=head2 out

 Title   : out
 Usage   : $obj->out('peptide');
 Function: Set and read the output coordinate system.
 Example :
 Returns : value of output system
 Args    : new value (optional), Bio::LocationI

=cut

sub out {
   my ($self,$value) = @_;

   $self->throw_not_implemented();
}

=head2 swap

 Title   : swap
 Usage   : $obj->swap;
 Function: Swap the direction of mapping: input <-> output)
 Example :
 Returns : 1
 Args    :

=cut

sub swap {
   my ($self) = @_;

   $self->throw_not_implemented();

}

=head2 test

 Title   : test
 Usage   : $obj->test;
 Function: test that both components are of same length
 Example :
 Returns : ( 1 | undef )
 Args    :

=cut

sub test {
   my ($self) = @_;

   $self->throw_not_implemented();
}

=head2 map

 Title   : map
 Usage   : $newpos = $obj->map($loc);
 Function: Map the location from the input coordinate system
           to a new value in the output coordinate system.
 Example :
 Returns : new value in the output coordiante system
 Args    : Bio::LocationI

=cut

sub map {
   my ($self,$value) = @_;

   $self->throw_not_implemented();

}

=head2 return_match

 Title   : return_match
 Usage   : $obj->return_match(1);
 Function: A flag to turn on the simplified mode of
           returning only one joined Match object or undef
 Example :
 Returns : boolean
 Args    : boolean (optional)

=cut

sub return_match {
   my ($self,$value) = @_;
   if( defined $value) {
       $value ? ( $self->{'_return_match'} = 1 ) :
                ( $self->{'_return_match'} = 0 );
   }
   return $self->{'_return_match'} || 0 ;
}

1;
