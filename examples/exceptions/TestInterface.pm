=head1 NAME

TestInterface - A simple subclass of Interface

=head1 DESCRIPTION

This module demonstrates how to use the generic Bio::Root::Interface
superclass.

=head1 AUTHOR

Steve Chervitz <sac@bioperl.org>

=cut

package TestInterface;

use Bio::Root::Interface;

@ISA = qw( Bio::Root::Interface );

sub data {
    my $self = shift;
    $self->throw_not_implemented;
}

sub foo {
    my $self = shift;
    $self->throw_not_implemented;
}

sub bar {
    my $self = shift;
    $self->throw_not_implemented;
}


1;
