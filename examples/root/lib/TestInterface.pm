=head1 NAME

TestInterface - A simple subclass of Interface

=head1 DESCRIPTION

This module demonstrates how to use the generic Bio::Root::RootI
superclass.

=head1 AUTHOR

Steve Chervitz E<lt>sac@bioperl.orgE<gt>

=cut

package TestInterface;


use base qw(Bio::Root::RootI);

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
