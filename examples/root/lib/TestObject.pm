=head1 NAME

TestObject - An implementation of TestInterface

=head1 DESCRIPTION

This module attempts to provide an implementation of TestInterface and
is used for illustrating exception throwing using Graham Barr's
Error.pm object.

=head1 AUTHOR

Steve Chervitz E<lt>sac@bioperl.orgE<gt>

=cut

#'

package TestObject;

use strict;


# Define a special type of error "Bio::TestException" as a subclass of Error.
# Note two things:
#   1. The ISA declaration effectively defines our new Exception object.
#   2. This declaration doesn't have to be located in the Bio directory.
#   3. We don't have to use Bio::Root::Exception in this module.
#   4. If Error.pm isn't available this statement doesn't matter.
@Bio::TestException::ISA = qw( Bio::Root::Exception );

use base qw(Bio::Root::Root TestInterface);


# Note that we're not implementing foo(), so calling it
# will result in a Bio::Root::NotImplemented exception.

sub data {
    my ($self, $data) = @_;
    print "Setting test data ($data)\n" if $data && $self->verbose;
    $self->{'data'} = $data if $data;

   return $self->{'data'}
}

sub bar {

    my $self = shift;

    print "\nExecuting method bar() in TestObject\n" if $self->verbose;
    print "Throwing a Bio::TestException\n" if $self->verbose;

    my $message = "A Test error";

    # Bio::Root::Root::throw() will make use of Error.pm if present.
    # The type of Error is specified with a -class parameter.
    # If -class is not supplied, a Bio::Root::Exception is throw.
    # In this case, the argument can consist of a simple string.

    $self->throw( -class => 'Bio::TestException',
                  -text  => $message );

    print "Code within bar() below the throw that shouldn't be executed.\n" if $self->verbose;

}

1;
