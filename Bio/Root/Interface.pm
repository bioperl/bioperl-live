#-----------------------------------------------------------------
# $Id$
#
# BioPerl module Bio::Root::Interface
#
# Cared for by Steve Chervitz <steve_chervitz@affymetrix.com>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

=head1 NAME

Bio::Root::Interface - Superclass for an abstract interface 

=head1 SYNOPSIS

    package MyInterface;

    require Bio::Root::Interface;

    @MyInterface::ISA = qw( Bio::Root::Interface );

    sub foo {
    	my $self = shift;
    	$self->throw_not_implemented;
    }


=head1 DESCRIPTION

B<Bio::Root::Interface> defines routine methods that is intended for
use by a module defining an abstract interface. Such a module
specifies the public methods that any implementation must provide.
The interface itself cannot be instantiated (i.e, a C<new()> method is
not provided).

B<Bio::Root::Interface> provides a single method,
L</throw_not_implemented()> that is intended to be called by all
methods within an interface so that if an implementation fails to
override them, an exception will be thrown. 

For example, say the C<MyInterface> module defines a method called C<foo()>.
Since this method is considered abstract and should be implemented by any
module implementing C<MyInterface>, the C<MyInterface::foo()> method 
would consist of the following:

    sub foo {
    	my $self = shift;
    	$self->throw_not_implemented;
    }

So, if an implementer of C<MyInterface> forgets to implement C<foo()>
and a user of the implementation calls C<foo()>, a
B<Bio::Exception::NotImplemented> exception will result.

Unfortunately, failure to implement a method can only be determined at
run time (i.e., you can't verify that an implementation is complete by
running C<perl -wc> on it). So it should be standard practice for a test
of an implementation to check each method and verify that it doesn't
throw a B<Bio::Exception::NotImplemented>.

=head2 Exception Throwing

When throwing an exception, this module will use Graham Barr's
B<Error.pm> module if it is installed. If not, it will attempt to call
the object method C<throw()> on itself, which would be available if
the implementing object also inherits from B<Bio::Root::RootI>. If the
object method C<throw()> is not available, C<Carp::confess()> is used.

There's a way to suppress the use of B<Error.pm> even if it is
installed.  Do this by defining the C<$DONT_USE_ERROR> within a BEGIN
{} block at the top of your main script before loading any other
module. This is intended for testing purposes.  Under normal
circumstances, Error.pm will be used if it is installed and will not
be used if it is not installed.

=head1 SEE ALSO

See the C<examples/exceptions> directory of the Bioperl distribution for 
working demo code.

B<Error.pm> available from CPAN

B<Bio::Root::Exception>

=head1 AUTHOR - Steve Chervitz 

steve_chervitz@affymetrix.com

=head1 COPYRIGHT

Copyright (c) 2001 Steve Chervitz. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 METHODS

=cut

# Let the code begin...'

package Bio::Root::Interface;

use strict;
use Carp;

=head2 B<throw_not_implemented()>

 Purpose : Throws a Bio::Root::NotImplemented exception.
           Intended for use in the method definitions of 
           abstract interface modules where methods are defined
           but are intended to be overridden by subclasses.
 Usage   : $object->throw_not_implemented();
 Example : sub method_foo { 
             $self = shift; 
             $self->throw_not_implemented();
           }
 Returns : n/a
 Args    : n/a
 Throws  : A Bio::Root::NotImplemented exception (via Error.pm)
           The message of the exception contains
             - the name of the method 
             - the name of the interface 
             - the name of the implementing class 
           The name of the method is also placed within the 'value'
           field of the Error object.

           This method will work even when Error.pm isn't installed,
  	   in which case $self->throw will be used.
           If the object doesn't have a throw() method, 
           Carp::confess() will be used.


=cut   #'

sub throw_not_implemented {
    my $self = shift;
    my $package = ref $self;
    my $iface = caller(0);
    my @call = caller(1);
    my $meth = $call[3];

    my $message = "Abstract method \"$meth\" is not implemented by package $package.\n" .
		   "This is not your fault - author of $package should be blamed!\n";

    # Checking if Error.pm is available in case the object isn't decended from
    # Bio::Root::Root, which knows how to check for Error.pm.

    if( $Bio::Root::Root::ERRORLOADED ) {
	throw Bio::Root::NotImplemented ( -text   => $message,
                                          -value  => $meth
                                        );
    }
    elsif( $self->can('throw') ) {
	 $self->throw( $message );
    }
    else {
	confess $message ;
    }
}


=head2 B<warn_not_implemented()>

 Purpose : Generates a warning that a method has not been implemented.
           Intended for use in the method definitions of 
           abstract interface modules where methods are defined
           but are intended to be overridden by subclasses.
           Generally, throw_not_implemented() should be used,
           but warn_not_implemented() may be used if the method isn't
           considered essential and convenient no-op behavior can be 
           provided within the interface.
 Usage   : $object->warn_not_implemented( method-name-string );
 Example : $self->warn_not_implemented( "get_foobar" );
 Returns : Calls $self->warn on this object, if available.
           If the object doesn't have a warn() method,
           Carp::carp() will be used.
 Args    : n/a


=cut   #'

sub warn_not_implemented {
    my $self = shift;
    my $package = ref $self;
    my $iface = caller(0);
    my @call = caller(1);
    my $meth = $call[3];

    my $message = "Abstract method \"$meth\" is not implemented by package $package.\n" .
		   "This is not your fault - author of $package should be blamed!\n";

    if( $self->can('warn') ) {
        $self->warn( $message );
    }
    else {
	carp $message ;
    }
}


1;




