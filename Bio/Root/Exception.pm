#-----------------------------------------------------------------
# $Id$
#
# BioPerl module Bio::Root::Exception
#
# Cared for by Steve Chervitz <steve_chervitz@affymetrix.com>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

=head1 NAME

Bio::Root::Exception - Generic exception objects for Bioperl

=head1 SYNOPSIS

* Exceptions defined in this module, all inherit from Bio::Root::Exception:

=over 8

=item B<Bio::Root::Exception>

=item B<Bio::Root::NotImplemented>

=item B<Bio::Root::IOException>

=item B<Bio::Root::FileOpenException>

=item B<Bio::Root::SystemException>

=item B<Bio::Root::BadParameter>

=item B<Bio::Root::OutOfRange>

=item B<Bio::Root::NoSuchThing>

=back

* Throwing and catching a Bio::Root::FileOpenException:

    use Bio::Root::Exception;
    use Error qw(:try);

    $file = shift;

    try {
        open (IN, $file) || throw Bio::Root::FileOpenException (
	  -text => "Can't open file $file for reading",
          -value => $!,
	  -object => $self );
    } 

    catch Bio::Root::FileOpenException with {
        my $err = shift;
        print STDERR "Using default input file: $default_file\n";
        open (IN, $default_file) || die "Can't open $default_file";
    }

    otherwise {
        my $err = shift;
    	print STDERR "An unexpected exception occurred: \n$err";
    
	# By placing an the error object reference within double quotes,
	# you're invoking its stringify() method.

    };   # IMPORTANT: Don't forget this final semicolon!


* Throwing Bio::Root::Exceptions via B<Bio::Root::RootI::throw()>:

     $obj->throw(-class => 'Bio::Root::BadParameter',
                 -text  => "The value you supplied looks bad: $value",
                 -value => $value );


* Defining a new Exception type as a subclass of Bio::Root::Exception:

    @Bio::TestException::ISA = qw( Bio::Root::Exception );


=head1 DESCRIPTION

B<Bio::Root::Exception> defines some generic exceptions intended for
use within Bioperl modules and scripts. These Exceptions all derive
from Bio::Root::Exception which can be used as a base class for all
Bioperl exceptions. This would allow a user to write a handler for all
Bioperl-derived exceptions as follows:

           use Bio::Whatever;
           use Error qw(:try);

           try {
                # some code that depends on Bioperl
           }
           catch Bio::Root::Exception with {
               my $err = shift;
               print "A Bioperl exception occurred:\n$err\n";
           };

Bio::Root::Exceptions can also be thrown indirectly via
B<Bio::Root::RootI::throw()> which automatically uses Error.pm 
when it's available. 

To take advantage of this capability, you must specify the type 
of exception you want to throw by using named parameters 
in your throw call. The C<-class> parameter is a string containing
the name of the Bio::Root::Exception subclass you want to throw.
The C<-text> parameter should contain an informative error message.
You can optionally define a C<-value> parameter to contain the value that
caused the exception. This can be processed when the exception is caught.

Note that Bio::Root::Exception does not need to be imported into
your module (or script) namespace in order to throw exceptions
via Bio::Root::RootI::throw().

The exceptions in Bio::Root::Exception are extensions of Graham Barr's
B<Error.pm> module available from CPAN.  Despite this dependency, the
Bio::Root::Exception module does not explicitly C<require Error>.
This permits Bio::Root::Exception to be loaded even when
Error.pm is not available.

=head1 SEE ALSO

See the C<examples/exceptions> directory of the Bioperl distribution for 
working demo code.

B<Error.pm> (available from CPAN)

Error.pm is helping to guide the design of exception handling in Perl 6. 
See these RFC's: 

     http://dev.perl.org/rfc/63.pod 

     http://dev.perl.org/rfc/88.pod


=head1 AUTHOR - Steve Chervitz 

steve_chervitz@affymetrix.com

=head1 COPYRIGHT

Copyright (c) 2001 Steve A. Chervitz. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 EXCEPTIONS

=cut

# Define some generic exceptions.'

package Bio::Root::Exception;

use strict;

my $debug = $Error::Debug;  # Prevents the "used only once" warning.

=head2 B<Bio::Root::Exception>

 Purpose : A generic base class for all BioPerl exceptions.
           By including a "catch Bio::Root::Exception" block, you
           should be able to trap all BioPerl exceptions.
 Example : throw Bio::Root::Exception( 
               -text   => "A generic exception",
               -object => $self );

=cut

#---------------------------------------------------------
@Bio::Root::Exception::ISA = qw( Error );
#---------------------------------------------------------

=head2 Methods defined by Bio::Root::Exception

=over 4

=item B< pretty_format() >

 Purpose : Get a nicely formatted string containing information about the 
           exception. Format is similar to that produced by 
           Bio::Root::RootI::throw(), with the addition of the name of
           the exception class in the EXCEPTION line and some other
           data available via the Error object.
 Example : print $error->pretty_format;

=cut

sub pretty_format {
    my $self = shift;
    my $msg = $self->text;
    my $stack = '';
    if( $Error::Debug ) {
      $stack = $self->_reformat_stacktrace();
    }
    my $value_string = $self->value ? "VALUE: ".$self->value."\n" : "";
    my $class = ref($self);

    my $title = "------------- EXCEPTION: $class -------------";
    my $footer = "\n" . '-' x CORE::length($title);
    my $out = $title . "\n" .
       "MSG: ".$msg."\n". $value_string. $stack. $footer . "\n";
    return $out;
}


# Modifications:
#   1. Shift the file:line data in line i to line i+1.
#   2. change xxx::__ANON__() to "try{} block"
#   3. skip the "require" and "Error::subs::try" stack entries (boring)
# This means that the first line in the stack won't have any file:line data
# But this isn't a big issue since it's for a Bio::Root::-based method 
# that doesn't vary from exception to exception.
sub _reformat_stacktrace {
    my $self = shift;
    my $msg = $self->text;
    my $stack = $self->stacktrace();
    $stack =~ s/$msg//;
    my @stack = split( /\n/, $stack);
    my @new_stack = ();
    my ($method, $file, $linenum, $prev_file, $prev_linenum);
    my $stack_count = 0;
    foreach my $i( 0..$#stack ) {
     # print "STACK-ORIG: $stack[$i]\n";
      if( ($stack[$i] =~ /^\s*([^(]+)\s*\(.*\) called at (\S+) line (\d+)/) ||
	   ($stack[$i] =~ /^\s*(require 0) called at (\S+) line (\d+)/)) {
	($method, $file, $linenum) = ($1, $2, $3);
	$stack_count++;
      }
      else{
	next;
      }
      if( $stack_count == 1 ) {
	push @new_stack, "STACK: $method";
       ($prev_file, $prev_linenum) = ($file, $linenum);
	next;
      }
      
      if( $method =~ /__ANON__/ ) {
	$method = "try{} block";
      }
      if( ($method =~ /^require/ and $file =~ /Error\.pm/ ) ||
	  ($method =~ /^Error::subs::try/ ) )   {
	last;
      }
      push @new_stack, "STACK: $method $prev_file:$prev_linenum";
      ($prev_file, $prev_linenum) = ($file, $linenum);
    }
    return join "\n", @new_stack;
}

=item B< stringify() >

 Purpose : Overrides Error::stringify() to call pretty_format(). 
           This is called automatically when an exception object 
           is placed between double quotes.
 Example : catch Bio::Root::Exception with {
              my $error = shift;
              print "$error";
           }

See Also: L<pretty_format()>

=cut

sub stringify {
    my ($self, @args) = @_;
    return $self->pretty_format( @args );
}

=back

=head1 Subclasses of Bio::Root::Exception 


=head2 B<Bio::Root::NotImplemented>

 Purpose : Indicates that a method has not been implemented.
 Example : throw Bio::Root::NotImplemented( 
               -text   => "Method \"foo\" not implemented in module FooBar.",
               -value  => "foo" );

=cut

#---------------------------------------------------------
@Bio::Root::NotImplemented::ISA = qw( Bio::Root::Exception );
#---------------------------------------------------------

=head2 B<Bio::Root::IOException>

 Purpose : Indicates that some input/output-related trouble has occurred.
 Example : throw Bio::Root::IOException( 
               -text   => "Can't save data to file $file.",
	       -value  => $! );

=cut

#---------------------------------------------------------
@Bio::Root::IOException::ISA = qw( Bio::Root::Exception );
#---------------------------------------------------------


=head2 B<Bio::Root::FileOpenException>

 Purpose : Indicates that a file could not be opened.
 Example : throw Bio::Root::FileOpenException( 
               -text   => "Can't open file $file for reading.",
	       -value  => $! );

=cut

#---------------------------------------------------------
@Bio::Root::FileOpenException::ISA = qw( Bio::Root::IOException );
#---------------------------------------------------------


=head2 B<Bio::Root::SystemException>

 Purpose : Indicates that a system call failed.
 Example : unlink($file) or throw Bio::Root::SystemException( 
               -text   => "Can't unlink file $file.",
	       -value  => $! );

=cut

#---------------------------------------------------------
@Bio::Root::SystemException::ISA = qw( Bio::Root::Exception );
#---------------------------------------------------------


=head2 B<Bio::Root::BadParameter>

 Purpose : Indicates that one or more parameters supplied to a method 
           are invalid, unspecified, or conflicting.
 Example : throw Bio::Root::BadParameter( 
               -text   => "Required parameter \"-foo\" was not specified",
               -value  => "-foo" );

=cut

#---------------------------------------------------------
@Bio::Root::BadParameter::ISA = qw( Bio::Root::Exception );
#---------------------------------------------------------


=head2 B<Bio::Root::OutOfRange>

 Purpose : Indicates that a specified (start,end) range or 
           an index to an array is outside the permitted range.
 Example : throw Bio::Root::OutOfRange( 
               -text   => "Start coordinate ($start) cannot be less than zero.",
               -value  => $start  );

=cut

#---------------------------------------------------------
@Bio::Root::OutOfRange::ISA = qw( Bio::Root::Exception );
#---------------------------------------------------------


=head2 B<Bio::Root::NoSuchThing>

 Purpose : Indicates that a requested thing cannot be located 
           and therefore could possibly be bogus.
 Example : throw Bio::Root::NoSuchThing( 
               -text   => "Accession M000001 could not be found.",
               -value  => "M000001"  );

=cut

#---------------------------------------------------------
@Bio::Root::NoSuchThing::ISA = qw( Bio::Root::Exception );
#---------------------------------------------------------


1;

