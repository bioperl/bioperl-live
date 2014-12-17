package Bio::Root::Exception;
use strict;

# ABSTRACT: generic exception objects for Bioperl
# AUTHOR:   Steve Chervitz <sac@bioperl.org>
# OWNER:    2001 Steve Chervitz
# LICENSE:  Perl_5

=head1 SYNOPSIS

=head2 Throwing exceptions using L<Error.pm throw|Error::throw>:

    use Bio::Root::Exception;
    use Error;

    # Set Error::Debug to include stack trace data in the error messages
    $Error::Debug = 1;

    $file = shift;
    open my $IN, '<', $file
        or Bio::Root::FileOpenException->throw("Could not read file '$file': $!");

=head2 Throwing exceptions using L<Bioperl throw|Bio::Root::Root/throw>:

    # Here we have an object that ISA Bio::Root::Root, so it inherits throw().

    open my $IN, '<', $file
        or $object->throw(-class => 'Bio::Root::FileOpenException',
                          -text  => "Could not read file '$file'",
                          -value => $!);

=head2 Catching and handling exceptions using L<Error.pm try|Error/try>:

    use Bio::Root::Exception;
    use Error qw(:try);

    # Note that we need to import the 'try' tag from Error.pm

    # Set Error::Debug to include stack trace data in the error messages
    $Error::Debug = 1;

    my $file = shift;
    my $IN;
    try {
        open $IN, '<', $file
            or Bio::Root::FileOpenException->throw("Could not read file '$file': $!");
    }
    catch Bio::Root::FileOpenException with {
        my $err = shift;
        print STDERR "Using default input file: $default_file\n";
        open $IN, '<', $default_file or die "Could not read file '$default_file': $!";
    }
    otherwise {
        my $err = shift;
        print STDERR "An unexpected exception occurred: \n$err";

        # By placing an the error object reference within double quotes,
        # you're invoking its stringify() method.
    }
   finally {
       # Any code that you want to execute regardless of whether or not
       # an exception occurred.
   };
   # the ending semicolon is essential!


=head2 Defining a new Exception type as a subclass of Bio::Root::Exception:

    @Bio::TestException::ISA = qw( Bio::Root::Exception );

=head1 DESCRIPTION

=head2 Exceptions defined in L<Bio::Root::Exception>

These are generic exceptions for typical problem situations that could arise
in any module or script.

=for :list
* C<Bio::Root::Exception()>
* C<Bio::Root::NotImplemented()>
* C<Bio::Root::IOException()>
* C<Bio::Root::FileOpenException()>
* C<Bio::Root::SystemException()>
* C<Bio::Root::BadParameter()>
* C<Bio::Root::OutOfRange()>
* C<Bio::Root::NoSuchThing()>

Using defined exception classes like these is a good idea because it
indicates the basic nature of what went wrong in a convenient,
computable way.

If there is a type of exception that you want to throw
that is not covered by the classes listed above, it is easy to define
a new one that fits your needs. Just write a line like the following
in your module or script where you want to use it (or put it somewhere
that is accessible to your code):

    @NoCanDoException::ISA = qw( Bio::Root::Exception );

All of the exceptions defined in this module inherit from a common
base class exception, Bio::Root::Exception. This allows a user to
write a handler for all Bioperl-derived exceptions as follows:

           use Bio::Whatever;
           use Error qw(:try);

           try {
                # some code that depends on Bioperl
           }
           catch Bio::Root::Exception with {
               my $err = shift;
               print "A Bioperl exception occurred:\n$err\n";
           };

So if you do create your own exceptions, just be sure they inherit
from Bio::Root::Exception directly, or indirectly by inheriting from a
Bio::Root::Exception subclass.

The exceptions in Bio::Root::Exception are extensions of Graham Barr's
L<Error> module available from CPAN.  Despite this dependency, the
L<Bio::Root::Exception> module does not explicitly C<require Error>.
This permits Bio::Root::Exception to be loaded even when
Error.pm is not available.

=head2 Throwing exceptions within Bioperl modules

Error.pm is not part of the Bioperl distibution, and may not be
present within  any given perl installation. So, when you want to
throw an exception in a Bioperl module, the safe way to throw it
is to use L<Bio::Root::Root/throw> which can use Error.pm
when it's available. See documentation in Bio::Root::Root for details.

=head1 SEE ALSO

See the C<examples/exceptions> directory of the Bioperl distribution for
working demo code.

L<Bio::Root::Root/throw> for information about throwing
L<Bio::Root::Exception>-based exceptions.

L<Error> (available from CPAN, author: GBARR)

Error.pm is helping to guide the design of exception handling in Perl 6.
See these RFC's:

     http://dev.perl.org/rfc/63.pod

     http://dev.perl.org/rfc/88.pod

=head1 EXCEPTIONS

=cut

my $debug = $Error::Debug;  # Prevents the "used only once" warning.
my $DEFAULT_VALUE = "__DUMMY__";  # Permits eval{} based handlers to work

=head2 L<Bio::Root::Exception>

 Purpose : A generic base class for all BioPerl exceptions.
           By including a "catch Bio::Root::Exception" block, you
           should be able to trap all BioPerl exceptions.
 Example : throw Bio::Root::Exception("A generic exception", $!);

=cut

#---------------------------------------------------------
@Bio::Root::Exception::ISA = qw( Error );
#---------------------------------------------------------

=head1 Methods defined by Bio::Root::Exception

=head2 new

 Purpose : Guarantees that -value is set properly before
           calling Error::new().

 Arguments: key-value style arguments same as for Error::new()

     You can also specify plain arguments as ($message, $value)
     where $value is optional.

     -value, if defined, must be non-zero and not an empty string
     in order for eval{}-based exception handlers to work.
     These require that if($@) evaluates to true, which will not
     be the case if the Error has no value (Error overloads
     numeric operations to the Error::value() method).

     It is OK to create Bio::Root::Exception objects without
     specifying -value. In this case, an invisible dummy value is used.

     If you happen to specify a -value of zero (0), it will
     be replaced by the string "The number zero (0)".

     If you happen to specify a -value of empty string (""), it will
     be replaced by the string "An empty string ("")".

=cut

sub new {
    my ($class, @args) = @_;
    my ($value, %params);
    if( @args % 2 == 0 && $args[0] =~ /^-/) {
        %params = @args;
        $value = $params{'-value'};
    }
    else {
        $params{-text} = $args[0];
        $value = $args[1];
    }

    if( defined $value ) {
        $value = "The number zero (0)" if $value =~ /^\d+$/ && $value == 0;
        $value = "An empty string (\"\")" if $value eq "";
    }
    else {
        $value ||= $DEFAULT_VALUE;
    }
    $params{-value} = $value;

    my $self = $class->SUPER::new( %params );
    return $self;
}

=head2 pretty_format()

 Purpose : Get a nicely formatted string containing information about the
           exception. Format is similar to that produced by
           Bio::Root::Root::throw(), with the addition of the name of
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
    my $value_string = $self->value ne $DEFAULT_VALUE ? "VALUE: ".$self->value."\n" : "";
    my $class = ref($self);

    my $title = "------------- EXCEPTION: $class -------------";
    my $footer = "\n" . '-' x CORE::length($title);
    my $out = "\n$title\n"
            . "MSG: $msg\n". $value_string. $stack. $footer . "\n";
    return $out;
}


=head2 _reformat_stacktrace

Reformatting of the stack performed by  _reformat_stacktrace:
for :list
1. Shift the file:line data in line i to line i+1.
2. change xxx::__ANON__() to "try{} block"
3. skip the "require" and "Error::subs::try" stack entries (boring)

This means that the first line in the stack won't have any file:line data
But this isn't a big issue since it's for a Bio::Root::-based method
that doesn't vary from exception to exception.

=cut

sub _reformat_stacktrace {
    my $self = shift;
    my $msg = $self->text;
    my $stack = $self->stacktrace();
    $stack =~ s/\Q$msg//;
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
    push @new_stack, "STACK: $prev_file:$prev_linenum";

    return join "\n", @new_stack;
}

=head2 stringify()

 Purpose : Overrides Error::stringify() to call pretty_format().
           This is called automatically when an exception object
           is placed between double quotes.
 Example : catch Bio::Root::Exception with {
              my $error = shift;
              print "$error";
           }

See Also: L<pretty_format()|pretty_format>

=cut

sub stringify {
    my ($self, @args) = @_;
    return $self->pretty_format( @args );
}

=head1 Subclasses of Bio::Root::Exception

=head2 L<Bio::Root::NotImplemented>

 Purpose : Indicates that a method has not been implemented.
 Example : throw Bio::Root::NotImplemented(
               -text   => "Method \"foo\" not implemented in module FooBar.",
               -value  => "foo" );

=cut

#---------------------------------------------------------
@Bio::Root::NotImplemented::ISA = qw( Bio::Root::Exception );
#---------------------------------------------------------

=head2 L<Bio::Root::IOException>

 Purpose : Indicates that some input/output-related trouble has occurred.
 Example : throw Bio::Root::IOException(
               -text   => "Can't save data to file $file.",
               -value  => $! );

=cut

#---------------------------------------------------------
@Bio::Root::IOException::ISA = qw( Bio::Root::Exception );
#---------------------------------------------------------


=head2 L<Bio::Root::FileOpenException>

 Purpose : Indicates that a file could not be opened.
 Example : throw Bio::Root::FileOpenException(
               -text   => "Can't open file $file for reading.",
               -value  => $! );

=cut

#---------------------------------------------------------
@Bio::Root::FileOpenException::ISA = qw( Bio::Root::IOException );
#---------------------------------------------------------


=head2 L<Bio::Root::SystemException>

 Purpose : Indicates that a system call failed.
 Example : unlink($file) or throw Bio::Root::SystemException(
               -text   => "Can't unlink file $file.",
               -value  => $! );

=cut

#---------------------------------------------------------
@Bio::Root::SystemException::ISA = qw( Bio::Root::Exception );
#---------------------------------------------------------


=head2 L<Bio::Root::BadParameter>

 Purpose : Indicates that one or more parameters supplied to a method
           are invalid, unspecified, or conflicting.
 Example : throw Bio::Root::BadParameter(
               -text   => "Required parameter \"-foo\" was not specified",
               -value  => "-foo" );

=cut

#---------------------------------------------------------
@Bio::Root::BadParameter::ISA = qw( Bio::Root::Exception );
#---------------------------------------------------------


=head2 L<Bio::Root::OutOfRange>

 Purpose : Indicates that a specified (start,end) range or
           an index to an array is outside the permitted range.
 Example : throw Bio::Root::OutOfRange(
               -text   => "Start coordinate ($start) cannot be less than zero.",
               -value  => $start  );

=cut

#---------------------------------------------------------
@Bio::Root::OutOfRange::ISA = qw( Bio::Root::Exception );
#---------------------------------------------------------


=head2 L<Bio::Root::NoSuchThing>

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
