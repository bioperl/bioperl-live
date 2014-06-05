package Bio::Root::Root;
use strict;
use Bio::Root::IO;
use Scalar::Util qw(blessed reftype);
use base qw(Bio::Root::RootI);

# ABSTRACT: hash-based implementation of L<Bio::Root::RootI>
# AUTHOR:   Steve Chervitz <sac@bioperl.org>
# AUTHOR:   Ewan Birney
# AUTHOR:   Lincoln Stein
# OWNER:    Steve Chervitz
# OWNER:    Ewan Birney
# OWNER:    Lincoln Stein
# LICENSE:  Perl_5

=head1 SYNOPSIS

  # Any Bioperl-compliant object is a RootI compliant object

  # Here's how to throw and catch an exception using the eval-based syntax.

  $obj->throw("This is an exception");

  eval {
      $obj->throw("This is catching an exception");
  };

  if( $@ ) {
      print "Caught exception";
  } else {
      print "no exception";
  }

  # Alternatively, using the new typed exception syntax in the throw() call:

  $obj->throw( -class => 'Bio::Root::BadParameter',
               -text  => "Can not open file $file",
               -value  => $file );

  # Want to see debug() outputs for this object

  my $obj = Bio::Object->new(-verbose=>1);

  my $obj = Bio::Object->new(%args);
  $obj->verbose(2);

  # Print debug messages which honour current verbosity setting

  $obj->debug("Boring output only to be seen if verbose > 0\n");

  # Deep-object copy

  my $clone = $obj->clone;

=head1 DESCRIPTION

This is a hashref-based implementation of the Bio::Root::RootI
interface.  Most Bioperl objects should inherit from this.

See the documentation for L<Bio::Root::RootI> for most of the methods
implemented by this module.  Only overridden methods are described
here.

=head2 Throwing Exceptions

One of the functionalities that L<Bio::Root::RootI> provides is the
ability to L<throw>() exceptions with pretty stack traces. Bio::Root::Root
enhances this with the ability to use L<Error> (available from CPAN)
if it has also been installed.

If L<Error> has been installed, L<throw>() will use it. This causes an
Error.pm-derived object to be thrown. This can be caught within a
C<catch{}> block, from wich you can extract useful bits of
information. If L<Error> is not installed, it will use the
L<Bio::Root::RootI>-based exception throwing facilty.

=head2 Typed Exception Syntax

The typed exception syntax of L<throw>() has the advantage of plainly
indicating the nature of the trouble, since the name of the class
is included in the title of the exception output.

To take advantage of this capability, you must specify arguments
as named parameters in the L<throw>() call. Here are the parameters:

=over 4

=item -class

name of the class of the exception.
This should be one of the classes defined in L<Bio::Root::Exception>,
or a custom error of yours that extends one of the exceptions
defined in L<Bio::Root::Exception>.

=item -text

a sensible message for the exception

=item -value

the value causing the exception or $!, if appropriate.

=back

Note that Bio::Root::Exception does not need to be imported into
your module (or script) namespace in order to throw exceptions
via Bio::Root::Root::throw(), since Bio::Root::Root imports it.

=head2 Try-Catch-Finally Support

In addition to using an eval{} block to handle exceptions, you can
also use a try-catch-finally block structure if L<Error> has been
installed in your system (available from CPAN).  See the documentation
for Error for more details.

Here's an example. See the L<Bio::Root::Exception> module for
other pre-defined exception types:

   my $IN;
   try {
    open $IN, '<', $file or $obj->throw( -class => 'Bio::Root::FileOpenException',
                                         -text  => "Cannot read file '$file'",
                                         -value => $!);
   }
   catch Bio::Root::BadParameter with {
       my $err = shift;   # get the Error object
       # Perform specific exception handling code for the FileOpenException
   }
   catch Bio::Root::Exception with {
       my $err = shift;   # get the Error object
       # Perform general exception handling code for any Bioperl exception.
   }
   otherwise {
       # A catch-all for any other type of exception
   }
   finally {
       # Any code that you want to execute regardless of whether or not
       # an exception occurred.
   };
   # the ending semicolon is essential!

=cut

our ($DEBUG, $ID, $VERBOSITY, $ERRORLOADED, $CLONE_CLASS);

BEGIN {
    $ID        = 'Bio::Root::Root';
    $DEBUG     = 0;
    $VERBOSITY = 0;
    $ERRORLOADED = 0;

    # Check whether or not Error.pm is available.

    # $main::DONT_USE_ERROR is intended for testing purposes and also
    # when you don't want to use the Error module, even if it is installed.
    # Just put a INIT { $DONT_USE_ERROR = 1; } at the top of your script.
    if( not $main::DONT_USE_ERROR ) {
        if ( eval "require Error; 1;"  ) {
            import Error qw(:try);
            require Bio::Root::Exception;
            $ERRORLOADED = 1;
            $Error::Debug = 1; # enable verbose stack trace
        }
    }
    if( !$ERRORLOADED ) {
        require Carp; import Carp qw( confess );
    }

    # set up _dclone()
    for my $class (qw(Clone Storable)) {
        eval "require $class; 1;";
        if (!$@) {
            $CLONE_CLASS = $class;
            if ($class eq 'Clone') {
                *Bio::Root::Root::_dclone = sub {shift; return Clone::clone(shift)};
            } else {
                *Bio::Root::Root::_dclone = sub {
                    shift;
                    local $Storable::Deparse = 1;
                    local $Storable::Eval = 1;
                    return Storable::dclone(shift);
                };
            }
            last;
        }
    }
    if (!defined $CLONE_CLASS) {
        *Bio::Root::Root::_dclone = sub {
            my ($self, $orig, $level) = @_;
            my $class = Scalar::Util::blessed($orig) || '';
            my $reftype = Scalar::Util::reftype($orig) || '';
            my $data;
            if (!$reftype) {
                $data = $orig
            } elsif ($reftype eq "ARRAY") {
                $data = [map $self->_dclone($_), @$orig];
            } elsif ($reftype eq "HASH") {
                $data = { map { $_ => $self->_dclone($orig->{$_}) } keys %$orig };
            } elsif ($reftype eq 'CODE') { # nothing, maybe shallow copy?
                $self->throw("Code reference cloning not supported; install Clone or Storable from CPAN");
            } else { $self->throw("What type is $_?")}
            if ($class) {
                bless $data, $class;
            }
            $data;
        }
    }

    $main::DONT_USE_ERROR;  # so that perl -w won't warn "used only once"
}

=head2 new

 Purpose   : generic instantiation function can be overridden if
             special needs of a module cannot be done in _initialize

=cut

sub new {
#    my ($class, %param) = @_;
    my $class = shift;
    my $self = {};
    bless $self, ref($class) || $class;

    if(@_ > 1) {
        # if the number of arguments is odd but at least 3, we'll give
        # it a try to find -verbose
        shift if @_ % 2;
        my %param = @_;
        ## See "Comments" above regarding use of _rearrange().
        $self->verbose($param{'-VERBOSE'} || $param{'-verbose'});
    }
    return $self;
}


=head2 clone

 Title   : clone
 Usage   : my $clone = $obj->clone();
           or
           my $clone = $obj->clone( -start => 110 );
 Function: Deep recursion copying of any object via Storable dclone()
 Returns : A cloned object.
 Args    : Any named parameters provided will be set on the new object.
           Unnamed parameters are ignored.
 Comments: Where possible, faster clone methods are used, in order:
           Clone::Fast::clone(), Clone::clone(), Storable::dclone.  If neither
           is present, a pure perl fallback (not very well tested) is used
           instead. Storable dclone() cannot clone CODE references.  Therefore,
           any CODE reference in your original object will remain, but will not
           exist in the cloned object.  This should not be used for anything
           other than cloning of simple objects. Developers of subclasses are
           encouraged to override this method with one of their own.

=cut

sub clone {
    my ($orig, %named_params) = @_;

    __PACKAGE__->throw("Can't call clone() as a class method") unless
        ref $orig && $orig->isa('Bio::Root::Root');

    # Can't dclone CODE references...
    # Should we shallow copy these? Should be harmless for these specific
    # methods...

    my %put_these_back = (
       _root_cleanup_methods => $orig->{'_root_cleanup_methods'},
    );
    delete $orig->{_root_cleanup_methods};

    # call the proper clone method, set lazily above
    my $clone = __PACKAGE__->_dclone($orig);

    $orig->{_root_cleanup_methods} = $put_these_back{_root_cleanup_methods};

    foreach my $key (grep { /^-/ } keys %named_params) {
        my $method = $key;
        $method =~ s/^-//;
        if ($clone->can($method)) {
            $clone->$method($named_params{$key})
        } else {
            $orig->warn("Parameter $method is not a method for ".ref($clone));
        }
    }
    return $clone;
}

=head2 _dclone

 Title   : clone
 Usage   : my $clone = $obj->_dclone($ref);
           or
           my $clone = $obj->_dclone($ref);
 Function: Returns a copy of the object passed to it (a deep clone)
 Returns : clone of passed argument
 Args    : Anything
 NOTE    : This differs from clone significantly in that it does not clone
           self, but the data passed to it.  This code may need to be optimized
           or overridden as needed.
 Comments: This is set in the BEGIN block to take advantage of optimized
           cloning methods if Clone or Storable is present, falling back to a
           pure perl kludge. May be moved into a set of modules if the need
           arises. At the moment, code ref cloning is not supported.

=cut

=head2 verbose

 Title   : verbose
 Usage   : $self->verbose(1)
 Function: Sets verbose level for how ->warn behaves
           -1 = no warning
            0 = standard, small warning
            1 = warning with stack trace
            2 = warning becomes throw
 Returns : The current verbosity setting (integer between -1 to 2)
 Args    : -1,0,1 or 2


=cut

sub verbose {
    my ($self,$value) = @_;
    # allow one to set global verbosity flag
    return $DEBUG  if $DEBUG;
    return $VERBOSITY unless ref $self;

    if (defined $value || ! defined $self->{'_root_verbose'}) {
        $self->{'_root_verbose'} = $value || 0;
    }
    return $self->{'_root_verbose'};
}

=head2 _register_for_cleanup

=cut

sub _register_for_cleanup {
    my ($self,$method) = @_;
    if ($method) {
        if(! exists($self->{'_root_cleanup_methods'})) {
            $self->{'_root_cleanup_methods'} = [];
        }
        push(@{$self->{'_root_cleanup_methods'}},$method);
    }
}

=head2 _unregister_for_cleanup

=cut

sub _unregister_for_cleanup {
    my ($self,$method) = @_;
    my @methods = grep {$_ ne $method} $self->_cleanup_methods;
    $self->{'_root_cleanup_methods'} = \@methods;
}

=head2 _cleanup_methods

=cut

sub _cleanup_methods {
    my $self = shift;
    return unless ref $self && $self->isa('HASH');
    my $methods = $self->{'_root_cleanup_methods'} or return;
    @$methods;
}

=head2 throw

 Title   : throw
 Usage   : $obj->throw("throwing exception message");
           or
           $obj->throw( -class => 'Bio::Root::Exception',
                        -text  => "throwing exception message",
                        -value => $bad_value  );
 Function: Throws an exception, which, if not caught with an eval or
           a try block will provide a nice stack trace to STDERR
           with the message.
           If Error.pm is installed, and if a -class parameter is
           provided, Error::throw will be used, throwing an error
           of the type specified by -class.
           If Error.pm is installed and no -class parameter is provided
           (i.e., a simple string is given), A Bio::Root::Exception
           is thrown.
 Returns : n/a
 Args    : A string giving a descriptive error message, optional
           Named parameters:
           '-class'  a string for the name of a class that derives
                     from Error.pm, such as any of the exceptions
                     defined in Bio::Root::Exception.
                     Default class: Bio::Root::Exception
           '-text'   a string giving a descriptive error message
           '-value'  the value causing the exception, or $! (optional)

           Thus, if only a string argument is given, and Error.pm is available,
           this is equivalent to the arguments:
                 -text  => "message",
                 -class => Bio::Root::Exception
 Comments : If Error.pm is installed, and you don't want to use it
            for some reason, you can block the use of Error.pm by
            Bio::Root::Root::throw() by defining a scalar named
            $main::DONT_USE_ERROR (define it in your main script
            and you don't need the main:: part) and setting it to
            a true value; you must do this within a BEGIN subroutine.

=cut

sub throw {
    my ($self, @args) = @_;

    my ($text, $class, $value) = $self->_rearrange( [qw(TEXT
                                                        CLASS
                                                        VALUE)], @args);
    $text ||= $args[0] if @args == 1;

    if ($ERRORLOADED) {
        # Enable re-throwing of Error objects.
        # If the error is not derived from Bio::Root::Exception,
        # we can't guarantee that the Error's value was set properly
        # and, ipso facto, that it will be catchable from an eval{}.
        # But chances are, if you're re-throwing non-Bio::Root::Exceptions,
        # you're probably using Error::try(), not eval{}.
        # TODO: Fix the MSG: line of the re-thrown error. Has an extra line
        # containing the '----- EXCEPTION -----' banner.
        if (ref($args[0])) {
            if( $args[0]->isa('Error')) {
                my $class = ref $args[0];
                $class->throw( @args );
            }
            else {
                my $text .= "\nWARNING: Attempt to throw a non-Error.pm object: " . ref$args[0];
                my $class = "Bio::Root::Exception";
                $class->throw( '-text' => $text, '-value' => $args[0] );
            }
        }
        else {
            $class ||= "Bio::Root::Exception";

            my %args;
            if( @args % 2 == 0 && $args[0] =~ /^-/ ) {
                %args = @args;
                $args{-text} = $text;
                $args{-object} = $self;
            }

            $class->throw( scalar keys %args > 0 ? %args : @args ); # (%args || @args) puts %args in scalar context!
        }
    }
    else {
        $class ||= '';
        $class = ': '.$class if $class;
        my $std = $self->stack_trace_dump();
        my $title = "------------- EXCEPTION$class -------------";
        my $footer = ('-' x CORE::length($title))."\n";
        $text ||= '';

        die "\n$title\n", "MSG: $text\n", $std, $footer, "\n";
    }
}

=head2 debug

 Title   : debug
 Usage   : $obj->debug("This is debugging output");
 Function: Prints a debugging message when verbose is > 0
 Returns : none
 Args    : message string(s) to print to STDERR

=cut

sub debug {
    my ($self, @msgs) = @_;

    # using CORE::warn doesn't give correct backtrace information; we want the
    # line from the previous call in the call stack, not this call (similar to
    # cluck).  For now, just add a stack trace dump and simple comment under the
    # correct conditions.
    if (defined $self->verbose && $self->verbose > 0) {
        if (!@msgs || $msgs[-1] !~ /\n$/) {
            push @msgs, "Debugging comment:" if !@msgs;
            push @msgs, sprintf("%s %s:%s", @{($self->stack_trace)[2]}[3,1,2])."\n";
        }
        CORE::warn @msgs;
    }
}

=head2 _load_module

 Title   : _load_module
 Usage   : $self->_load_module("Bio::SeqIO::genbank");
 Function: Loads up (like use) the specified module at run time on demand.
 Example :
 Returns : TRUE on success. Throws an exception upon failure.
 Args    : The module to load (_without_ the trailing .pm).

=cut

sub _load_module {
    my ($self, $name) = @_;
    my ($module, $load, $m);
    $module = "_<$name.pm";
    return 1 if $main::{$module};

    # untaint operation for safe web-based running (modified after
    # a fix by Lincoln) HL
    if ($name !~ /^([\w:]+)$/) {
        $self->throw("$name is an illegal perl package name");
    } else {
        $name = $1;
    }

    $load = "$name.pm";
    my $io = Bio::Root::IO->new();
    # catfile comes from IO
    $load = $io->catfile((split(/::/,$load)));
    eval {
        require $load;
    };
    if ( $@ ) {
        $self->throw("Failed to load module $name. ".$@);
    }
    return 1;
}

=head2 DESTROY

=cut

sub DESTROY {
    my $self = shift;
    my @cleanup_methods = $self->_cleanup_methods or return;
    for my $method (@cleanup_methods) {
        $method->($self);
    }
}

1;
