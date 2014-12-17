package Bio::Root::RootI;
use strict;
use Carp 'confess','carp';

# ABSTRACT: abstract interface to root object code
# AUTHOR:   Steve Chervitz <sac@bioperl.org>
# AUTHOR:   Ewan Birney <birney@ebi.ac.uk>
# AUTHOR:   Lincoln Stein
# OWNER:    Steve Chervitz
# OWNER:    Ewan Birney
# OWNER:    Lincoln Stein
# LICENSE:  Perl_5

# CONTRIBUTOR: Sendu Bala <bix@sendu.me.uk>
# CONTRIBUTOR: Jason Stajich

=head1 SYNOPSIS

  # any bioperl or bioperl compliant object is a RootI
  # compliant object

  $obj->throw("This is an exception");

  eval {
      $obj->throw("This is catching an exception");
  };

  if( $@ ) {
      print "Caught exception";
  } else {
      print "no exception";
  }

  # Using throw_not_implemented() within a RootI-based interface module:

  package Foo;
  use base qw(Bio::Root::RootI);

  sub foo {
      my $self = shift;
      $self->throw_not_implemented;
  }


=head1 DESCRIPTION

This is just a set of methods which do not assume B<anything> about the object
they are on. The methods provide the ability to throw exceptions with nice
stack traces.

This is what should be inherited by all Bioperl compliant interfaces, even
if they are exotic XS/CORBA/Other perl systems.

=head2 Using throw_not_implemented()

The method L<throw_not_implemented()|throw_not_implemented> should be
called by all methods within interface modules that extend RootI so
that if an implementation fails to override them, an exception will be
thrown.

For example, say there is an interface module called C<FooI> that
provides a method called C<foo()>. Since this method is considered
abstract within FooI and should be implemented by any module claiming to
implement C<FooI>, the C<FooI::foo()> method should consist of the
following:

    sub foo {
        my $self = shift;
        $self->throw_not_implemented;
    }

So, if an implementer of C<FooI> forgets to implement C<foo()>
and a user of the implementation calls C<foo()>, a
L<Bio::Exception::NotImplemented> exception will result.

Unfortunately, failure to implement a method can only be determined at
run time (i.e., you can't verify that an implementation is complete by
running C<perl -wc> on it). So it should be standard practice for a test
of an implementation to check each method and verify that it doesn't
throw a L<Bio::Exception::NotImplemented>.

=cut

use vars qw($DEBUG $ID $VERBOSITY);
BEGIN {
    $ID        = 'Bio::Root::RootI';
    $DEBUG     = 0;
    $VERBOSITY = 0;
}

=head2 new

=cut

sub new {
    my $class = shift;
    my @args = @_;
    unless ( $ENV{'BIOPERLDEBUG'} ) {
        carp("Use of new in Bio::Root::RootI is deprecated.  Please use Bio::Root::Root instead");
    }
    eval "require Bio::Root::Root";
    return Bio::Root::Root->new(@args);
}

# for backwards compatibility
sub _initialize {
    my($self,@args) = @_;
    return 1;
}


=head2 throw

 Title   : throw
 Usage   : $obj->throw("throwing exception message")
 Function: Throws an exception, which, if not caught with an eval brace
           will provide a nice stack trace to STDERR with the message
 Returns : nothing
 Args    : A string giving a descriptive error message


=cut

sub throw{
    my ($self,$string) = @_;

    my $std = $self->stack_trace_dump();

    my $out = "\n-------------------- EXCEPTION --------------------\n"
            . "MSG: " . $string . "\n"
            . $std."-------------------------------------------\n";
    die $out;
}

=head2 warn

 Title   : warn
 Usage   : $object->warn("Warning message");
 Function: Places a warning. What happens now is down to the
           verbosity of the object  (value of $obj->verbose)
            verbosity 0 or not set => small warning
            verbosity -1 => no warning
            verbosity 1 => warning with stack trace
            verbosity 2 => converts warnings into throw
 Returns : n/a
 Args    : string (the warning message)

=cut

sub warn {
    my ($self,$string) = @_;

    my $verbose = $self->verbose;

    my $header = "\n--------------------- WARNING ---------------------\nMSG: ";
    my $footer =   "---------------------------------------------------\n";

    if ($verbose >= 2) {
        $self->throw($string);
    }
    elsif ($verbose <= -1) {
        return;
    }
    elsif ($verbose == 1) {
        CORE::warn $header, $string, "\n", $self->stack_trace_dump, $footer;
        return;
    }

    CORE::warn $header, $string, "\n", $footer;
}

=head2 deprecated

 Title   : deprecated
 Usage   : $obj->deprecated("Method X is deprecated");
           $obj->deprecated("Method X is deprecated", 1.007);
           $obj->deprecated(-message => "Method X is deprecated");
           $obj->deprecated(-message => "Method X is deprecated",
                            -version => 1.007);
 Function: Prints a message about deprecation unless verbose is < 0
           (which means be quiet)
 Returns : none
 Args    : Message string to print to STDERR
           Version of BioPerl where use of the method results in an exception
 Notes   : The method can be called two ways, either by positional arguments:

           $obj->deprecated('This module is deprecated', 1.006);

           or by named arguments:

           $obj->deprecated(
                -message => 'use of the method foo() is deprecated, use bar() instead',
                -version => 1.006  # throw if $VERSION is >= this version
                );

           or timed to go off at a certain point:

           $obj->deprecated(
                -message => 'use of the method foo() is deprecated, use bar() instead',
                -warn_version    => 1.006 # warn if $VERSION is >= this version
                -throw_version   => 1.007 # throw if $VERSION is >= this version
                );

           Using the last two named argument versions is suggested and will
           likely be the only supported way of calling this method in the future
           Yes, we see the irony of deprecating that particular usage of
           deprecated().

           The main difference between usage of the two named argument versions
           is that by designating a 'warn_version' one indicates the
           functionality is officially deprecated beginning in a future version
           of BioPerl (so warnings are issued only after that point), whereas
           setting either 'version' or 'throw_version' (synonyms) converts the
           deprecation warning to an exception.

           For proper comparisons one must use a version in lines with the
           current versioning scheme for Perl and BioPerl, (i.e. where 1.006000
           indicates v1.6.0, 5.010000 for v5.10.0, etc.).

=cut

sub deprecated{
    my ($self) = shift;

    my $class = ref $self || $self;
    my $class_version = do {
        no strict 'refs';
        ${"${class}::VERSION"}
    };

    if( $class_version && $class_version =~ /set by/ ) {
        $class_version = 0.0001;
    }

    my ($msg, $version, $warn_version, $throw_version) =
        $self->_rearrange([qw(MESSAGE VERSION WARN_VERSION THROW_VERSION)], @_);

    $throw_version ||= $version;
    $warn_version  ||= $class_version;

    for my $v ( $warn_version, $throw_version) {
        no warnings 'numeric';
        $self->throw("Version must be numerical, such as 1.006000 for v1.6.0, not $v")
            unless !defined $v || $v + 0 eq $v;
    }

    # below default insinuates we're deprecating a method and not a full module
    # but it's the most common use case
    $msg ||= "Use of ".(caller(1))[3]."() is deprecated.";

    if( $throw_version && $class_version && $class_version >= $throw_version ) {
        $self->throw($msg)
    }
    elsif( $warn_version && $class_version && $class_version >= $warn_version ) {

        $msg .= "\nTo be removed in $throw_version." if $throw_version;

        # passing this on to warn() should deal properly with verbosity issues
        $self->warn($msg);
    }
}

=head2 stack_trace_dump

 Title   : stack_trace_dump
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


=cut

sub stack_trace_dump{
    my ($self) = @_;

    my @stack = $self->stack_trace();

    shift @stack;
    shift @stack;
    shift @stack;

    my $out;
    my ($module,$function,$file,$position);


    foreach my $stack ( @stack) {
        ($module,$file,$position,$function) = @{$stack};
        $out .= "STACK $function $file:$position\n";
    }

    return $out;
}


=head2 stack_trace

 Title   : stack_trace
 Usage   : @stack_array_ref= $self->stack_trace
 Function: gives an array to a reference of arrays with stack trace info
           each coming from the caller(stack_number) call
 Returns : array containing a reference of arrays
 Args    : none


=cut

sub stack_trace{
    my ($self) = @_;

    my $i = 0;
    my @out = ();
    my $prev = [];
    while( my @call = caller($i++)) {
        # major annoyance that caller puts caller context as
        # function name. Hence some monkeying around...
        $prev->[3] = $call[3];
        push(@out,$prev);
        $prev = \@call;
    }
    $prev->[3] = 'toplevel';
    push(@out,$prev);
    return @out;
}


=head2 _rearrange

 Usage     : $object->_rearrange( array_ref, list_of_arguments)
 Purpose   : Rearranges named parameters to requested order.
 Example   : $self->_rearrange([qw(SEQUENCE ID DESC)],@param);
           : Where @param = (-sequence => $s,
           :                 -desc     => $d,
           :                 -id       => $i);
 Returns   : @params - an array of parameters in the requested order.
           : The above example would return ($s, $i, $d).
           : Unspecified parameters will return undef. For example, if
           :        @param = (-sequence => $s);
           : the above _rearrange call would return ($s, undef, undef)
 Argument  : $order : a reference to an array which describes the desired
           :          order of the named parameters.
           : @param : an array of parameters, either as a list (in
           :          which case the function simply returns the list),
           :          or as an associative array with hyphenated tags
           :          (in which case the function sorts the values
           :          according to @{$order} and returns that new array.)
           :          The tags can be upper, lower, or mixed case
           :          but they must start with a hyphen (at least the
           :          first one should be hyphenated.)
 Source    : This function was taken from CGI.pm, written by Dr. Lincoln
           : Stein, and adapted for use in Bio::Seq by Richard Resnick and
           : then adapted for use in Bio::Root::Object.pm by Steve Chervitz,
           : then migrated into Bio::Root::RootI.pm by Ewan Birney.
 Comments  :
           : Uppercase tags are the norm,
           : (SAC)
           : This method may not be appropriate for method calls that are
           : within in an inner loop if efficiency is a concern.
           :
           : Parameters can be specified using any of these formats:
           :  @param = (-name=>'me', -color=>'blue');
           :  @param = (-NAME=>'me', -COLOR=>'blue');
           :  @param = (-Name=>'me', -Color=>'blue');
           :  @param = ('me', 'blue');
           : A leading hyphenated argument is used by this function to
           : indicate that named parameters are being used.
           : Therefore, the ('me', 'blue') list will be returned as-is.
           :
           : Note that Perl will confuse unquoted, hyphenated tags as
           : function calls if there is a function of the same name
           : in the current namespace:
           :    -name => 'foo' is interpreted as -&name => 'foo'
           :
           : For ultimate safety, put single quotes around the tag:
           : ('-name'=>'me', '-color' =>'blue');
           : This can be a bit cumbersome and I find not as readable
           : as using all uppercase, which is also fairly safe:
           : (-NAME=>'me', -COLOR =>'blue');
           :
           : Personal note (SAC): I have found all uppercase tags to
           : be more manageable: it involves less single-quoting,
           : the key names stand out better, and there are no method naming
           : conflicts.
           : The drawbacks are that it's not as easy to type as lowercase,
           : and lots of uppercase can be hard to read.
           :
           : Regardless of the style, it greatly helps to line
           : the parameters up vertically for long/complex lists.
           :
           : Note that if @param is a single string that happens to start with
           : a dash, it will be treated as a hash key and probably fail to
           : match anything in the array_ref, so not be returned as normally
           : happens when @param is a simple list and not an associative array.

=cut

sub _rearrange {
    my ($self, $order, @args) = @_;

    return @args unless $args[0] && $args[0] =~ /^\-/;

    push @args, undef unless $#args % 2;

    my %param;
    for( my $i = 0; $i < @args; $i += 2 ) {
        (my $key = $args[$i]) =~ tr/a-z\055/A-Z/d; #deletes all dashes!
        $param{$key} = $args[$i+1];
    }
    return @param{map uc, @$order};
}

=head2 _set_from_args

 Usage     : $object->_set_from_args(\%args, -methods => \@methods)
 Purpose   : Takes a hash of user-supplied args whose keys match method names,
           : and calls the method supplying it the corresponding value.
 Example   : $self->_set_from_args(\%args, -methods => [qw(sequence id desc)]);
           : Where %args = (-sequence    => $s,
           :                -description => $d,
           :                -ID          => $i);
           :
           : the above _set_from_args calls the following methods:
           : $self->sequence($s);
           : $self->id($i);
           : ( $self->description($i) is not called because 'description' wasn't
           :   one of the given methods )
 Argument  : \%args | \@args : a hash ref or associative array ref of arguments
           :                   where keys are any-case strings corresponding to
           :                   method names but optionally prefixed with
           :                   hyphens, and values are the values the method
           :                   should be supplied. If keys contain internal
           :                   hyphens (eg. to separate multi-word args) they
           :                   are converted to underscores, since method names
           :                   cannot contain dashes.
           : -methods => []  : (optional) only call methods with names in this
           :                   array ref. Can instead supply a hash ref where
           :                   keys are method names (of real existing methods
           :                   unless -create is in effect) and values are array
           :                   refs of synonyms to allow access to the method
           :                   using synonyms. If there is only one synonym it
           :                   can be supplied as a string instead of a single-
           :                   element array ref
           : -force => bool  : (optional, default 0) call methods that don't
           :                   seem to exist, ie. let AUTOLOAD handle them
           : -create => bool : (optional, default 0) when a method doesn't
           :                   exist, create it as a simple getter/setter
           :                   (combined with -methods it would create all the
           :                   supplied methods that didn't exist, even if not
           :                   mentioned in the supplied %args)
           : -code => '' | {}: (optional) when creating methods use the supplied
           :                   code (a string which will be evaulated as a sub).
           :                   The default code is a simple get/setter.
           :                   Alternatively you can supply a hash ref where
           :                   the keys are method names and the values are
           :                   code strings. The variable '$method' will be
           :                   available at evaluation time, so can be used in
           :                   your code strings. Beware that the strict pragma
           :                   will be in effect.
           : -case_sensitive => bool : require case sensitivity on the part of
           :                           user (ie. a() and A() are two different
           :                           methods and the user must be careful
           :                           which they use).
 Comments  :
           : The \%args argument will usually be the args received during new()
           : from the user. The user is allowed to get the case wrong, include
           : 0 or more than one hyphens as a prefix, and to include hyphens as
           : multi-word arg separators: '--an-arg' => 1, -an_arg => 1 and
           : An_Arg => 1 are all equivalent, calling an_arg(1). However, in
           : documentation users should only be told to use the standard form
           : -an_arg to avoid confusion. A possible exception to this is a
           : wrapper module where '--an-arg' is what the user is used to
           : supplying to the program being wrapped.
           :
           : Another issue with wrapper modules is that there may be an
           : argument that has meaning both to Bioperl and to the program, eg.
           : -verbose. The recommended way of dealing with this is to leave
           : -verbose to set the Bioperl verbosity whilst requesting users use
           : an invented -program_verbose (or similar) to set the program
           : verbosity. This can be resolved back with
           : Bio::Tools::Run::WrapperBase's _setparams() method and code along
           : the lines of:
           : my %methods = map { $_ => $_ } @LIST_OF_ALL_ALLOWED_PROGRAM_ARGS
           : delete $methods{'verbose'};
           : $methods{'program_verbose'} = 'verbose';
           : my $param_string = $self->_setparams(-methods => \%methods);
           : system("$exe $param_string");

=cut

sub _set_from_args {
    my ($self, $args, @own_args) = @_;
    $self->throw("a hash/array ref of arguments must be supplied") unless ref($args);

    my ($methods, $force, $create, $code, $case);
    if (@own_args) {
        ($methods, $force, $create, $code, $case) =
            $self->_rearrange([qw(METHODS
                                  FORCE
                                  CREATE
                                  CODE
                                  CASE_SENSITIVE)], @own_args);
    }
    my $default_code = 'my $self = shift;
                        if (@_) { $self->{\'_\'.$method} = shift }
                        return $self->{\'_\'.$method};';

    my %method_names = ();
    my %syns = ();
    if ($methods) {
        my @names;
        if (ref($methods) eq 'HASH') {
            @names = keys %{$methods};
            %syns = %{$methods};
        }
        else {
            @names = @{$methods};
            %syns = map { $_ => $_ } @names;
        }
        %method_names = map { $case ? $_ : lc($_) => $_ } @names;
    }

    # deal with hyphens
    my %orig_args = ref($args) eq 'HASH' ? %{$args} : @{$args};
    my %args;
    while (my ($method, $value) = each %orig_args) {
        $method =~ s/^-+//;
        $method =~ s/-/_/g;
        $args{$method} = $value;
    }

    # create non-existing methods on request
    if ($create) {
        unless ($methods) {
            %syns = map { $_ => $case ? $_ : lc($_) } keys %args;
        }

        foreach my $method (keys %syns) {
            $self->can($method) && next;

            my $string = $code || $default_code;
            if (ref($code) && ref($code) eq 'HASH') {
                $string = $code->{$method} || $default_code;
            }

            my $sub = eval "sub { $string }";
            $self->throw("Compilation error for $method : $@") if $@;

            no strict 'refs';
            *{ref($self).'::'.$method} = $sub;
        }
    }

    # create synonyms of existing methods
    while (my ($method, $syn_ref) = each %syns) {
        my $method_ref = $self->can($method) || next;

        foreach my $syn (@{ ref($syn_ref) ? $syn_ref : [$syn_ref] }) {
            next if $syn eq $method;
            $method_names{$case ? $syn : lc($syn)} = $syn;
            next if $self->can($syn);
            no strict 'refs';
            *{ref($self).'::'.$syn} = $method_ref;
        }
    }

    # set values for methods
    while (my ($method, $value) = each %args) {
        $method = $method_names{$case ? $method : lc($method)} || ($methods ? next : $method);
        $self->can($method) || next unless $force;
        $self->$method($value);
    }
}


=head2 _rearrange_old

=cut

#----------------'
sub _rearrange_old {
#----------------
    my($self,$order,@param) = @_;

    # JGRG -- This is wrong, because we don't want
    # to assign empty string to anything, and this
    # code is actually returning an array 1 less
    # than the length of @param:

    ## If there are no parameters, we simply wish to return
    ## an empty array which is the size of the @{$order} array.
    #return ('') x $#{$order} unless @param;

    # ...all we need to do is return an empty array:
    # return unless @param;

    # If we've got parameters, we need to check to see whether
    # they are named or simply listed. If they are listed, we
    # can just return them.

    # The mod test fixes bug where a single string parameter beginning with '-' gets lost.
    # This tends to happen in error messages such as: $obj->throw("-id not defined")
    return @param unless (defined($param[0]) && $param[0]=~/^-/o && ($#param % 2));

    # Tester
#    print "\n_rearrange() named parameters:\n";
#    my $i; for ($i=0;$i<@param;$i+=2) { printf "%20s => %s\n", $param[$i],$param[$i+1]; }; <STDIN>;

    # Now we've got to do some work on the named parameters.
    # The next few lines strip out the '-' characters which
    # preceed the keys, and capitalizes them.
    for (my $i=0;$i<@param;$i+=2) {
        $param[$i]=~s/^\-//;
        $param[$i]=~tr/a-z/A-Z/;
    }

    # Now we'll convert the @params variable into an associative array.
    # local($^W) = 0;  # prevent "odd number of elements" warning with -w.
    my(%param) = @param;

    # my(@return_array);

    # What we intend to do is loop through the @{$order} variable,
    # and for each value, we use that as a key into our associative
    # array, pushing the value at that key onto our return array.
    # my($key);

    #foreach (@{$order}) {
    # my($value) = $param{$key};
    # delete $param{$key};
    #push(@return_array,$param{$_});
    #}

    return @param{@{$order}};

#    print "\n_rearrange() after processing:\n";
#    my $i; for ($i=0;$i<@return_array;$i++) { printf "%20s => %s\n", ${$order}[$i], $return_array[$i]; } <STDIN>;

    # return @return_array;
}

=head2 _register_for_cleanup

 Title   : _register_for_cleanup
 Usage   : -- internal --
 Function: Register a method to be called at DESTROY time. This is useful
           and sometimes essential in the case of multiple inheritance for
           classes coming second in the sequence of inheritance.
 Returns :
 Args    : a code reference

The code reference will be invoked with the object as the first
argument, as per a method.  You may register an unlimited number of
cleanup methods.

=cut

sub _register_for_cleanup {
    my ($self,$method) = @_;
    $self->throw_not_implemented();
}

=head2 _unregister_for_cleanup

 Title   : _unregister_for_cleanup
 Usage   : -- internal --
 Function: Remove a method that has previously been registered to be called
           at DESTROY time.  If called with a method to be called at DESTROY time.
           Has no effect if the code reference has not previously been registered.
 Returns : nothing
 Args    : a code reference

=cut

sub _unregister_for_cleanup {
    my ($self,$method) = @_;
    $self->throw_not_implemented();
}

=head2 _cleanup_methods

 Title   : _cleanup_methods
 Usage   : -- internal --
 Function: Return current list of registered cleanup methods.
 Returns : list of coderefs
 Args    : none

=cut

sub _cleanup_methods {
    my $self = shift;
    unless ( $ENV{'BIOPERLDEBUG'} || $self->verbose  > 0 ) {
        carp("Use of Bio::Root::RootI is deprecated.  Please use Bio::Root::Root instead");
    }
    return;
}

=head2 throw_not_implemented

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
 Throws  : A Bio::Root::NotImplemented exception.
           The message of the exception contains
             - the name of the method
             - the name of the interface
             - the name of the implementing class

           If this object has a throw() method, $self->throw will be used.
           If the object doesn't have a throw() method,
           Carp::confess() will be used.


=cut

#'

sub throw_not_implemented {
    my $self = shift;

    # Bio::Root::Root::throw() knows how to check for Error.pm and will
    # throw an Error-derived object of the specified class (Bio::Root::NotImplemented),
    # which is defined in Bio::Root::Exception.
    # If Error.pm is not available, the name of the class is just included in the
    # error message.

    my $message = $self->_not_implemented_msg;

    if ( $self->can('throw') ) {
        my @args;
        if ( $self->isa('Bio::Root::Root') ) {
            # Use Root::throw() hash-based arguments instead of RootI::throw()
            # single string argument whenever possible
            @args = ( -text  => $message, -class => 'Bio::Root::NotImplemented' );
        } else {
            @args = ( $message );
        }
        $self->throw(@args);

    } else {
        confess $message;
    }
}


=head2 warn_not_implemented

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


=cut

#'

sub warn_not_implemented {
    my $self = shift;
    my $message = $self->_not_implemented_msg;
    if( $self->can('warn') ) {
        $self->warn( $message );
    }else {
        carp $message ;
    }
}

=head2 _not_implemented_msg

Unify 'not implemented' message. -Juguang
=cut

sub _not_implemented_msg {
    my $self = shift;
    my $package = ref $self;
    my $meth = (caller(2))[3];
    my $msg =<<EOD_NOT_IMP;
Abstract method \"$meth\" is not implemented by package $package.
This is not your fault - author of $package should be blamed!
EOD_NOT_IMP
    return $msg;
}

1;
