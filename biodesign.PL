#!/usr/bin/perl
use strict;
use LocalConfig;
use File::Basename qw(basename dirname);
my $file;
($file = basename($0)) =~ s/\.PL$/\.pod/;
open OUT,">$file" or die "Can't create $file: $!";

print OUT<<'!NO!GROK!';
=head1 NAME

biodesign.pod - Bioperl Design Documentation

=head1 SYNOPSIS

Not appropiate. Read on...

=head1 DESCRIPTION

Bioperl is a coordinated project which has a number of design features
which allow it to be well used and extensible as well as collaborate with
other packages. biodesign focuses on a number of areas.

  Bioperl etiquette and learning about it
  Bioperl root object - exception throwing, exceptions etc.
  Bioperl interface design
  Bioperl sequence object design notes

=head1 AUTHOR

This was written by Ewan Birney in a variety of airports across the US.

Additions added by:
Brian Osborne
Jason Stajich

=head1 Reusing code and working in collaborative projects

The biggest problem often in reusing a code base like Bioperl is that
it requires both the people using it and the people contributing to
it to change their attitude towards code. Generally people in bioinformatics
are more likely to be self-taught programmers who put together most
of their scripts/programs working alone. Bioperl is a truly collaborative
project (the core code is the product of about 15 individuals) and anyone
will be only contributing some part of it in the future.

Here are some notes about how my coding style has changed to work in
collaborative projects.

=head2 Learn to read documentation

Reading documentation is sometimes as tough as writing the
documentation. Try to read documentation before you ask a question -
not only might it answer your question, but more importantly it will
give you idea why the person who wrote the module wrote it - and this
will be the frame work in which you can understand his or her answer.

You might also want to examine the models, or class diagrams, in the models
directory. These diagrams are not guaranteed to include every single class
but may help you understand the overall layout of Bioperl's modules.

Documentation on Bio::Root can also be found in the form of scripts - check
the scripts/root_object directory for a start. We have also created a 
"howto" directory that addresses specific design questions, not
necessarily related to Bio::Root. See doc/howto.

=head2 Respect people's code (in particular if it works)

If the code does what you want, the fact that it is not written the
way you would write should not be a big issue. Of course, if there is
some glaring error then that is worth pointing out to
someone. Dismissing a module on the basis of its coding style is a
tremendously wrong thing to do.

=head2 Learn how to provide good feedback

This ranges from giving very accurate bug reports (this script --E<gt>
makes this error, giving all data), through to pointing out design
issues in a constructive manner (not - this *sucks*). If you find
a problem, then providing a patch using diff or a work around is
a great thing to do - the author/maintainer of the module will
love you for it.

Providing "I used XXX and it did just what I wanted it to do" feedback
is also really great. Developers generally only hear about their mistakes.
To hear about successes gives everyone a warm glow.

One trick we have learnt is that when we download a new project/code or
use a new module we open up a fresh buffer in emacs and keep a mini diary
of everything that we did or thought when we started to use the package. After
we used it we could go back, edit the buffer and then send it to the author
either with "it was great - it did just what I wanted, but I found that
the documentation here was misleading" to "to get it to install I had
to incant the following things..."

=head2 Taking on a project

When you want to get involved, hopefully it will be because you want to
extend something or provide better facilities to something. The important
thing here is not to work in a vacuum. Providing the main list with
a good proposal before you start about what you are going to do (and listen
to the responses) is a must. We have been pulled up so many times by other
people looking at our designs that we can't imagine coding stuff now without
feedback.

=head2 Designing good tests

Sadly, you might think that you have written good code, but you don't know
that until you manage to test it! The CPAN style perl modules have a wonderful
test suite system (delve around into the t/ directories) and we have extended
the makefile system so that the test script which you write to test the module
can be part of the t/ system from the start. Once a test is in the t/ system it
will be run millions of times worldwide when Bioperl is downloaded, providing
incredible and continual regression testing of your module (for free!).

=head2 Having fun

The coding process should be enjoyable, and we get very proud of
people who tell us that they picked up Bioperl and it worked for them,
even if they don't use a single module that we wrote. There is a
brilliant sense of community in Bioperl about providing useful, stable
code and it should be a pleasure to contribute to it.

So - we are always looking forward to people posting on the guts list
with their feedback/questions/proposals. As well as the long standing fun we
have making new releases.

=head1 Bioperl Bio::Root::Root Object

All objects in Bioperl (but for interfaces) inherit from the Bio::Root::Root.
The Bioperl root object allows a number of very useful concepts to be provided.
In particular.

=over 3

=item exceptions, warning, and debugging

The Bioperl Root object allows exceptions to be throw on the object
with very nice debugging output.  These are throw by calling the
method throw() and passing in the message string.  This will cause the
execution of the script to die with a stack trace.

Similarly the warn() method can be called which will produce a
warning message - use this instead of print for warning messages to
the user because if the verbose flag is set to -1 warnings will be
skipped.  Additionally setting the verbose flag to 1 will print a
stack trace for every warning in addition to the message and setting
verbose to 2 will convert warnings into thrown exceptions.

Finally, the the debug() method prints messages to STDERR when the
verbose flag is set to 1.

=item rearrange

Bioperl root object have some helper methods, in particular
rearrange() to help functions which take hash inputs.  This allows one
to specify named arguments as a hash and map them to the expected
input parameters specified by an array.

=back

You can go to L<Bio::Root::Root> for more information. There are also
a number of useful example scripts in the scripts/root_object directory.


=head2 Using the root object

To use the root object, the object has to inherit from it. This means
the @ISA array should have (Bio::Root::Root) in it and that the module
has a "use Bio::Root::Root". (If you are an emacs user, consider using
the boilerplate methods in the bioperl.lisp to lay out your module
initially for you).  The root object provides a top level -E<gt>new
function.  You should inherit from this new method by calling the new
method of the superclass which is accessible by using SUPER.  This is
called chaining the constructors and allows a child class to utilize
the initialization procedure of the superclass in addition to
executing its own.  This is a very powerful technique and allows Bioperl 
to behave in an Object Oriented manner.

The full code is given below for a basic skeleton object that uses bioperl.

  # convention is that if you are using the Bio::Root object you
  # should put it inside the Bio namespace

  package Bio::MyNewObject;
  use vars qw(@ISA);
  use strict;

  use Bio::Root::Root;
  @ISA = qw(Bio::Root::Root);

  sub new {
     my($class,@args) = @_;
     # call superclasses initialize
     my $self = $class->SUPER::new(@args);

     # do your own argument processing here

     my ($arg1) = $self->_rearrange([qw(NAMEDARGUMENT1)], @args);

     # set default attributes etc...

     return $self;
  }

=head2 Throwing Exceptions

Exceptions are die functions, in which the $@ variable (a scalar) is
used to indicate how it died. The exceptions can be caught using the
eval {} system. The Bioperl root object has a method called "-E<gt>throw"
which calls die but also provides a full stack trace of where this
throw happened on. So an exception like

  $obj->throw("I am throwing an exception");

Provides the following output on STDERR if is not caught.

------------- EXCEPTION: Bio::Root::Exception -------------
MSG: I am throwing an exception
STACK: Error::throw
STACK: Bio::Root::Root::throw /home/jason/bioperl/core/Bio/Root/Root.pm:313
-----------------------------------------------------------

indicating that this exception was thrown at line 7 of subroutine
my_subroutine, in myscript.pl

Exceptions can be caught using an eval block, such as

 my $obj = Bio::SomeObject->new();
 my $obj2
 eval {
   $obj2 = $obj->method1();
   $obj2->method2(10);
 }

 if( $@ ) {
   # exception was thrown
   &tell_user("Exception was thrown, preventing whatever I wanted to do. Actual exception $@");
   exit(0);
 }

 # else - use $obj2

Notice that the eval block can have multiple statements in it, and
also that if you want to use variables outside of the eval block, they
must be declared with my outside of the eval block (you are planning
to "use strict" in your scripts, aren't you!).

This context is particularly useful when objects are produced from a
database. This is because some exceptions are really due to problems
with the data in an object rather than the code. These sort of
exceptions are better tracked down when you know where the object came
from, not where in the code the exception is thrown.

One of the drawbacks to this scheme is that the attribute -E<gt>name is
"special" from Bioperl's perspective. We believe it is best to stay
away from using $obj-E<gt>name() to mean anything from the object's
perspective (for example -E<gt>id() ), leaving it free to be used as a
context for debugging purposes. You might prefer to overload the name
attribute to be "useful" for the object.

See scripts/root_object/error.pl for demonstration code.


=head1 Bioperl Interface design

Bioperl has been moving to a split between B<interface> and
B<implementation> definitions.  An interface is solely the definition
of what methods one can call on an object, without any knowledge of
how it is implemented. An implementation is an actual, working
implementation of an object. In languages like Java, interface
definition is part of the language. In Perl, like many aspects of Perl
you have to roll your own.

In Bioperl, the interface names are called Bio::MyObjectI, with the
trailing I indicating it is an interface definition of an object. The
interface files (sometimes nicknamed the 'I files') provide mainly
documentation on what the interface is, and how to use and implement
it. All the functions which the implementation is expected to provide
are defined as subroutines, and then die with an informative
warning. The exception to this rule are the implementation independent
functions (see L<Implementation functions in Interface files>).

Objects which want to implement this interface should inherit the
Bio::MyObjectI file in their @ISA array. This means that if the
implementation does not provide a method which the interface defines,
rather than the user getting a "method not found error" it gets a
"mymethod() was not defined in MyObjectI, but should have been" which
makes it clearer that whoever provided the implementation was to
blame, and not the caller/script writer.

When people want to check they have valid objects being passed to
their functions they should test the presence of the interface, not
the implementation. for example

  sub my_sequence_routine {
    my($seq,$other_argument) = @_;

    $seq->isa('Bio::SeqI') || die "[$seq] is not a sequence. Cannot process";

    # do stuff
  }

This is in contrast to

  sub my_incorrect_sequence_routine {
    my($seq,$other_argument) = @_;

    # this line is INCORRECT
    $seq->isa('Bio::Seq') || die "[$seq] is not a sequence. Cannot process";

    # do stuff
  }


=head2 Rationale of interface design

Some people might justifiably argue "why do this?". The main reason is
to support external objects from Bioperl, and allow them to masquerade
as real Bioperl objects. For example you might have your own quite
intricate sequence object which you want to use in Bioperl functions,
but don't want to lose your own neat coding. One option would be to
have a function which built a Bioperl sequence object from your
object, but then you would be endlessly building temporary objects and
destroying them, in particular if the script yo-yoed between your code
and Bioperl code.

A better solution would be to implement the Bio::SeqI interface. You
would read L<Bio::SeqI>, and then provide the methods
which it required, and put Bio::SeqI in your @ISA array. Then you
could pass in your object into Bioperl routines and I<eh voila> - you
B<are> a Bioperl sequence object.

(A problem might arise if your object has the same methods as the
Bio::SeqI methods but use them differently - your $obj-E<gt>id() might
mean provide the raw memory location of the object, whereas the
documentation for Bio::SeqI $obj-E<gt>id() says it should return the
human-readable name. If so you need to look into providing an
'Adaptor' class, as suggested in the Gang-of-four).

Interface classes really come into their own when we start leaving
Perl and enter extensions wrapped over C or over databases, or
through systems like CORBA to other languages, like Java/Python
etc. Here the "object" is often a very thin wrapper over the
a DBI interface, or an XS interface, and how it stores the object
is really different. By providing a very clear, implementation free
interface with good documentation there is a very clear target
to hit.

Some people might complain that we are doing something very "un-perl-like"
by providing these separate interface files. They are 90% documentation,
and could be provided anywhere, in many ways they could be merged with
the actual implementation classes and just made clear that if someone
wants to mimic a class they should override the following methods. However,
we (and in particular myself - Ewan) prefers a clear separation of the
interface. It gives us a much clearer way of defining what is going on.
It is in many ways just "design sugar" (as opposed to syntactic sugar)
to help us, but it really helps, so that's good enough justification to me.

=head2 Implementation functions in Interface files

One of the issues we discovered early on in using Interface files was
that there were methods that we would like to provide for classes
which were independent of their implementation. A good example is
a "Range" interface, which might define the following methods

  $obj->start()
  $obj->end()

Now a client to the object might want to use a $obj-E<gt>length() method.
because it is much easier than retrieving the two attributes and
substracting them. However, the -E<gt>length() method is just a pain for
someone providing the implementation to provide - once start() and
end() is defined, length is. There seems to be a catch-22 here: to
make an object definition good for a B<client> one needs to have
additional, helper methods "on top of" the interface, however to make
life easier for the B<object implementation> one wants to have the
bare minimum of functions defined which the implementer has to
provide.

In the Range interface this became more than annoyance, as a lot of the
"smarts" of the Range system was that we wanted to have the ability to
say

  if( $range->intersection($someother_range) )

We wanted a generic RangeI interface that we could apply to many
objects, with definitions required only for -E<gt>start, -E<gt>end and
-E<gt>strand. However we wanted the -E<gt>intersection, and
-E<gt>union methods to be on all ranges, without us having to
reimplement this every time.

Our (Matt Pocock and Ewan Birney's) solution was to allow
implementation into the RangeI interface file, but only when these
implementations sat "on top" of the interface definition and therefore
provided helper client operations. In a language like Java, we would
clearly have two classes, with a composition/delegation method:

   MyPublicSomethingClass has-a MyInternalSomethingInterface, with

   ADifferentImplemtation implements MyInternalSomethingInterface

However this is really heavy handed in Perl (and people were
complaining about having different implementation/interface
classes). We were quite happy about merging the implementation
independent functions with the interface definition, and we have used
this in other interfaces since then. The documentation has to be clear
about what is going on, but we think in general it is.

=head1 A Note on performance

Since Object Oriented programming in perl is not as elegant as
intentionally object oriented programming languages we incur
some overhead when calling the chained new constructors.  For most
cases this is perfectly okay as the object creation is not a
significant portion of many of the procedures.  However in certain
cases - reading in a large number of sequences with features requires
the creation of many objects and can perform poorly.  One can work
around this by creating the hashes directly and NOT chaining the new
calls.  An example of this is implemented in the Bio::SeqIO::FTHelper
objects in the treatment of Location objects for features. Please see
L<Bio::SeqIO::FTHelper> for details.

=head1 Getting Started

The Bioperl effort has managed to strike that balance between creating
reliable and powerful software and providing outlets for individual 
creativity and imagination. In this spirit of open-mindedness we invite
you to share your ideas with us at bioperl-l@bioperl.org.
!NO!GROK!
