#-----------------------------------------------------------------------------
# PACKAGE : Bio::Root::Object.pm
# AUTHOR  : Steve Chervitz (sac@bioperl.org)
# CREATED : 23 July 1996
# REVISION: $Id$
# STATUS  : Alpha
#
# For documentation, run this module through pod2html
# (preferably from Perl v5.004 or better).
#
# MODIFICATION NOTES: See bottom of file.
#
# Copyright (c) 1996-2000 Steve Chervitz. All Rights Reserved.
#           This module is free software; you can redistribute it and/or
#           modify it under the same terms as Perl itself.
#           Retain this notice and note any modifications made.
#-----------------------------------------------------------------------------

package Bio::Root::Object;
use strict;

require 5.002;
use Bio::Root::Global qw(:devel $AUTHORITY $CGI);

use Exporter ();

#use AutoLoader;
#*AUTOLOAD = \&AutoLoader::AUTOLOAD;

use vars qw(@EXPORT_OK %EXPORT_TAGS);
@EXPORT_OK = qw(&find_object &stack_trace &containment &_rearrange);
%EXPORT_TAGS = ( std => [qw(&stack_trace &containment)] );

use vars qw($ID %Objects_created);

use base qw(Bio::Root::Root);


# %Objects_created can be used for tracking all objects created.
# See _initialize() for details.

$ID       = 'Bio::Root::Object';

### POD Documentation:

=head1 NAME

Bio::Root::Object - A core Perl 5 object.

=head1 SYNOPSIS

  # Use this module as the root of your inheritance tree.

=head2 Object Creation

    require Bio::Root::Object;

    $dad = new Bio::Root::Object();
    $son = new Bio::Root::Object(-name    => 'Junior',
			         -parent  => $dad,
			         -make    => 'full');


See the L<new()|new> method for a complete description of parameters.
See also L<the USAGE section | USAGE>.

=head1 DESCRIPTION

L<Bio::Root::Object> attempts to encapsulate the "core" Perl5
object: What are the key data and behaviors ALL (or at least most) Perl5
objects should have?

=head2 Rationale

Use of L<Bio::Root::Object> within the Bioperl framework facilitates
operational consistency across the different modules defined within
the C<Bio::> namespace.  Not all objects need to derive from
L<Bio::Root::Object>. However, when generating lots of different types
of potentially complex objects which should all conform to a set of
basic expectations, this module may be handy.

At the very least, this module saves you from re-writing the L<new()|new>
method for each module you develop. It also permits consistent and
robust handling of C<-tag =E<gt> value> method arguments via the
L<Bio::Root::RootI::_rearrange()|Bio::Root::RootI> method and provides a
object-oriented way handle exceptions and warnings via the L<Bio::Root::Root::throw()|Bio::Root::Root> and L<Bio::Root::Root::warn()|Bio::Root::Root> methods.

See L<the APPENDIX section | APPENDIX> for some other handy methods.

=head2 Fault-Tolerant Objects

A major motivation for this module was to promote the creation of robust,
fault-tolerant Perl5 objects. The L<Bio::Root::Root::throw()|Bio::Root::Root> method relies on Perl's built-in
C<eval{}/die> exception mechanism to generate fatal exceptions.
The data comprising an exception is managed by the L<Bio::Root::Err>
module, which essentially allows the data thrown by a C<die()> event to be
wrapped into an object that can be easily examined and possibly re-thrown.

The intent here is three-fold:

=over 4

=item 1 Detailed error reporting.

Allow objects to report detailed information about the error condition
(who, what, where, why, how).

=item 2 Handle complex errors in objects.

The goal is to make it relatively painless to detect and handle the wide
variety of errors possible with a complex Perl object.
Perl's error handling mechanism is a might clunky when it comes to
handling complex errors within complex objects, but it is improving.

=item 3 Efficient & easy exception handling.

To enable robust exception handling without incurring a significant
performance penalty in the resulting code. Ideally, exception handling
code should be transparent to the cpu until and unless an exception
arises.

=back

These goals may at times be at odds and we are not claiming
to have achieved the perfect balance. Ultimately, we want self-
sufficient object-oriented systems able to deal with their own errors.
This area should improve as the module, and Perl, evolve.
One possible modification might be to utilize Graham Barr's L<Error>
module or Torsten Ekedahl's L<Experimental::Exception> module
(see L<Other Exception Modules>).
Technologies such as these may eventually be
incorporated into future releases of Perl. The exception handling
used by L<Bio::Root::Object> can be expected to change as Perl's
exception handling mechanism evolves.

B<TERMINOLOGY NOTE:> In this discussion and elsewhere in this module,
the terms "Exception" and "Error" are used interchangeably to mean
"something unexpected occurred" either as a result of incorrect user
input or faulty internal processing.

=head1 USAGE

=head2 Basic Exception handling

Object construction is a common place for exceptions to occur. By wrapping
the construction in an C<eval{ }> block, we can prevent the exception from
crashing the script and attempt to recover gracefully:

    # Package Foo.pm IS-A Bio::Root::Object.pm

    $obj = eval { new Foo(@data) };  # ending semicolon required.
    if($@) {
        print STDERR "\nTrouble creating Foo object: $@\n";
        recover_gracefully($@);
    }

A common strategy when generating lots of objects is to collect
data about which objects failed to build but still permit
the successfully created ones get processed:

    @errs = ();
    foreach $thing ( @stuff ) {
        my $obj = eval { new Foo($thing) };
        if($@) {
            push @err, [$thing, $@];
        }
        else {
            process_obj($obj);
        }
    }

Post-mortem reporting, logging, or analysis of the problems ensues:

    if(@errs) {
        printf "\n%d things failed:\n", scalar(@errs);
        foreach(@errs) { print "$err->[0], ";}

        print "\n\nTrapped exceptions:\n";
        foreach(@errs) { print "$err->[1]\n";}
    }

New with C<Perl 5.005> is the ability to C<die()> with an object
reference in C<$@> instead of just a string. This feature is not yet
exploited in Bio::Root::Object.pm but may be in future versions.
Bio::Root::Err.pm objects can be reconstructed from the contents of C<$@>:

    eval{ # exception-prone code here... };
    if($@) {
        $err = new Bio::Root::Err($@);
        printf "Trouble: %s\n". $err->msg;
        printf "Stack trace: %s\n". $err->stack;
    }


=head2 Demo Scripts

Some demo script that illustrate working with Bio::Root::Objects
are included with the distribution in the examples/root_object directory.


=head1 STRICTNESS & VERBOSITY

There are two global variables that can be used to control sensitivity to
exceptions/warnings and the amount of reporting for all objects within a process.
These are accessed via functions B<strictness()> and B<verbosity()> exported by
Bio::Root::Global (see L<Bio::Root::Global>).

  $STRICTNESS  - Regulates the sensitivity of the object to exceptions and warnings.

  $VERBOSITY   - Regulates the amount of reporting by an object.


The L<strict()|strict> and L<verbose()|verbose> methods of L<Bio::Root::Object>
originally operated at the the object level, to permit individual
strictness and verbosity levels for different objects. This level of
control is not usually required and can often be inconvenient; one
typically wants to set these properties globally for a given
script. While this sacrifices some flexibility, it saves time and
memory when working with lots of objects. For instance, child objects
don't have to worry about checking their parents to determine their
strictness/verbosity levels. Strictness and verbosity are
globally-defined values, but different classes of objects can be
differentially sensitive to these values depending on design criteria.

Strictness and verbosity can be positive or negative. Negative
verbosity equals terseness; negative strictness equals permissiveness.
In L<Bio::Root::Object> only the Bio::Root::Root::throw() and
Bio::Root::Root::warn() methods (see L<Bio::Root::Root>) are sensitive to
these values as indicated in the tables below:

    +---------+
    | throw() |         v e r b o s i t y
    +---------+ -------------------------------------
                   -1             0            1
    s           ----------   -----------   ----------
    t
    r   -2   --     throw() converted into warn()
    i
    c   -1   |   Exception    Exception      Exception
    t    0   |_  printed      printed        printed
    n    1   |   without      with           with stack
    e    2   |   stack trace  stack trace    trace and
    s        |                               sysbeep
    s


    +---------+
    | warn()  |         v e r b o s i t y
    +---------+ --------------------------------------
                   -1             0            1
    s           ----------   -----------   -----------
    t
    r   -2   |   Warning      Warning        Warning
    i   -1   |_  not          printed        printed
    c    0   |   printed      without        with stack
    t    1   |   but          stack trace    trace and
    n        |   attached*                   sysbeep
    e
    s    2   --      warn() converted into throw()
    s

     (*) Warnings will be attached to an object if the
     -record_err =>1 flag is set when constructing the object
     or if $object->record_err(1) is called subsequent to creation.

See the methods L<verbose()|verbose>, L<strict()|strict>, L<record_err()|record_err>,
Bio::Root::Root::throw(), and Bio::Root::Root::warn() in
L<Bio::Root::Root> for more details.


=head1 DEPENDENCIES

As the L<Bio::Root::Object> does not inherit from any modules
but wraps (i.e., provides an interface and delegates
functionality to) other modules in the Bio::Root:: hierarchy:

   Module                    Purpose
   --------------------      ------------------------------------
   Bio::Root::Err.pm         Exception handling
   Bio::Root::IOManager.pm   Input/output of object data or error data
   Bio::Root::Xref.pm        Arbitrary links between objects

All of these modules are loaded only when necessary.
L<Bio::Root::Err> is an object representing an exception.
L<Bio::Root::IOManager> and B<Bio::Root::Xref> are more experimental. They are
utilized via delegation, which permits them to be developed and utilized
independently of L<Bio::Root::Object>.

Since this module is at the root of potentially many different objects
in a particular application, efficiency is important. Bio::Root::Object.pm is
intended to be a lightweight, lean and mean module.


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other Bioperl modules.
Send your comments and suggestions preferably to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR

Steve Chervitz E<lt>sac@bioperl.orgE<gt>

See L<the FEEDBACK section | FEEDBACK> for where to send bug reports and comments.

=head1 VERSION

Bio::Root::Object.pm, 0.041


=head1 TODO

=over 2

=item * Experiment with other Exception classes.

Consider incorporating a more widely-used Error/Exception module
(see L<Other Exception Modules>).

=item * Think about integration with Data::Dumper.pm for persisting objects.

=back

=head1 SEE ALSO

L<Bio::Root::Err>       - Error/Exception object
L<Bio::Root::IOManager> - Input/Output manager object
L<Bio::Root::Vector>    - Manages dynamic lists of objects
L<Bio::Root::Xref>      - Cross-reference object
L<Bio::Root::Global>    - Manages global variables/constants

http://bio.perl.org/                       - Bioperl Project Homepage

=head2 Other Exception Modules

Experimental::Exception.pm   - ftp://ftp.matematik.su.se/pub/teke/
Error.pm                     - http://www.cpan.org/authors/id/GBARR/
Throwable.pm                 - mailto:kstevens@globeandmail.ca

http://genome-www.stanford.edu/perlOOP/exceptions.html

=head1 ACKNOWLEDGEMENTS

This module was developed under the auspices of the Saccharomyces Genome
Database:

http://genome-www.stanford.edu/Saccharomyces

Other Bioperl developers contributed ideas including Ewan Birney, Ian Korf,
Chris Dagdigian, Georg Fuellen, and Steven Brenner.

=head1 COPYRIGHT

Copyright (c) 1996-98 Steve Chervitz. All Rights Reserved.
This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

=cut


#
##
###
#### END of main POD documentation. '
###
##
#


=head1 APPENDIX

Methods beginning with a leading underscore are considered private
and are intended for internal use by this module. They are
B<not> considered part of the public interface and are described here
for documentation purposes only.

=cut

#
# This object is deprecated as the root of the inheritance tree, but some
# modules depend on it as a legacy. We issue a deprecation warning for all
# other modules.
#
my @inheriting_modules = ('Bio::Root::Object',
			  'Bio::Root::IOManager');


#######################################################
#               CONSTRUCTOR/DESTRUCTOR                #
#######################################################


=head2 new

 Purpose   : Creates a blessed object reference (hash) for the indicated class
           : and calls _initialize() for the class passing it all parameters.
 Usage     : new CLASS_NAME [ %named_parameters];
 Example   : $obj = new Bio::Root::Object 'george';
           : $obj = Bio::Root::Object->new(-name    => 'L56163',
	   : 			           -parent => $obj2 );
           : $obj = Bio::Root::Object->new();
 Returns   : Blessed hash reference.
 Argument  : Named parameters:  (PARAMETER TAGS CAN BE UPPER OR LOWERCASE).
           : (all are optional)
           :  -NAME       => arbitrary string to identify an object;
           :                 should be unique within its class.
           :  -PARENT    => blessed reference for an object that
           :                 is responsible for the present object
           :                 (e.g., a container).
           :  -MAKE       => string to specify special constructor option.
           :  -OBJ        => object reference for an object to be cloned.
           :  -RECORD_ERR => boolean (if true, attach all Err.pm objects generated by
	   :                 warn() or throw() calls to the present object;
	   :		     default = false).
           :
           : The use of STRICT and VERBOSE in constructors is no longer
           : necessary since there is no object-specific strict or verbose setting.
           : Use the strictness() and verbosity() functions exported by
           : Bio::Root::Global.pm. These options are still provided
           : in the constructor but the will affect *all* objects within a
           : given process.
           :
           :  -STRICT     => integer (level of strictness: -2, -1, 0, 1, 2).
           :  -VERBOSE    => integer (level of verbosity: -1, 0, 1)
           :                 Verbosity can be used to control how much reporting
           :                 an object should do generally. In this module,
           :                 verbosity affects the behavior of throw() and warn()
           :                 only.
           :
           :
 Comments  : This method creates blessed HASH references.
           : An object is free to define its own strict, and verbose
           : behavior as well as its own make (constructor) options.

See Also   : L<_initialize()|_initialize>, L<name()|name>, L<parent()|parent>, L<make()|make>, L<strict()|strict>, L<verbose()|verbose>, L<record_err()|record_err>, and Bio::Root::Root::throw() and Bio::Root::Root::warn() in L<Bio::Root::Root>

=cut

#----------
sub new {
#----------
    my($class, @param) = @_;
    my $self = {};
    bless $self, ref($class) || $class;
    $DEBUG==2 && print STDERR "CREATING $self";
    $self->_initialize(@param);
    $self;
}


=head2 _initialize

 Purpose   : Initializes key Bio::Root::Object.pm data (name, parent, make, strict).
           : Called by new().
 Usage     : n/a; automatically called by Bio::Root::Object::new()
 Returns   : String containing the -MAKE constructor option or 'default'
           : if none defined (if a -MAKE parameter is defined, the value
           : returned will be that obtained from the make() method.)
           : This return value saves any subclass from having to call
           : $self->make() during construction. For example, within a
           : subclass _initialize() method, invoke the Bio::Root::Object::
           : initialize() method as follows:
           :    my $make = $self->SUPER::_initialize(@param);
 Argument  : Named parameters passed from new()
           :  (PARAMETER TAGS CAN BE ALL UPPER OR ALL LOWER CASE).
 Comments  : This method calls name(), make(), parent(), strict(), index()
           : and thus enables polymorphism on these methods. To save on method
           : call overhead, these methods are called only if the data need
           : to be set.
           :
           : The _set_clone() method is called if the -MAKE option includes
           : the string 'clone' (e.g., -MAKE => 'clone').
           :
           : The index() method is called if the -MAKE option includes
           : the string 'index'. (This is an experimental feature)
           : (Example: -MAKE => 'full_index').
           :
           : NOTE ON USING _rearrange():
           :
           : _rearrange() is a handy method for working with tagged (named)
           : parameters and it permits case-insensitive in tag names
           : as well as handling tagged or un-tagged parameters.
           : _initialize() does not currently call _rearrange() since
           : there is a concern about performance when setting many objects.
           : One issue is that _rearrange() could be called with many elements
           : yet the caller is interested in only a few. Also, derived objects
           : typically invoke _rearrange() in their constructors as well.
           : This could particularly degrade performance when creating lots
           : of objects with extended inheritance hierarchies and lots of tagged
           : parameters which are passes along the inheritance hierarchy.
           :
           : One thing that may help is if _rearrange() deleted all parameters
           : it extracted. This would require passing a reference to the param list
           : and may add excessive dereferencing overhead.
           : It also would cause problems if the same parameters are used by
           : different methods or objects.

See Also   : L<new()|new>, L<make()|make>, L<name()|name>, L<parent()|parent>, L<strict()|strict>, L<index()|index>, L<verbose()|verbose>

=cut

#----------------
sub _initialize {
#----------------
    local($^W) = 0;
    my($self, %param) = @_;

    if(! grep { ref($self) =~ /$_/; } @inheriting_modules) {
	$self->warn("Class " . ref($self) .
		    " inherits from Bio::Root::Object, which is deprecated. ".
		    "Try changing your inheritance to Bio::Root::Root.");
    }
    my($name, $parent, $make, $strict, $verbose, $obj, $record_err) = (
	($param{-NAME}||$param{'-name'}), ($param{-PARENT}||$param{'-parent'}),
	($param{-MAKE}||$param{'-make'}), ($param{-STRICT}||$param{'-strict'}),
	($param{-VERBOSE}||$param{'-verbose'}),
        ($param{-OBJ}||$param{'-obj'}, $param{-RECORD_ERR}||$param{'-record_err'})
					  );

    ## See "Comments" above regarding use of _rearrange().
#	$self->_rearrange([qw(NAME PARENT MAKE STRICT VERBOSE OBJ)], %param);

    $DEBUG and do{ print STDERR ">>>> Initializing $ID (${\ref($self)}) ",$name||'anon';<STDIN>};

    if(defined($make) and $make =~ /clone/i) {
	$self->_set_clone($obj);

    } else {
	$name ||= ($#_ == 1 ? $_[1] : '');  # If a single arg is given, use as name.

	## Another performance issue: calling name(), parent(), strict(), make()
	## Any speed diff with conditionals to avoid method calls?

	$self->name($name) if $name;
	$self->parent($parent) if $parent;
	$self->{'_strict'}  = $strict  || undef;
	$self->{'_verbose'} = $verbose || undef;
	$self->{'_record_err'} = $record_err || undef;

	if($make) {
	    $make = $self->make($make);

	    # Index the Object in the global object hash only if requested.
	    # This feature is not used much. If desired, an object can always
	    # call Bio::Root::Object::index()  any time after construction.
	    $self->index() if $make =~ /index/;
	}
    }

    $DEBUG and print STDERR "---> Initialized $ID (${\ref($self)}) ",$name,"\n";

    ## Return data of potential use to subclass constructors.
#    return (($make || 'default'), $strict);   # maybe (?)
    return $make || 'default';
}



=head2 DESTROY

 Purpose   : Provides indication that the object is being reclaimed
           : by the GC for debugging purposes only.
 Usage     : n/a; automatically called by Perl when the ref count
           : on the object drops to zero.
 Argument  : n/a
 Comments  : Setting the global $DEBUG to 2 will print messages upon
           : object destruction.
           : Subclasses should override this method to
           : clean up any resources (open file handles, etc.)
           : The overridden method should end with a call to
           : SUPER::DESTROY;

See Also   : L<destroy()|destroy>

=cut

#-----------
sub DESTROY {
#-----------
    my $self=shift;

    $DEBUG==2 && print STDERR "DESTROY called in $ID for ${\$self->to_string} ($self)\n";
}


=head2 destroy

 Purpose   : Clean up any resources allocated by the object and
           : remove links to all objects connected to the present
           : object with the ultimate aim of signaling the GC to
           : reclaim all memory allocated for the object.
           : This method breaks links to any Err, IOManager, and Xref objects
           : and drops the present object as a child from any parent objects.
 Usage     : $object->destroy(); undef $object;
           : undef-ing the object reference signals the GC to reclaim
           : the object's memory.
 Returns   : undef
 Argument  : n/a
 Comments  : Circular reference structures are problematic for garbage
           : collection schemes such as Perl's which are based on reference
           : counting. If you create such structures outside of
           : the parent-child relationship, be sure to properly break
           : the circularity when destroying the object.
           : Subclasses should override this method to call destroy()
           : on any contained child objects. The overridden method
           : should end with a call to SUPER::destroy().
 Bugs      : Bio::Root::Xref.pm objects have not been tested and
           : may not be handled properly here.
           : Bio::Root::Vector.pm objects are also not yet handled
           : properly so beware of crunching lots of Vector objects.

=cut

#-------------'
sub destroy {
#-------------
## Note: Cannot delete parent and xref object refs since they are not
##       owned by this object, merely associated with it.
    my $self = shift;

    if(ref($self->{'_parent'})) {
	$self->{'_parent'}->_drop_child($self);
	undef $self->{'_parent'};
    }

    if(ref($self->{'_io'})) {
	$self->{'_io'}->destroy;
	undef $self->{'_io'};
    }

    if(ref($self->{'_err'})) {
	$self->{'_err'}->remove_all;
	undef $self->{'_err'};
    }

    if(ref($self->{'_xref'})) {
	$self->{'_xref'}->remove_all;
	undef $self->{'_xref'};
    }

    $self->_remove_from_index if scalar %Objects_created;
}


=head2 _drop_child

 Usage     : $object->_drop_child(object_ref)
           : Used internally by destroy().
 Purpose   : To remove a parent-to-child inter-object relationship.
           : The aim here is to break cyclical object refs to permit Perl's
           : GC to reclaim the object's memory. The expectation is that
           : a child object requests of its parent that the parent drop the
           : child object making the request. Parents do not drop children
           : unless requested by the child in question.
 Example   : $self->parent->_drop_child($self);
 Returns   : undef
 Argument  : Object reference for the child object to be dropped
 Throws    : Exception if an object ref is not provided as an argument.
 Comments  : This is a simplistic version that systematically checks every
           : data member, searching all top-level array, hash, and scalar
           : data members.
           : It does not recurse through all levels of complex data members.
           : Subclasses could override this method to handle complex child
           : data members for more optimal child searching. However, the
           : version here is probably sufficient for most situations.
           :
           : _drop_child() is called by Bio::Root::Object::destroy() for
           : all objects with parents.
 Status    : Experimental

See Also   : L<destroy()|destroy>

=cut

#---------------'
sub _drop_child {
#---------------
    my ($self, $child) = @_;
    my ($member, $found);

    $self->throw("Child not defined or not an object ($child).") unless ref $child;

    local($^W = 0);
    foreach $member (keys %{$self}) {
	next unless ref($self->{$member});
	# compare references.
	if (ref($self->{$member}) eq 'ARRAY') {
	    my ($i);
	    for($i=0; $i < @{$self->{$member}}; $i++) {
		if ($self->{$member}->[$i] eq $child) {
		    $DEBUG==2 && print STDERR "Removing array child $child\n";
		    undef $self->{$member}->[$i];
		    $found = 1; last;
		}
	    }
	} elsif(ref($self->{$member}) eq 'HASH') {
	    foreach(keys %{$self->{$member}}) {
		if ($self->{$member}->{$_} eq $child) {
		    $DEBUG==2 && print STDERR "Removing hash child $child\n";
		    undef $self->{$member}->{$_};
		    $found = 1; last;
		}
	    }
	} else {
	    if ($self->{$member} eq $child) {
		$DEBUG==2 && print STDERR "Removing child $child\n";
		undef $self->{$member};
		$found = 1; last;
	    }
	}
    }
    # Child not found:
    #   It is possible that a child object has a parent but has not yet been added to
    #   the parent due to a failure during construction of the child. Not warning.
    #$self->warn(sprintf "Child %s not found in Parent %s.", $child->to_string, $self->to_string) unless $found;

    undef;
}


#################################################################
#                    ACCESSORS & INSTANCE METHODS
#################################################################



=head2 name

 Usage     : $object->name([string]);
 Purpose   : Set/Get an object's common name.
 Example   : $myName = $myObj->name;
           : $myObj->name('fred');
 Returns   : String consisting of the object's name or
           : "anonymous <CLASSNAME>" if name is not set.
           : Thus, this method ALWAYS returns some string.
 Argument  : String to be used as the common name of the object.
           : Should be unique within its class.

See also   : L<has_name()|has_name>

=cut

#---------
sub name {
#---------
    my $self = shift;

#    $DEBUG and do{ print STDERR "\n$ID: name(@_) called.";<STDIN>; };

    if (@_) { $self->{'_name'} = shift }
    return defined $self->{'_name'} ? $self->{'_name'} : 'anonymous '.ref($self);
}


=head2 to_string

 Usage     : $object->to_string();
 Purpose   : Get an object as a simple string useful for debugging purposes.
 Example   : print $myObj->to_string;  # prints: Object <PACKAGE NAME> "<OBJECT NAME>"
 Returns   : String consisting of the package name + object's name
           : Object's name is obtained by calling the name() method.
 Argument  : n/a
 Throws    : n/a

See also   : L<name()|name>

=cut

#-------------
sub to_string {
#-------------
    my $self = shift;
    return sprintf "Object %s \"%s\"", ref($self), $self->name;
}


=head2 parent

 Usage     : $object->parent([object | 'null']);
 Purpose   : Set/Get the current object's source object.
           : An object's source object (parent) is defined as the object
           : that is responsible for creating the current object (child).
           : The parent object may also have a special mechanism for
           : destroying the child object. This should be included
           : in the parent object's DESTROY method which should end with a
           : call to $self->SUPER::DESTROY.
 Example   : $myObj->parent($otherObject);
 Returns   : Object reference for the parent object or undef if none is set.
 Argument  : Blessed object reference (optional) or the string 'null'.
           :  'null' = sets the object's _parent field to undef,
           :           breaking the child object's link to its parent.
 Throws    : Exception if argument is not an object reference or 'null'.
 Comments  : This method may be renamed 'parent' in the near future.
           : When and if this happens, parent() will still be supported but
           : will be deprecated.

See also   : L<destroy()|destroy>

=cut

#------------'
sub parent {
#------------
    my ($self) = shift;
    if (@_) {
	my $arg = shift;
	if(ref $arg) {
	    $self->{'_parent'} = $arg;
	} elsif($arg =~ /null/i) {
	    $self->{'_parent'} = undef;
	} else {
	    $self->throw("Can't set parent using $arg: Not an object");
	}
    }
    $self->{'_parent'};
}


=head2 src_obj

 Usage     : $object->src_obj([object | 'null']);
           : THIS METHOD IS NOW DEPRECATED. USE parent() INSTEAD.
 Purpose   : Set/Get the current object's source object (parent).

See also   : L<parent()|parent>

=cut

#------------'
sub src_obj {
#------------
    my ($self) = shift;
    $self->warn("DEPRECATED METHOD src_obj() CALLED. USE parent() INSTEAD.\n");
    $self->parent(@_);
}


=head2 has_name

 Usage     : $object->has_name();
 Purpose   : To determine if an object has a name.
 Returns   : True (1) if the object's {'Name'} data member is defined.
           : False otherwise.
 Comments  : One may argue, why not just use the name() method as a
           : combination setter/getter? has_name() is necessary for
           : the following reasons:
           :   (1) If an object's name is not defined, name() returns
           :        "anonymous <CLASSNAME>".
           :   (2) If an object's name is 0 (zero) or '' (empty string),
           : conditionals that simply check name() would fail incorrectly.

See also   : L<name()|name>

=cut

#--------------'
sub has_name { my $self = shift; return defined $self->{'_name'}; }
#--------------



=head2 make

 Usage     : $object->make([string]);
 Purpose   : Set/Get an object's constructor option.
           : make() is intended for use during object construction
           : to essentially permit alternate constructors since
           : Perl doesn't have a built-in mechanism for this.
 Example   : $make = $object->make();
           : $object->make('optionA');
 Returns   : String consisting of the object's make option
           : or 'default' if make is not set.
           : Thus, this method ALWAYS returns some string.
 Argument  : String to be used as an option during object construction.
 Comments  : A typical use of a make option is when cloning an object
           : from an existing object. In this case, the new() method
           : is called with -MAKE => 'clone'.

See also   : L<_initialize()|_initialize>, L<clone()|clone>

=cut

#----------'
sub make {
#----------
    my $self = shift;
    if(@_) { $self->{'_make'} = shift; }
    $self->{'_make'} || 'default';
}


=head2 err

 Usage     : $self->err([$data], [$delimit])
 Purpose   : Check for exceptions/warnings and get data about them.
           : (object validation and error data retrieval)
 Example   : $self->err && print "has err";
           : $errCount = $self->err('count');
           : $errMsgs  = $self->err('msg',"\t");
           : @errNotes = $self->err('note');
 Returns   : One of the following:
           :   1. If no arguments are given
           :          a. If the object has an error, the err data member is
           :             returned (this is an Bio::Root::Err.pm object),
           :          b. otherwise, undef is returned.
           :   2. The number of Errs in the object's err data member (if $data eq 'count').
           :   3. A string containing data from a specific field from an object's err member.
           :      -- If the object contains multiple errors, data for all errors will be
           :         strung together in reverse chronological order with each error's data
           :         preceeded by "Error #n\n" and followed by two delimiters.
           :   4. A list containing data from a specific field from an object's err member.
           :      -- If the object contains multiple errors, data for all errors will be
           :         added in reverse chronological order as separate elements in the list
           :         with NO "Error #n\n" identifier. Individual err list data
           :         (note,tech,stack) will be tab-delimited.
 Arguments : $data    = The name of a specific Err data member (see %Bio::Root::Err::ERR_FIELDS)
           :            OR 'count'.
           : $delimit = The delimiter separating a single Err's list data member's elements.
           :            Default is "\n". For multi-error objects, two of these
           :            delimiters separate data from different errors.
           :            If wantarray is true or delimiter is 'list', data from multiple
           :            errors will be returned as a list
           :
 Comments  : Since Err objects are now fatal and are not attached to the object by default,
           : this method is largely moot. It is a relic from the former
           : error "polling" days.
           : It is handy for accessing non-fatal warnings thrown by the object,
           : or in situations where fatal errors are converted to warnings
           : as when $self->strict is -1 or $WARN_ON_FATAL is true.
           : (Note: an object now only attaches Err objects to itself when
	   : constructed with -RECORD_ERR =>1 or if the global $RECORD_ERR is true).
           :
           : This method is intended mainly to test whether or not an object
           : has any Err objects associated with it and if so, obtaining the
           : Err object or specific data about it.
           : For obtaining ALL data about an error, use err_string().
           : For more detailed manipulations with the Err data, retrieve the
           : Err object and process its data as necessary.

See also   : L<err_string()|err_string>, L<print_err()|print_err>, L<Bio::Root::Err::get_all|Bio::Root::Err>

=cut

#----------
sub err {
#----------
    my( $self, $data, $delimit) = @_;

    return unless defined $self->{'_err'};

    $data    ||= 'member';
#    $delimit ||= (wantarray ? 'list' : "\n");
    $delimit ||= "\n";

    $data eq 'member' and return $self->{'_err'};
    $data eq 'count'  and return $self->{'_err'}->size();

    return $self->{'_err'}->get_all($data, $delimit );
}


=head2 record_err

 Usage     : $object->record_err([0|1]);
 Purpose   : Set/Get indicator for whether an object should save
           : the Bio::Root::Err.pm objects it generates via calls
           : to throw() or warn().
 Example   : $myObj->record_err(1)
 Returns   : Boolean (0|1)
 Argument  : Boolean (0|1)
 Comments  : Record_err is generally useful only for examining
           : warnings produced by an object, since calls to throw()
           : are normally fatal (unless strictness is set to -2).
           : To turn on recording of errors for all objects in a process,
           : use Bio::Root::Global::record_err().
 Status    : Experimental

See also   : L<err()|err>, and record_err() in L<Bio::Root::Err>

=cut

#---------------
sub record_err {
#---------------
    my $self = shift;

    if (@_) { $self->{'_record_err'} = shift }
    return $self->{'_record_err'} || 0;
}


=head2 err_state

 Usage     : $object->err_state();
 Purpose   : To assess the status of the object's Err object (if any).
 Returns   : A string: 'EXCEPTION' | 'WARNING' | 'FATAL' | 'OKAY'
           : (OKAY is returned if there are no Errors)
 Status    : Experimental

=cut

#-------------'
sub err_state {
#-------------
    my $self = shift;
    return 'OKAY' if not defined $self->{'_err'};
    $self->{'_errState'} || 'OKAY';
}


=head2 clear_err

 Purpose   : To remove any error associated with the given object.
 Usage     : $myObj->clear_err;

See  also  : L<err()|err>

=cut

#-------------
sub clear_err {
#-------------
    my $self = shift;
    undef $self->{'_err'};
}





=head2 containment

 Usage     : $aref = $object->containment();
           : Since this method can be exported, the following can be used:
           : $aref = containment($object);
 Purpose   : To determine the containment hierarchy of a object.
 Returns   : An array reference in which each element is a string
           : containing the class and name of
           : the object in which this object is contained.
           : Indentation increases progressively as the
           : hierarchy is traversed.
           : E.g.,  Object MyClass "Foo"
           :         Contained in object YourClass "Bar"
           :          Contained in object HisClass "Moo"
 Comments  : This method will report only one object at each level
           : since an object can currently have only one source object.
 Status    : Exported

See also   : L<err()|err>

=cut

#------------------
sub containment {
#------------------
    my( $self) = @_;
    my(@hierarchy);

#    print "$ID: getting err hierarchy.\n";
    push @hierarchy, $self->to_string;
    my $obj = $self;
    my $count = 0;

    while( ref $obj->parent) {
	$obj = $obj->parent;
	push @hierarchy, sprintf "%sContained in %s", ' ' x ++$count, $obj->to_string;
    }
    return \@hierarchy;
}


=head2 set_stats

 Usage     : $object->set_stats(KEY => DATA [,KEY2 => DATA2])
 Purpose   : To declare and initialize a set of statistics germain
           : to an object. Each statistic name becomes a data member
           : prefixed with an underscore (if not already) and first
           : character after the underscore is lowercased.
 Example   : $object->set_stats('num_A' =>1,
           :                    'Num_B' =>10 ):
           : This sets :
           :     $object->{'_num_A'} = 1
           :     $object->{'_num_B'} = 10;
 Returns   : n/a
 Comments  : This method implements a convention for naming Perl
           : object data members with a leading underscore,
           : consistent with the naming convention of private methods.
           : Data members should not be part of an object's public
           : interface. The leading underscore helps flag the members
           : as private and also prevents inadvertant clobbering.

=cut

#--------------'
sub set_stats {
#--------------
    my( $self, %param ) = @_;

    my ($val);
    foreach (keys %param) {
	$val = $param{$_};;
	s/^(\w)/_\l$1/;
	$self->{$_} = $val;
    }
}


=head2 strict

 Usage     : $object->strict( [-2|-1|0|1|2] );
           : warn $message if $object->strict > 0;
 Purpose   : To make the object hyper- or hyposensitive to exceptions & warnings.
           : Strict = 2  : extremely hyper-sensitive, converts warn() into throw().
           : Strict = 1  : hyper-sensitive, but calls to warn are not converted.
           : Strict = 0  : no change (throw() = fatal, warn() = non-fatal).
           : Strict = -1 : hypo-sensitive, but calls to throw are not converted.
           : Strict = -2 : extremely hypo-sensitive, converts throw() into warn()
           :
           : Two degrees of positive and negative values for strict permit
           : the following functionality:
           :   1. Setting strict to 2 or -2 leads to more dramatic strictness
           :      or permissiveness, respectively. With 2, all calls to warn()
           :      become calls to throw() and are therefore fatal. With -2,
           :      the opposite is true and calls to throw become non-fatal.
           :      A strict value of 2 is thus an object-level version of
           :      Perl's "use strict" pragma.
           :
           :   2. Setting strict to 1 or -1 does not affect the behavior of
           :      throw() and warn(). This allows an object to implement its
           :      its own strictness policy. A strict value of 1 is thus an
           :      an object-level version of Perl's -w flag.
           :
 Returns   : Integer between -2 to 2.
 Comments  : This method no longer accesses an object-specific strictness
           : level but rather the global $STRICTNESS variable
           : defined in Bio::Root::Global.pm and accessed via the
           : strictness() method exported by that package.
           : Thus, all objects share the same strictness which
           : is generally more convenient.
 Status    : Experimental

See also   : warn() and throw() in L<Bio::Root::Root>, L<STRICTNESS & VERBOSITY>, strictness() in L<Bio::Root::Global>

=cut

#------------
sub strict {
#------------
    my $self = shift;

    # Use global strictness?
    if( $self->{'_use_global_strictness'}) {
	return &strictness(@_);
    }
    else {
        # Object-specific strictness
        if (@_) { $self->{'_strict'} = shift; }
        defined($self->{'_strict'})
            ? return $self->{'_strict'}
            : (ref $self->{'_parent'} ? $self->{'_parent'}->strict : 0);
    }
}

=head2 use_global_strictness

 Usage     : $object->use_global_strictnness( [1|0] );
 Purpose   : Set/Get accessor for a flag indicating whether or not
           : to use the global strictness setting or to instead use
           : object-specific strictness.
 Returns   : Boolean
 Comments  :
 Status    : Experimental

See also   : L<strict()|strict>, L<STRICTNESS & VERBOSITY>, strictness() in L<Bio::Root::Global>

=cut

sub use_global_strictness {
    my ($self, $value) = @_;

    if( defined $value ) {
	$self->{'_use_global_strictness'} = $value;
    }

    return $self->{'_use_global_strictness'};
}


=head2 clone

 Purpose   : To deeply copy an object.
           : Creates a new object reference containing an exact
           : copy of an existing object and all its data members.
 Usage     : $myClone = $myObj->clone;
 Comments  : This method only clones the Bio::Root::Object data members.
           : To fully clone an object that has data members beyond
           : those inherited from Bio::Root::Object, you must provide a
           : constructor in your class to copy all data of an object
           : data into the clone. For an example, see how _set_clone()
           : is called by _initialize() in this class.
           :
           : clone() will pass the named parameters {-MAKE=>'clone'}
           : and {-OBJ=>$self} to the object's constructor. The
           : constructor should then either check the -MAKE parameter
           : directly or should check the return value from
           : a call to the superclass constructor (see _initialize()
           : for an example) and then copy the required data members from OBJ
           : into the new object, bypassing the normal construction process.
           : Cloning of objects has not been extensively tested.
           : USE WITH CAUTION.
 Status    : Experimental

See Also   : L<_set_clone()|_set_clone>, L<_initialize()|_initialize>

=cut

#-------------'
sub clone {
#-------------
    my($self) = shift;

#    warn sprintf "\nCloning %s \"%s\"\n\n", ref($self),$self->name;

    my $clone = $self->new(-MAKE    =>'clone',
			   -OBJ     =>$self);
    if($self->err()) { $clone->err($self->err); }
    $clone;
}



=head2 _set_clone

 Usage     : n/a; internal method used by _initialize()
           : $self->_set_clone($object_to_be_cloned)
 Purpose   : Deep copy all Bio::Root::Object.pm data members
           : into a new object reference.
           : (This is basically a copy constructor).
 Argument  : object ref for object to be cloned.
 Throws    : Exception if argument is not an object reference.
 Comments  : Data members which are objects are cloned (parent, io, err).
           : Cloning of objects has not been extensively tested.
           : USE WITH CAUTION.

See Also   : L<_initialize()|_initialize>

=cut

#----------------
sub _set_clone {
#----------------
    my($self, $obj) = @_;

    ref($obj) || throw($self, "Can't clone $ID object: Not an object ref ($obj)");

    local($^W) = 0;  # suppress 'uninitialized' warnings.

    $self->{'_name'}     = $obj->{'_name'};
    $self->{'_strict'}   = $obj->{'_strict'};
    $self->{'_make'}     = $obj->{'_make'};
    $self->{'_verbose'}  = $obj->{'_verbose'};
    $self->{'_errState'} = $obj->{'_errState'};
    ## Better to use can() with Perl 5.004.
    $self->{'_parent'}   = ref($obj->{'_parent'}) and $obj->{'_parent'}->clone;
    $self->{'_io'}       = ref($obj->{'_io'}) and $obj->{'_io'}->clone;
    $self->{'_err'}      = ref($obj->{'_err'}) and $obj->{'_err'}->clone;
}



=head2 verbose

 Usage     : $object->verbose([-1|0|1]);
 Purpose   : Set/Get an indicator for how much ruporting an object should do.
 Argument  : integer (-1, 0, or 1)
 Returns   : integer (-1, 0, or 1)
           : Returns 0 if verbosity has not been defined.
           : Verbosity > 0 indicates extra reporting.
           : Verbosity < 0 indicates minimal reporting.
           : Verbosity = 0 or undefined indicates default reporting.
 Comments  : This method no longer accesses an object-specific verbosity
           : level but rather the global $VERBOSITY variable
           : defined in Bio::Root::Global.pm and accessed via the
           : verbosity() method exported by that package.
           : Thus, all objects share the same verbosity which
           : is generally more convenient.
 Status    : Experimental

See Also   : L<strict()|strict>, L<STRICTNESS & VERBOSITY>, verbosity() in L<Bio::Root::Global>

=cut

#------------
sub verbose {
#------------
    my $self = shift;

    # Using global verbosity
    return &verbosity(@_);

    # Object-specific verbosity (not used unless above code is commented out)
    if(@_) { $self->{'_verbose'} = shift; }
    defined($self->{'_verbose'})
	? return $self->{'_verbose'}
	: (ref $self->{'_parent'} ? $self->{'_parent'}->verbose : 0);
}



=head1 I/O-RELATED METHODS (Delegated to L<Bio::Root::IOManager>)

=head2 _io

 Usage     : $object->_io()
 Purpose   : Get the Bio::Root::IOManager.pm object for the current object.

See also   : L<display()|display>, L<read()|read>, L<file()|file>

=cut

#-------
sub _io  {  my $self = shift;   return $self->{'_io'}; }
#-------



=head2 _set_io

 Usage     : n/a; internal use only.
 Purpose   : Sets a new Bio::Root::IOManager.pm object for the current object.

See also   : L<display()|display>, L<read()|read>, L<file()|file>

=cut

#------------
sub _set_io {
#------------
    my $self = shift;

    require Bio::Root::IOManager;

# See PR#192.
#    $self->{'_io'} = new Bio::Root::IOManager(-PARENT=>$self, @_);
    $self->{'_io'} = new Bio::Root::IOManager(-PARENT=>$self);
}



=head2 set_display

 Usage     : $object->set_display( %named_parameters).
           : See Bio::Root::IOManager::set_display() for a description of parameters.
 Purpose   : Sets the output stream for displaying data associated with an object.
           : Delegates to Bio::Root::IOManager::set_display().
 Argument  : Named parameters (optional).
           : See Bio::Root::IOManager::set_display() for a
           : description of arguments.
 Status    : Experimental
 Comments  : Sets the IOManager.pm object if it is not set.
           : I'm not satisfied with the current display()/set_display() strategy.

See also   : set_display() in L<Bio::Root::IOManager>

=cut

#----------------'
sub set_display {
#----------------
    my($self, @param) = @_;

    $self->_set_io(@param) if !ref($self->{'_io'});

    eval { $self->{'_io'}->set_display(@param);  };

    if($@) {
	my $er = $@;
	$self->throw(-MSG=>$er, -NOTE=>"Can't set_display for ${\$self->name}");
    }

   return $self->{'_io'}->fh;
}


=head2 display

 Usage     : $object->display( named parameters)
           : See Bio::Root::IOManager::display() for a description of parameters.
 Purpose   : Output information about an object's data.
           : Delegates this task to Bio::Root::IOManager::display()
 Argument  : Named parameters for IOManager::set_display()
 Status    : Experimental
 Comments  : Sets the IOManager.pm object if it is not set.
           : IOManager::set_display()is then called on the new IOManager object.
           :
           : The motivation behind the display() method and IOManager.pm
           : is to allow for flexible control over output of an
           : object's data to/from filehandles, pipes, or STDIN/STDOUT,
           : and for passing file handles between objects. Currently,
           : it is used mainly for output to STDOUT.
           :
           : There is some concern whether this much functionality is
           : actually necessary, hence the "Experimental" status of this
           : method.
           :
           : -------
           : It might be worthwhile to also have a string() method
           : that will put an object's data into a string that can be
           : further processed as desired. Stringification for persistence
           : issues might be best handled by Data::Dumper.pm.
           :
           : When overriding this method, use the following syntax:
           :
           : sub display {
           :    my ($self, %param) = @_;
           :    $self->SUPER::display(%param);
           :    my $OUT = $self->fh();
           :    print $OUT "\nSome data...\n";
           :  ...
           : }
           : Now $OUT holds a filhandle reference (or the string 'STDOUT')
           : which can be passed to other methods to display different
           : data for the object.
           : _set_display() is automatically called with $OUT as the sole
           : argument (after $self) by IOManager.pm::display()
           : if the -SHOW parameter is set to 'stats' or 'default'.
           :
 Bugs      : Because the $OUT variable can be a FileHandle or a string,
           : it is necessary to include the line before using $OUT in
           : print statements:
           : I am considering a cleaner way of dealing with this.
           : Setting $OUT to a glob (*main::STDOUT) was unsuccessful.
           :
           : I'm not satisfied with the current display()/set_display() strategy.

See also   : display() in L<Bio::Root::IOManager>

=cut

#-------------
sub display {
#-------------
    my( $self, @param ) = @_;
    $self->{'_io'} || $self->set_display(@param);
    $self->{'_io'}->display(@param);
}




=head2 _display_stats

 Usage     : n/a; called automatically by Bio::Root::Object::display(-SHOW=>'stats');
 Purpose   : Display stereotypical data for an object.
           : Automatically called via display().
 Argument  : Filehandle reference or string 'STDOUT' 'STDIN' 'STDERR'
 Status    : Experimental

See also   : L<display()|display>

=cut

#-------------------
sub _display_stats {
#-------------------
    my($self, $OUT) = @_;


    printf ( $OUT "%-15s: %s\n","NAME", $self->name());
    printf ( $OUT "%-15s: %s\n","MAKE", $self->make());
    if($self->parent) {
	printf ( $OUT "%-15s: %s\n","PARENT", $self->parent->to_string);
    }
    printf ( $OUT "%-15s: %d\n",'ERRORS', (defined $self->err('count') ? $self->err('count') : 0)); ###JES###
    printf ( $OUT "%-15s: %s\n","ERR STATE", $self->err_state());
    if($self->err()) {
	print $OUT "ERROR:\n";
	$self->print_err();
    }
}



=head2 read

 Usage     : $object->read( named parameters)
           : See Bio::Root::IOManager::read() for a description of parameters.
 Purpose   : Inputs data from an arbitrary source (file or STDIN).
           : Delegates this task to Bio::Root::IOManager::read().
 Status    : Experimental
 Comments  : Sets the IOManager.pm object if it is not set.
           : See the comments for the display() method for some comments
           : about IO issues for objects.
           : Note that the read() method uses a different strategy than
           : the display() method.
           : IO issues are considered experimental.

See also   : L<display()|display>, read() in L<Bio::Root::IOManager>

=cut

#--------
sub read {
#--------
    my $self = shift;

    $self->_set_io(@_) if not defined $self->{'_io'};

    $self->{'_io'}->read(@_);
}



=head2 fh

 Usage     : $object->fh(['name'])
           : See Bio::Root::IOManager::fh() for a complete usage description.
 Purpose   : Get an object's current FileHandle object or IO stream indicator.
           : Delegates to Bio::Root::IOManager.pm.
 Status    : Experimental
 Comments  : Sets the IOManager.pm object if it is not set.

See also   : fh() in L<Bio::Root::IOManager>

=cut

#--------'
sub fh      {
#--------
    my $self = shift;
    $self->_set_io(@_) if !defined $self->{'_io'};
    $self->{'_io'}->fh(@_);
}


=head2 show

 Usage     : $object->show()
           : See Bio::Root::IOManager::show() for details.
 Purpose   : Get the string used to specify what to display
           : using the display() method.
           : Delegates to Bio::Root::IOManager.pm.
 Status    : Experimental
 Comments  : Sets the IOManager.pm object if it is not set.

See also   : show() in L<Bio::Root::IOManager>, set_display() in L<Bio::Root::IOManager>

=cut

#-----------
sub show    {
#-----------
    my $self = shift;
    $self->_set_io(@_) if !defined $self->{'_io'};
    $self->{'_io'}->show;
}



=head2 file

 Usage     : $object->file()
           : See Bio::Root::IOManager::file() for details.
 Purpose   : Set/Get name of a file associated with an object.
           : Delegates to Bio::Root::IOManager.pm.
 Status    : Experimental
 Comments  : Sets the IOManager.pm object if it is not set.

See also   : file() in L<Bio::Root::IOManager>

=cut

#---------
sub file    {
#---------
    my $self = shift;
    $self->_set_io(@_) if !defined $self->{'_io'};
    $self->{'_io'}->file(@_);
}


=head2 compress_file

 Usage     : $object->compress_file([filename])
           : See Bio::Root::IOManager::compress_file() for details.
 Purpose   : Compress a file associated with the current object.
           : Delegates to Bio::Root::IOManager.pm.
 Throws    : Propagates exceptions thrown by Bio::Root::IOManager.pm
 Status    : Experimental
 Comments  : Sets the IOManager.pm object if it is not set.

See also   : L<file()|file>, compress_file() in L<Bio::Root::IOManager>

=cut

#-------------------
sub compress_file {
#-------------------
    my $self = shift;
    $self->_set_io(@_) if !defined $self->{'_io'};
    $self->{'_io'}->compress_file(@_);
}



=head2 uncompress_file

 Usage     : $object->uncompress_file([filename])
           : Delegates to Bio::Root::IOManager.pm.
 Purpose   : Uncompress a file associated with the current object.
 Throws    : Propagates exceptions thrown by Bio::Root::IOManager.pm
 Status    : Experimental
 Comments  : Sets the IOManager.pm object if it is not set.

See also   : L<file()|file>, uncompress_file() in L<Bio::Root::IOManager>

=cut

#--------------------
sub uncompress_file {
#--------------------
    my $self = shift;
    $self->_set_io(@_) if !defined $self->{'_io'};
    $self->{'_io'}->uncompress_file(@_);
}


=head2 delete_file

 Usage     : $object->delete_file([filename])
           : See Bio::Root::IOManager::delete_file() for details.
 Purpose   : Delete a file associated with the current object.
           : Delegates to Bio::Root::IOManager.pm.
 Throws    : Propagates exceptions thrown by Bio::Root::IOManager.pm
 Status    : Experimental
 Comments  : Sets the IOManager.pm object if it is not set.

See also   : L<file()|file>, delete_file() in L<Bio::Root::IOManager>

=cut

#-----------------
sub delete_file {
#-----------------
    my $self = shift;
    $self->_set_io(@_) if !defined $self->{'_io'};
    $self->{'_io'}->delete_file(@_);
}


=head2 file_date

 Usage     : $object->file_date( %named_parameters )
           : See Bio::Root::IOManager::file_date() for details.
 Purpose   : Obtain the last modified data of a file.
           : Delegates to Bio::Root::IOManager.pm.
 Example   : $object->file_date('/usr/home/me/data.txt');
 Throws    : Propagates exceptions thrown by Bio::Root::IOManager.pm
 Status    : Experimental
 Comments  : Sets the IOManager.pm object if it is not set.

See also   : L<file()|file>, file_date() in L<Bio::Root::IOManager>

=cut

#---------------
sub file_date {
#---------------
    my $self = shift;
    $self->_set_io(@_) if !defined $self->{'_io'};
    $self->{'_io'}->file_date(@_);
}



=head1 EXPERIMENTAL METHODS


=head2 xref

 Usage     : $object->xref([object | 'null']);
 Purpose   : Sets/Gets an object(s) cross-referenced
           : to the current object.
 Example   : $myObj->xref('null');       #remove all xrefs
           : $myObj->xref($otherObject); #add a cross referenced object
 Argument  : Object reference or 'null' ('undef' also accepted).
 Returns   : Object reference or undef if the object has no xref set.
 Throws    : fatal error if argument is not an object reference or 'null'.
 Comments  : An Xref.pm object is a vectorized wrapper for an object.
           : Thus, the number of objects cross-referenced can grow
           : and shrink at will.
 Status    : Experimental
 WARNING   : NOT FULLY TESTED.

See Also   : L<Bio::Root::Xref>

=cut

#---------
sub xref  {
#---------
    my $self = shift;
    if(@_) {
	my $arg = shift;
	if(ref $arg) {
	    require Bio::Root::Xref;

	    if( !defined $self->{'_xref'}) {
		$self->{'_xref'} = new Bio::Root::Xref(-PARENT =>$self,
						       -OBJ     =>$arg);
	    } else {
		$self->{'_xref'}->add($arg);
	    }
	} elsif($arg =~ /null|undef/i) {
	    undef $self->{'_xref'};
	} else {
	    $self->throw("Can't set Xref using $arg: Not an object");
	}
    }

    $self->{'_xref'};
}



=head2 index

 Purpose   : To add an object to a package global hash of objects
           : for tracking or rapid retrieval.
 Usage     : $self->index();
 Status    : Experimental
 Comments  : The object's name is used to index it into a hash. Objects in
           : different classes (packages) will be indexed in different hashes.
           : An object's name should thus be unique within its class.
           : To find an object, use find_object().
           : Uses the package global %Objects_created.

See also   : L<find_object()|find_object>

=cut

#----------
sub index {
#----------
    my $self    = shift;
    my $class   = ref $self;
    my $objName = $self->{'_name'};

    if( not defined $objName ) {
	$self->throw("Can't index $class object \"$objName\".");
    }

    $DEBUG and do{ print STDERR "$ID: Indexing $class object \"$objName\"."; <STDIN>; };

    $Objects_created{ $class }->{ $objName } = $self;
}

#----------------------
sub _remove_from_index {
#----------------------
    my $self    = shift;
    my $class   = ref $self;
    my $objName = $self->{'_name'};

    undef $Objects_created{$class}->{$objName} if exists $Objects_created{$class}->{$objName};
}



=head2 find_object

 Purpose   : To obtain any object reference based on its unique name
           : within its class.
 Usage     : $myObj = &find_object('fred');
           : No need to specify the class (package) name of the object.
 Comments  : To use this method, the object must be previously
           : indexed by Bio::Root::Object.pm. This can be accomplished
           : by including 'index' in the -MAKE parameter during object
           : construction OR by calling the index() method on the
           : the object at any point after construction.
           : This is not an instance method.
 Status    : Experimental

See  also  : L<index()|index>

=cut

#---------------
sub find_object {
#---------------
    my $name   = shift;   # Assumes name has been validated.
    my $class  = undef;
    my $object = undef;

    foreach $class ( keys %Objects_created ) {
	if( exists $Objects_created{ $class }->{ $name } ) {
	    $object = $Objects_created{ $class }->{ $name };
	    last;
	}
    }
    $object;
}



=head2 has_warning

 Purpose   : Test whether or not an object has a non-fatal error (warning).
 Usage     : $self->has_warning;
 Comments  : This method is not usually needed. Checking err() is
           : sufficient since throw()ing an exception is a fatal event
           : and must be handled when it occurs.
 Status    : Experimental

See also   : L<err()|err>, warn() in L<Bio::Root::Root>, throw() in L<Bio::Root::Root>

=cut

#----------------
sub has_warning {
#----------------
    my $self = shift;
    my $errData = $self->err('type');
    return 1 if $errData =~ /WARNING/;
    0;
}



=head2 print_err

 Usage     : print_err([-WHERE=>FileHandle_object [,-SHOW=>msg|note|tech|stack] or any combo])
 Purpose   : Reports error data for any errors an object may have
           : as a string. This will only print warnings since exceptions
           : are fatal (unless a strictness of -2 is used).
 Example   : $myObj->print_err;
           : $myObj->print_err(-WHERE=>$myObj->fh('err'), -SHOW=>'msgtechstack');
 Argument  : SHOW parameter  : specify a sub-set of the err data.
           : WHERE parameter : specify a filehandle for printing.
 Returns   : n/a
 Status    : Experimental

See also   : L<err_string()|err_string>, L<strict()|strict>

=cut

#-------------
sub print_err {
#-------------
    my( $self, %param ) = @_;

#    print "$ID: print_err()\n";

    my $OUT = $self->set_display(%param);

#    print "$ID: OUT = $OUT\n";

    print $OUT $self->err_string( %param );

#    print "$ID: done print_err()\n";
}



=head2 err_string

 Usage     : err_string([-SHOW =>msg|note|tech|stack])
           : err_string([-SHOW =>'msgnote'] or other combos)
 Purpose   : Reports all warnings generated by the object as a string.
 Example   : $errData = $myObj->err_string;
           : print MYHANDLE $myObj->err_string();
 Argument  : SHOW parameter : return a specific sub-set of the err data.
 Returns   : A string containing the error data of the object.
 Comments  : This method is provided as a safer and slightly easier to type
           : alternative to $self->err->string.
 Status    : Experimental

See also   : L<print_err()|print_err>, string() in L<Bio::Root::Err>

=cut

#----------------
sub err_string {
#----------------
    my( $self, %param ) = @_;
    my($out);
    my $errCount = $self->err('count');

#    print "$ID: err_string(): count = $errCount\n";

    if( $errCount) {
	$out = sprintf("\n%d error%s in %s \"%s\"\n",
		       $errCount, $errCount>1?'s':'', ref($self), $self->name);
	$out .= $self->err->string( %param );
    } else {
	$out = sprintf("\nNo errors in %s \"%s\"\n", ref($self), $self->name);
    }
    $out;
}




#################################################################
#            DEPRECATED or HIGHLY EXPERIMENTAL METHODS
#################################################################

=head1 HIGHLY EXPERIMENTAL/DEPRECATED METHODS

=head2 terse

 Usage     : $object->terse([0|1]);
 Purpose   : Set/Get an indicator to report less than the normal amount.
 Argument  : Boolean (0|1)
 Returns   : Boolean (0|1)
 Comments  : This method is for reducing the amount of reporting
           : an object will do.
           : terse can be set during object construction with the
           : -TERSE => 1 flag.
           : Not putting this method in IOManager.pm since that class
           : is concerned with "where" to report, not "what" or "how much".
 Status    : Deprecated
           : Use verbose() with a negative value instead.

See also   : L<verbose()|verbose>

=cut

#----------
sub terse {
#----------
    my $self = shift;
    if(@_) { $self->{'_verbose'} = -1 * shift; }

    $self->warn("Deprecated method 'terse()'. Use verbose(-1) instead.");

    my $verbosity = $self->{'_verbose'} or
	(ref $self->{'_parent'} and $self->{'_parent'}->verbose) or 0;

    return $verbosity * -1;
}


#----------------------
=head2 set_err_data()
#----------------------

 Usage     : $object->set_err_data( field, data);
 Purpose   : Alters data within the last error set by the object.
           : Interface to Bio::Root::Err::set().
 Returns   : Calls Bio::Root::Err::set()
 Argument  : field = string, name of Bio::Root::Err.pm data field to set.
           : data  = string, data to set it to.
 Throws    : Exception if object has no errors.
 Status    : Deprecated

See Also   : set() in L<Bio::Root::Err>

=cut

#-----------------
sub set_err_data {
#-----------------
    my( $self, $field, $data) = @_;

    $self->throw("Object has no errors.") if !$self->{'_err'};

#    print "$ID: set_err_data($field)  with data = $data\n  in object ${\$self->name}:\n", $self->err->last->string(-CURRENT=>1); <STDIN>;

    $self->{'_err'}->last->set( $field, $data );
}

=head2 set_read

 Usage     : see Bio::Root::IOManager::set_read()
 Purpose   : Sets an input stream for importing data associated with an object.
           : Delegates to Bio::Root::IOManager::set_read().
 Status    : Experimental
 WARNING!  : This method has not been tested.

See also   : set_read() in L<Bio::Root::IOManager>

=cut

#--------------
sub set_read {
#--------------
    my($self,%param) = @_;

    $self->_set_io(%param) if !defined $self->{'_io'};

    $self->{'_io'}->set_read(%param);
}



=head2 set_log_err

 Usage     : see Bio::Root::IOManager::set_log_err()
 Purpose   : Sets the output stream for logging information about
           : an object's errors.
           : Delegates to Bio::Root::IOManager::set_log_err().
 Status    : Experimental
 WARNING!  : This method has not been tested.

See also   : set_log_err() in L<Bio::Root::IOManager>

=cut

#---------------'
sub set_log_err {
#---------------
    my($self,%param) = @_;

    $self->_set_io(%param) if !defined $self->{'_io'};

    $self->{'_io'}->set_log_err(%param);
}


1;
__END__


#####################################################################################
#                                  END OF CLASS                                     #
#####################################################################################

=head1 FOR DEVELOPERS ONLY

=head2 Data Members

Information about the various data members of this module is provided for those
wishing to modify or understand the code. Two things to bear in mind:

=over 4

=item 1 Do NOT rely on these in any code outside of this module.

All data members are prefixed with an underscore to signify that they are private.
Always use accessor methods. If the accessor doesn't exist or is inadequate,
create or modify an accessor (and let me know, too!).

=item 2 This documentation may be incomplete and out of date.

It is easy for this documentation to become obsolete as this module is still evolving.
Always double check this info and search for members not described here.

=back

An instance of Bio::Root::Object.pm is a blessed reference to a hash containing
all or some of the following fields:

 FIELD          VALUE
 ------------------------------------------------------------------------
  _name         Common name for an object useful for indexing.
 	        Should be unique within its class.

  _parent       The object which created and is responsible for this object.
                When a parent is destroyed, it takes all of its children with it.

  _err          Bio::Root::Err.pm object reference. Undefined if the object has no error
                or if the _record_err member is false (which is the default).
 	        If object has multiple errors, err becomes a linked
 	        list of Err objects and the err member always points to latest err.
 	        In theory, an object should care only about whether or not it HAS
 	        an Err not how many it has. I've tried to make the management of
 	        multiple errors as opaque as possible to Bio::Root::Object.

 _errState      One of @Bio::Root::Err::ERR_TYPES. Allows an object to quickly determine the
 	        the type of error it has (if any) without having to examine
 	        potentially multiple Err object(s).

  _xref         Bio::Root::Xref object (Vector) for tracking other object(s) related to the
 	        present object not by inheritance or composition but by some arbitrary
 	        criteria. This is a new, experimental feature and is not fully implemented.

  _make         Used as a switch for custom object initialization. Provides a
 	        mechanism for alternate constructors. This is somewhat experimental.
 	        It may be useful for contruction of complex objects and may be of
 	        use for determining how an object was constructed post facto.

  _io           Bio::Root::IOManager.pm object reference. Used primarily for handling the
	        display of an object's data.

  _strict       Integer flag to set the sensitivity to exceptions/warnings
 	        for a given object.

  _verbose      Boolean indicator for reporting more or less than the normal amount.

  _record_err   Boolean indicator for attaching all thrown exception objects
                to the current object. Default = false (don't attach exceptions).

=cut


MODIFICATION NOTES:
-----------------------
0.041, sac --- Thu Feb  4 03:50:58 1999
 * warn() utilizes the Global $CGI indicator to supress output
   when script is running as a CGI.

0.04, sac --- Tue Dec  1 04:32:01 1998
 *  Incorporated the new globals $STRICTNESS and $VERBOSITY
    and eliminated WARN_ON_FATAL, FATAL_ON_WARN and DONT_WARN.
 *  Deprecated terse() since it is better to think of terseness
    as negative verbosity.
 *  Removed autoloading-related code and comments.

0.035, 28 Sep 1998, sac:
  * Added _drop_child() method to attempt to break cyclical refs
    between parent and child objects.
  * Added to_string() method.
  * Err objects no longer know their parents (no need).

0.031, 2 Sep 1998, sac:
  * Documentation changes only. Wrapped the data member docs
    at the bottom in POD comments which fixes compilation bug
    caused by commenting out __END__.

0.03, 16 Aug 1998, sac:
  * Calls to warn() or throw() now no longer result in Err.pm objects
    being attached to the current object. For discussion about this
    descision, see comments under err().
  * Added the -RECORD_ERR constructor option and Global::record_err()
    method to enable the attachment of Err.pm object to the current
    object.
  * Minor bug fixes with parameter handling (%param -> @param).
  * Added note about AUTOLOADing.

0.023, 20 Jul 1998, sac:
  * Changes in Bio::Root::IOManager::read().
  * Improved memory management (destroy(), DESTROY(), and changes
    in Bio::Root::Vector.pm).

0.022, 16 Jun 1998, sac:
  * Changes in Bio::Root::IOManager::read().

0.021, May 1998, sac:
  * Touched up _set_clone().
  * Refined documentation in this and other Bio::Root modules
    (converted to use pod2html in Perl 5.004)


