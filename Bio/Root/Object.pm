#-----------------------------------------------------------------------------
# PACKAGE : Bio::Root::Object.pm
# AUTHOR  : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED : 23 July 1996
# REVISION: $Id$
# STATUS  : Alpha
#            
# For documentation, run this module through pod2html 
# (preferably from Perl v5.004 or better).
#
# MODIFICATION NOTES: See bottom of file.
#
# Copyright (c) 1996-98 Steve A. Chervitz. All Rights Reserved.
#           This module is free software; you can redistribute it and/or 
#           modify it under the same terms as Perl itself.
#           Retain this notice and note any modifications made.
#-----------------------------------------------------------------------------

package Bio::Root::Object;

require 5.002;
use Bio::Root::Global qw(:devel $AUTHORITY $CGI);
use Exporter ();
#use AutoLoader; 
#*AUTOLOAD = \&AutoLoader::AUTOLOAD;

@EXPORT_OK = qw($VERSION &find_object &stack_trace &containment &_rearrange);  
%EXPORT_TAGS = ( std => [qw(&stack_trace &containment)] );

use strict;
use vars qw($ID $VERSION %Objects_created $Revision);

# %Objects_created can be used for tracking all objects created.
# See _initialize() for details.

$ID       = 'Bio::Root::Object';
$VERSION  = 0.041;
$Revision = '$Id$';  #'

### POD Documentation:

=head1 NAME

Bio::Root::Object - A core Perl 5 object.

=head1 SYNOPSIS

=head2 Object Creation

    require Bio::Root::Object;
 
    $dad = new Bio::Root::Object();
    $son = new Bio::Root::Object(-name    => 'Junior', 
			         -parent  => $dad,
			         -make    => 'full');


See the L<new>() method for a complete description of parameters.
See also L<USAGE>.


=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.

=head1 DESCRIPTION

B<Bio::Root::Object.pm> attempts to encapsulate the "core" Perl5
object: What are the key data and behaviors ALL (or at least most) Perl5
objects should have?

=head2 Rationale

Use of B<Bio::Root::Object.pm> within the Bioperl framework facilitates
operational consistency across the different modules defined within
the B<Bio::> namespace.  Not all objects need to derive from
B<Bio::Root::Object.pm>. However, when generating lots of different types
of potentially complex objects which should all conform to a set of
basic expectations, this module may be handy.

At the very least, this module saves you from re-writing the L<new>()
method for each module you develop. It also permits consistent and
robust handling of C<-tag =E<gt> value> method arguments via the
L<_rearrange>() method and provides a object-oriented way handle
exceptions and warnings via the L<throw>() and L<warn>() methods. 
See the L<APPENDIX> for some other handy methods.

=head2 Fault-Tolerant Objects

A major motivation for this module was to promote the creation of robust,
fault-tolerant Perl5 objects. The L<throw>() method relies on Perl's built-in 
C<eval{}/die> exception mechanism to generate fatal exceptions.
The data comprising an exception is managed by the B<Bio::Root::Err.pm> 
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

These goals may at times be at odds and I am not claiming 
to have achieved the perfect balance. Ultimately, we want self-
sufficient object-oriented systems able to deal with their own errors. 
This area should improve as the module, and Perl, evolve. 
One possible modification might be to utilize Graham Barr's B<Error.pm>
module or Torsten Ekedahl's B<Experimental::Exception.pm> module 
(see L<Other Exception Modules>).
Technologies such as these may eventually be
incorporated into future releases of Perl. The exception handling 
used by B<Bio::Root::Object.pm> can be expected to change as Perl's 
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

New with B<Perl 5.005> is the ability to C<die()> with an object 
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

Some demo script that illustrate working with Bio::Root::Objects are included with the distribution (see L<INSTALLATION>). These are also available at:

    http://bio.perl.org/Core/Examples/Root_object


=head1 STRICTNESS & VERBOSITY

There are two global variables that can be used to control sensitivity to
exceptions/warnings and the amount of reporting for all objects within a process.
These are accessed via functions C<strictness()> and C<verbosity()> exported by
B<Bio::Root::Global.pm>.

  $STRICTNESS  - Regulates the sensitivity of the object to exceptions and warnings.

  $VERBOSITY   - Regulates the amount of reporting by an object.


The L<strict>() and L<verbose>() methods of B<Bio::Root::Object.pm>
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
In B<Bio::Root::Object.pm> only the L<throw>() and L<warn>() methods
are sensitive to these values as indicated in the tables below:

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

See the methods L<verbose>(), L<strict>(), L<throw>(), L<warn>(), and
L<record_err>() for more details.


=head1 DEPENDENCIES

As the B<Bio::Root::Object.pm> does not inherit from any modules 
but wraps (i.e., provides an interface and delegates 
functionality to) other modules in the Bio::Root:: hierarchy:

   Module                    Purpose
   --------------------      ------------------------------------
   Bio::Root::Err.pm         Exception handling
   Bio::Root::IOManager.pm   Input/output of object data or error data
   Bio::Root::Xref.pm        Arbitrary links between objects

All of these modules are loaded only when necessary.
B<Bio::Root::Err.pm> is an object representing an exception.
B<Bio::Root::IOManager.pm> and B<Bio::Root::Xref.pm> are more experimental. They are
utilized via delegation, which permits them to be developed and utilized
independently of B<Bio::Root::Object.pm>.

Since this module is at the root of potentially many different objects
in a particular application, efficiency is important. Bio::Root::Object.pm is 
intended to be a lightweight, lean and mean module. 


=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other Bioperl modules.
Send your comments and suggestions preferably to one of the Bioperl mailing lists.
Your participation is much appreciated.

    vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
    vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
    http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track the bugs and 
their resolution. Bug reports can be submitted via email or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           

=head1 AUTHOR

Steve A. Chervitz, sac@genome.stanford.edu

See the L<FEEDBACK> section for where to send bug reports and comments.

=head1 VERSION

Bio::Root::Object.pm, 0.041


=head1 TODO

=over 0

=item * Experiment with other Exception classes.

Consider incorporating a more widely-used Error/Exception module 
(see L<Other Exception Modules>).

=item * Think about integration with Data::Dumper.pm for persisting objects.

=back

=head1 SEE ALSO

  Bio::Root::Err.pm       - Error/Exception object
  Bio::Root::IOManager.pm - Input/Output manager object
  Bio::Root::Vector.pm    - Manages dynamic lists of objects
  Bio::Root::Xref.pm      - Cross-reference object
  Bio::Root::Global.pm    - Manages global variables/constants

  http://bio.perl.org/Projects/modules.html  - Online module documentation
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

Copyright (c) 1996-98 Steve A. Chervitz. All Rights Reserved.
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

See Also   : L<_initialize>(), L<name>(), L<parent>(), L<make>(), L<strict>(), L<verbose>(), L<throw>(), L<warn>(), L<record_err>()

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

See Also   : L<new>(), L<make>(), L<name>(), L<parent>(), L<strict>(), L<index>(), L<_rearrange>(), L<_set_clone>(), L<verbose>()

=cut

#----------------
sub _initialize {
#----------------
    local($^W) = 0;
    my($self, %param) = @_;
    
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

See Also   : L<destroy>() 

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
           : counting. If you create such such structures outside of
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
           : data member, searching all top-level array, hash, and scalar data members. 
           : It does not recurse through all levels of complex data members.
           : Subclasses could override this method to handle complex child data members 
           : for more optimal child searching. However, the version here is
           : probably sufficient for most situations.
           :
           : _drop_child() is called by Bio::Root::Object::destroy() for all objects
           : with parents.
 Status    : Experimental

See Also   : L<destroy>()

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


=head2 _rearrange

 Usage     : $object->_rearrange( array_ref, list_of_arguments)
 Purpose   : Rearranges named parameters to requested order.
 Example   : $self->_rearrange([qw(SEQUENCE ID DESC)],@param);
           : Where @param = (-sequence => $s, 
	   :                 -id       => $i, 
	   :	             -desc     => $d);
 Returns   : @params - an array of parameters in the requested order.
           : The above example would return ($s, $i, $d)
 Argument  : $order : a reference to an array which describes the desired
           :          order of the named parameters.
           : @param : an array of parameters, either as a list (in
           :          which case the function simply returns the list),
           :          or as an associative array with hyphenated tags
           :          (in which case the function sorts the values 
           :          according to @{$order} and returns that new array.)
	   :	      The tags can be upper, lower, or mixed case
           :          but they must start with a hyphen (at least the
           :          first one should be hyphenated.)
 Source    : This function was taken from CGI.pm, written by Dr. Lincoln
           : Stein, and adapted for use in Bio::Seq by Richard Resnick and
           : then adapted for use in Bio::Root::Object.pm by Steve A. Chervitz.
 Comments  : (SAC)
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
	   :    ('-name'=>'me', '-color' =>'blue');
           : This can be a bit cumbersome and I find not as readable
           : as using all uppercase, which is also fairly safe:
	   :    (-NAME=>'me', -COLOR =>'blue');
	   :
           : Personal note (SAC): I have found all uppercase tags to
           : be more managable: it involves less single-quoting,
           : the code is more readable, and there are no method naming conlicts.
           : Regardless of the style, it greatly helps to line
	   : the parameters up vertically for long/complex lists.

See Also   : L<_initialize>() 

=cut

#----------------'
sub _rearrange {
#----------------
    my($self,$order,@param) = @_;
    
    # If there are no parameters, we simply wish to return
    # an empty array which is the size of the @{$order} array.
    return ('') x $#{$order} unless @param;
    
    # If we've got parameters, we need to check to see whether
    # they are named or simply listed. If they are listed, we
    # can just return them. 

    return @param unless (defined($param[0]) && $param[0]=~/^-/); 
    
    # Tester
#    print "\n_rearrange() named parameters:\n";
#    my $i; for ($i=0;$i<@param;$i+=2) { printf "%20s => %s\n", $param[$i],$param[$i+1]; }; <STDIN>;

    # Now we've got to do some work on the named parameters.
    # The next few lines strip out the '-' characters which
    # preceed the keys, and capitalizes them.
    my $i;
    for ($i=0;$i<@param;$i+=2) {
	$param[$i]=~s/^\-//;
	$param[$i]=~tr/a-z/A-Z/;
    }
    
    # Now we'll convert the @params variable into an associative array.
    local($^W) = 0;  # prevent "odd number of elements" warning with -w.
    my(%param) = @param;
    
    my(@return_array);
    
    # What we intend to do is loop through the @{$order} variable,
    # and for each value, we use that as a key into our associative
    # array, pushing the value at that key onto our return array.
    my($key);
    
    foreach $key (@{$order}) {
	my($value) = $param{$key};
	delete $param{$key};
	push(@return_array,$value);
    }
    
#    print "\n_rearrange() after processing:\n";
#    my $i; for ($i=0;$i<@return_array;$i++) { printf "%20s => %s\n", ${$order}[$i], $return_array[$i]; } <STDIN>;

    return (@return_array);
}



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

See also   : L<has_name>()

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

See also   : L<name>()

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

See also   : L<destroy>()

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
           :       conditionals that simply check name() would fail incorrectly.

See also   : L<name>()

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

See also   : L<_initialize>(), L<clone>()

=cut

#----------'
sub make { 
#----------
    my $self = shift;
    if(@_) { $self->{'_make'} = shift; }
    $self->{'_make'} || 'default'; 
}



=head2 throw

 Purpose   : Generate, report, and set a fatal error on the object.
           : Uses Perl's die() function to report error data.
           : This does not invalidate the object but will crash the script
           : unless it is trapped with eval{}. 
           : (fatal = un-recoverable)'
 Usage     : $object->throw([arguments for _set_err()])
 Returns   : die()s with the contents of the error object in a string.
           : This string is human-readable and can be used to reconstruct
           : the Bio::Root::Err.pm object (e.g., new  Bio::Root::Err($@) ).
           :
           : The behavior of throw() is affected by the current verbosity
           : and strictness settings:
           : If verbosity is < 0, the stack trace is not printed.
           : If verbosity is > 0, all data including stack trace is shown and a
           :   system beep is issued.
           : If verbosity = 0, print all data but no beep (message, note, tech note,
           :  containment hierarchy, stack).
           : If strictness is less than -1, the throw() call is converted
           : into a warn() call.
 Argument  : Arguments for _set_err() 
 Comments  : Calling $self->throw() method creates a Bio::Root::Err.pm object.
           : There are two ways to generate errors:
           :   1) $object->throw(<ERROR DATA>);
           :   2) &Bio::Root::Err::throw($object, ERROR_DATA);
           : To use the second option, include the line use Bio::Root::Err qw(:std);
           : in your script or module. ERROR_DATA = arguments for _set_err().
           :
           : Some auxilliary issues:
           :   * It would be great if Perl could throw an object reference with die().
           :     This would permit more intelligent exception handlers. For now the
           :     Err object is reconstructed from the output of Err::string().
           :
           : All errors are reported to STDERR. 
           : Redirection to an alternate location for storing errors
           : can be achieved by redirecting STDERR manually [ open(STDERR, ">>filename") ],
           : or by using set_log_err().

See also   : L<_set_err>(), L<warn>(), L<strict>(), L<verbose>(), L<set_log_err>(), L<STRICTNESS & VERBOSITY>, B<Bio::Root::Global:strictness()>, B<Bio::Root::Global:verbosity()>

=cut

#----------
sub throw {
#----------
    my($self,@param) = @_;

#    printf "Throwing exception for object %s \"%s\"\n",ref($self),$self->name;
#    print "param: @param";<STDIN>;

    my $verbosity = $self->verbose;

    if($self->strict < -1) {
	# Convert throw() to warn()
	return $self->warn(@param);

    } else {
	if($verbosity < 0) {
	    # Low verbosity: no stack trace.
	    die $self->_set_err(@param)->string(-SHOW=>'msgnotechcontext', -CURRENT=>1);
	} elsif($verbosity > 0) {
	    # Extra verbosity: all data and beep.
	    die $self->_set_err(@param)->string(-BEEP=>1, -CURRENT=>1);
	} else {
	    # Default: all data (msg, note, tech, context, stack trace) but no beep.
	    die $self->_set_err(@param)->string(-CURRENT=>1);
	}
    }
    0;
}

=head2 warn

 Usage     : $object->warn([arguments for _set_err()])
 Purpose   : Generate, report, and set a recoverable error on the object.
 Returns   : Prints the contents of the error to STDERR and returns false (0).
           : The behavior of warn() is affected by the current verbosity
           : and strictness settings:
           : If verbose() is < 0, nothing is printed but a warning is still set.
           : If verbose() is > 0, the full error listing is shown
           :  (message, note, tech note, containment hierarchy, stack).
           : If verbosity  = 0, the message, note, and tech note are shown.
           : If the strict() indicator is greater than 1, warn() calls are 
           : converted into throw() calls.
 Argument  : Arguments for _set_err() 
 Comments  : The return value is experimental. Typically, warnings are not
           : programatically trappable: a method will issue a warning and 
           : then go about its business. By allowing
           : warn() calls to evaluate to zero, a method can halt execution 
           : by returning a warning to signal the warning without setting a 
           : fatal error on itself. Still, returning 0 does not guarantee
           : the exception will be noticed. This sort of polling-based
           : exception handling is generally frowned upon. Using throw()
           : and trapping any exceptions is highly recommended unless
           : the condition is truly inconsequential.
           :
           : All errors are reported to STDERR. 
           : Redirection to an alternate location for storing errors
           : can be achieved by redirecting STDERR manually [ open(STDERR, ">>filename") ],
           : or by using set_log_err().

See also   : L<_set_warning>(), L<throw>(), L<strict>(), L<verbose>(), L<set_log_err>(),  L<STRICTNESS & VERBOSITY>, B<Bio::Root::Global:strictness()>, B<Bio::Root::Global:verbosity()>

=cut

#---------
sub warn {
#---------
    my($self,@param) = @_;

    my $verbosity = $self->verbose;

    if($self->strict > 1) {
	# Convert warn() to throw()
	$self->throw(@param);

    } else {
	if($verbosity < 0 || $CGI) {
	    # Low verbosity or script is a cgi: don't print anything but set warning.
	    $self->_set_warning(@param);
	} elsif($verbosity > 0) {
	    # Extra verbosity: print all data and beep
	    print STDERR $self->_set_warning(@param)->string(-BEEP=>1, -CURRENT=>1);
	} else {
	    # Default: message and notes only. No beep.
	    print STDERR $self->_set_warning(@param)->string(-SHOW=>'msgnotech', -CURRENT=>1);
	}
    }
    0;
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

See also   : L<err_string>(), L<print_err>(), B<Bio::Root::Err.pm>::get_all

=cut

#----------
sub err {
#----------
    my( $self, $data, $delimit) = @_;

    return undef unless defined $self->{'_err'};
    
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

See also   : L<_set_err>(), B<Bio::Root::Global::record_err()>

=cut

#---------------
sub record_err {
#---------------
    my $self = shift;

    if (@_) { $self->{'_record_err'} = shift }
    return $self->{'_record_err'} || 0;
}


=head2 _set_err

 Purpose   : To create a Bio::Root::Err.pm object and optionally attach it
           : it to the current object.
 Usage     : This is an internal method and should not be called directly
           : $object->_set_err( msg)  
           : $object->_set_err( msg, note) 
           : $object->_set_err( msg, note, tech)
           : $object->_set_err( -MSG  =>"main message", 
	   :	                -TECH =>"technical note only")
           : $object->_set_err($object->err())  # Transfers pre-existing err
           : $object->_set_err()                # Re-sets an object's error state 
           :                                    # (Public method: clear_err())
 Example   : $self->throw("Data not found.");
           : To throw an error:
           : $myData eq 'foo' || return $self->throw("Data is not 'foo'.")
 Returns   : Object reference to the newly created Bio::Root::Err.pm object
           : via call to _set_err_state().
           :
 Argument  : @param may be empty, or contain a single error object, 
           : named parameters, or a list of unnamed parameters for 
           : building an Bio::Root::Err object.
           :   msg  = string, basic description of the exception.
           :   note = string, additional note to indicate cause or exception
           :          or provide information about how to fix/report it.
           :   tech = string, addition note with technical information
           :          of interest to developer.
           :
           : When using unnamed parameters, the number of items in @param 
           : is used as a "syntactic sugar" to indicate which fields in the 
           : err object to set (1 = msg, 2 = msg + note, 3 = msg + note + tech)
           : Calling _set_err() with no arguments clears the {'_err'} and 
           : {'_errState'} data members and destroys the Err object.
           : 
 Comments  : NEW VERSION: NOT ATTACHING EXCEPTIONS TO THE OBJECT.
           : Since exceptions are fatal, it is more expedient for the calling code
           : to handle them as they arise. Attaching exceptions to the objects
           : that generated them implies that the object assumes responsibility for 
           : any error it might throw, which is not usually appropriate and is 
           : difficult to manage.
           :
           : The new code now by default will not attach Err objects to the 
           : object that. Attaching Err objects can be enabled using the -RECORD_ERR
           : constructor option or the record_err() method. Bio::Root::Global::record_err()
           : turns on Err attaching for all objects in a script.
           :
           : Attaching exceptions to the objects that produced them is considered
           : non-standard and must be explicitly requested. This behavior might be 
           : useful in situations where one runs some code in an unsupervised 
           : setting and needs a means for reporting all warnings/errors later.
           :
           : One problem with attaching Err objects is that if an object is contained 
           : within another object, the containing object will not know about the 
           : warning unless it polls all of it contained objects, (bad design).
           : One could propagate the warning through the containment hierarchy
           : but the hierarchy may not be accessible to the objects themselves:
           : a given object may not know where it is contained (i.e, it may not 
           : have a parent).
           :    
           : * To transfer an error between objects, you can use 
           :   $self->warn($object->err) or $self->throw($object->err) or
           :   $self->_set_err($object->err) to not generate a warning.

See also   : L<_set_warning>(), L<err>(), L<warn>(), L<throw>(), L<record_err>(), L<_set_err_state>(), B<Bio::Root::Err.pm>

=cut

#--------------
sub _set_err {   
#--------------
    my($self, @param) = @_;
    
    require Bio::Root::Err; import Bio::Root::Err qw(:data);

    ## Re-set the object's 'Err' member and clear the 'ErrState'
    ## if no data was provided. 
    if(!@param) {
	if(ref($self->{'_err'})) {
	    $self->{'_err'}->remove_all;
	    %{$self->{'_err'}} = ();
	}
	undef $self->{'_err'};
	undef $self->{'_errState'};
	return;
    }
    
    # If client specifies any error field (ie. -MSG=>'Bad data') the 'custom' constructor
    # is used. If no error field is specified by name, the constructor is selected based on
    # the number of arguments in @param.

    local($^W) = 0;  
    my %err_fields = %Bio::Root::Err::ERR_FIELDS;  # prevents warnings
    my $constructor = 'custom';
    if( !grep exists $Bio::Root::Err::ERR_FIELDS{uc($_)}, @param) {
	$constructor = scalar @param;
    }
    
    $DEBUG and do{ print STDERR "\n$ID: setting error ${\ref($self)} with constructor = $constructor"; <STDIN>; };
    
    ## Adjust the constructor number if STACK_NUM was given in @param (for warnings).
    ## This is a bit of a hack, but will do for now.
    ## The constructor needs to be adjusted since it included the "-STACK_NUM=>#" data
    ## which increases the number of arguments by 2 and is not to be included in 
    ## the exception's data.

    my $stackNum = 3;  
    if($constructor =~ /\d/ and grep (/STACK_NUM/i, @param)) { 
	$constructor -= 2;  ## Since STACK_NUM data was appended to @_.
	$stackNum=$param[$#param];
    }

    my ($err);
    eval {
	local $_ = $constructor;
	## Switching on the number of items in @param (now in $_).
	SWITCH: {
	    ## Single argument: $param[0] is either an Err object or a message.
	    /1/ && do{ 	if((ref($param[0]) =~ /Err/)) {
#		print "$ID: cloning err for object ${\ref $self}, \"${\$self->name}\"\n";<STDIN>;
		$err = $param[0]->clone(); ## Cloning error object.
	    } else {
		$err = new Bio::Root::Err(-MSG     =>$param[0], 
					  -STACK   =>scalar($self->stack_trace($stackNum)),
					  -CONTEXT =>$self->containment,
				       # -PARENT =>$self
					  );
	    }
			last SWITCH; };
	    
	    ## Two arguments: Message and note data.
	    /2/ && do{ $err = new Bio::Root::Err(-MSG     =>$param[0], 
						 -NOTE    =>$param[1], 
						 -STACK   =>scalar($self->stack_trace($stackNum)),
						 -CONTEXT =>$self->containment,
						# -PARENT =>$self 
						 );
		       last SWITCH; };
	    
	    ## Three arguments: Message, note, and tech data.
	    /3/ && do{ $err = new Bio::Root::Err(-MSG     =>$param[0], 
						 -NOTE    =>$param[1],  
						 -TECH    =>$param[2], 
						 -STACK   =>scalar($self->stack_trace($stackNum)),
						 -CONTEXT =>$self->containment,
					      #   -PARENT =>$self 
						 );
		       last SWITCH; };
	    
	    ## Default. Pass arguments to Err object for custom construction.
	    ## Note: Stack data is not added. Should be provided in @_.
	    $err = new Bio::Root::Err( @param, 
				       -STACK   =>scalar($self->stack_trace($stackNum)),
				     #  -PARENT =>$self
				       );
	}
    };
    if($@) {
	printf STDERR "%s \"%s\": Failed to create Err object: \n $@", ref($self),$self->name;<STDIN>;
	print STDERR "\nReturning $self->{'_err'}:";<STDIN>;
	return $self->{'_err'}; }  ## Err construction will fail if the err is a duplicate.
                                  ## In any event, the Err object creation error is
                                  ## simply reported to STDERR and 
                                  ## the current 'Err' member is returned since 
                                  ## there is no call to _set_err_state().

    ## ONLY ATTACH THE Err OBJECT IF DIRECTED TO.
    if($RECORD_ERR or $self->{'_record_err'}) {
	if( ref($self->{'_err'})) {
	    ## Expand the linked list of Err objects if necessary.
#	    print "Expanding linked list of errors for ${\$self->name}.\nError = ${\$err->msg}\n";
	    $self->{'_err'}->last->add($err);
	} else {
#	    print "Adding new error for ${\$self->name}.\nError = ${\$err->msg}\n";
	    $self->{'_err'} = $err;
	}

	## Propagate the error up the containment hierarchy to ensure it gets noticed.
	## No longer needed since exception is fatal.
#	$self->_propagate_err($err);
    }
#    else { print "NOT SETTING ERR OBJECT: \n${\$err->string}\n\n"; }

    $self->_set_err_state($err);
}



=head2 _set_err_state

 Usage     : n/a; called automatically by _set_err()
 Purpose   : Sets the {'_errState'} data member to one of @Bio::Root::Err::ERR_TYPES.
           : This method is called after setting a new error with _set_err().
 Returns   : An Err.pm object (the current {'_err'} data member)
 Argument  : An Err.pm object (the one jsut created by _set_err()).
 Comments  : Modifications to state are permitted only if the object:
           :   1. has only one error, OR
           :   2. has multiple errors and none of those errors are fatal.
           : This prevents an object from setting its state to warning
           : if it already has a fatal error.
           :
           : The unfatal() method circumvents this method since the conditions
           : under which unfatal() is called are different. _set_err_state() is
           : only called when setting new errors.

See also   : L<_set_err>(), L<_set_warning>() 

=cut

#--------------------
sub _set_err_state {  
#--------------------
    my( $self, $err ) = @_;
    my @state = ();
    
    require Bio::Root::Err; import Bio::Root::Err qw(:data);
    
    my $data = $err->type || 'EXCEPTION';

    if($self->{'_errState'} and $self->{'_errState'} !~ /EXCEPTION|FATAL/) {
	
	my @err_types = @Bio::Root::Err::ERR_TYPES; # prevents warnings
	if( @state = grep /$data/i, @Bio::Root::Err::ERR_TYPES ) {
	    $self->{'_errState'} = $state[0];
	} else {
	    $self->{'_errState'} = 'UNKNOWN STATE';
	}
    }
    $DEBUG and do{ print STDERR "$ID: Setting state to $self->{'_errState'} (arg=$data)\n"; <STDIN>; };

#    $self->{'_err'}->last;
    return $err;
}



=head2 err_state

 Usage     : $object->err_state();
 Purpose   : To assess the status of the object's Err object (if any).
 Returns   : A string: 'EXCEPTION' | 'WARNING' | 'FATAL' | 'OKAY'
           : (OKAY is returned if there are no Errors)
 Status    : Experimental

See also   : L<_set_err_state>()

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

See  also  : L<_set_err>()

=cut

#-------------
sub clear_err {
#-------------
    my $self = shift;
    $self->_set_err();
}


=head2 _set_warning

 Purpose   : To record data regarding recoverable error conditions.
 Usage     : n/a; called automatically by Bio::Root::Object::warn() 
 Arguments : Arguments passed as-is to _set_err().
 Comments  : An object with a warning should be considered 
           : completely operational, so use this type of error sparingly. 
           : These errors are intended for problem conditions which:
           :  1. Don't destroy the basic functionality of the object.
           :  2. Might be of incidental interest to the user.
           :  3. Are of interest to the programmer but not the end user.

See also   : L<warn>(), L<_set_err>(), L<err>()

=cut

#-----------------'
sub _set_warning {  
#-----------------
    my( $self, @data ) = @_;  
    $DEBUG and print STDERR "\n$ID: setting warning.\n";
    
    my $err = $self->_set_err(@data, -STACK_NUM=>4);
    $err->last->set('type','WARNING');
    $self->_set_err_state($err);
}



=head2 stack_trace

 Usage     : $stack_aref = $myObj->stack_trace([start_index, [end_index]]);
           : @stack_list = $myObj->stack_trace([start_index, [end_index]]);
           : @stack_list = stack_trace($object);  # As an exported method.
 Purpose   : Returns the contents of the current call stack
           : in a slightly modified, more intuitive form. 
           : Permits extraction of a portion of the stack. 
           : Call stack is obtained from the perl caller() function.
           : MODIFIED FORMAT: Line numbers are shifted down one
           : level in the stack entries so that they correspond 
           : to the location of the indicated method.
 Example   : @stackData = $self->stack_trace(2); 
 Argument  : start_index : number of the beginning entry in the
           :               desired stack trace. The call to stack_trace()
           :               is at index 0.
           : end_index   : number of the last entry in the
           :               desired stack trace.
 Returns   : A list or list reference (depending on wantarray)
           : consisting of the desired portion of the call stack.

See also   : B<Bio::Root::Err::format_stack_entry()>

=cut

#-----------------
sub stack_trace {
#-----------------
    my($self,$beg,$end) = @_;
    my(@call,@data);

    ## Set the complete stack trace.
    my $i = 0;
    while( @call = caller($i++)) { 
	my @callData = @call; 
#	print "CALL DATA $i: \n";
#	my $j = 0; for($j=0; $j<@callData; $j++) {print "$j: $callData[$j]\n"; }; <STDIN>;
	next if $callData[3] eq '(eval)';  ## Screening out the (eval) calls.
	push @data, \@callData;
    }

    ## Shift the line numbers down so that they correspond to 
    ## the location of the shown method. This is more intuitive.
    ## Processing stack in reverse.
    my( @base_call, $temp);
    for($i=$#data; $i > 0; $i--) {
	$temp = $data[$i]->[2];
	$data[$i]->[2] = $data[$i-1]->[2];
	if($i == $#data) { @base_call = @{$data[$i]}; 
			   $base_call[2] = $temp;
			   $base_call[3] = "$data[$i]->[0]::$data[$i]->[1]";
		       }
    }
    @data = (@data, \@base_call);
#    print "FULL STACK:\n";foreach(@data){print "@$_\n";};<STDIN>;

    ## Get everything but the call to stack_trace
    $beg ||= 1;
    $end ||= $#data;
    @data = @data[$beg..$end];

    wantarray ? @data : \@data;
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

See also   : L<_set_err>()

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

See also   : L<warn>(), L<throw>(), L<STRICTNESS & VERBOSITY>, B<Bio::Root::Global::strictness()>

=cut

#------------
sub strict {
#------------
    my $self = shift;

    # Using global strictness
    return &strictness(@_);

    # Object-specific strictness (not used unless above code is commented out)
    if (@_) { $self->{'_strict'} = shift; }
    defined($self->{'_strict'}) 
	? return $self->{'_strict'}
	: (ref $self->{'_parent'} ? $self->{'_parent'}->strict : 0);
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

See Also   : L<_set_clone>(), L<_initialize>()

=cut

#-------------'
sub clone {
#-------------
    my($self) = shift; 

#    warn sprintf "\nCloning %s \"%s\"\n\n", ref($self),$self->name;

    my $clone = $self->new(-MAKE    =>'clone', 
			   -OBJ     =>$self);  
    if($self->err()) { $clone->_set_err($self->err); } 
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

See Also   : L<_initialize>()

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

See Also   : L<strict>(), L<STRICTNESS & VERBOSITY>, B<Bio::Root::Global::verbosity()>

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



=head1 I/O-RELATED METHODS (Delegated to B<Bio::Root::IOManager.pm>)

=head2 _io

 Usage     : $object->_io()
 Purpose   : Get the Bio::Root::IOManager.pm object for the current object.

See also   : L<display>(), L<read>(), L<file>()

=cut

#-------
sub _io  {  my $self = shift;   return $self->{'_io'}; }
#-------



=head2 _set_io

 Usage     : n/a; internal use only.
 Purpose   : Sets a new Bio::Root::IOManager.pm object for the current object.

See also   : L<display>(), L<read>(), L<file>()

=cut

#------------
sub _set_io {
#------------
    my $self = shift;
    
    require Bio::Root::IOManager;

    $self->{'_io'} = new Bio::Root::IOManager(-PARENT=>$self, @_);
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

See also   : B<Bio::Root::IOManager.pm>::set_display

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

See also   : B<Bio::Root::IOManager::display()>

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

See also   : L<display>()

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
    printf ( $OUT "%-15s: %d\n",'ERRORS',    $self->err('count'));
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
 
See also   : L<display>(), B<Bio::Root::IOManager::read()>

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

See also   : B<Bio::Root::IOManager::fh()>

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

See also   : B<Bio::Root::IOManager::show()>, B<Bio::Root::IOManager::set_display()>

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

See also   : B<Bio::Root::IOManager::file()>

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

See also   : L<file>(), B<Bio::Root::IOManager::compress_file()>

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

See also   : L<file>(), B<Bio::Root::IOManager::uncompress_file()>

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

See also   : L<file>(), B<Bio::Root::IOManager.pm>::delete_file

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

See also   : L<file>(), B<Bio::Root::IOManager::file_date()>

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

See Also   : B<Bio::Root::Xref.pm>

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

See also   : L<find_object>()

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

See  also  : L<index>()

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

See also   : L<err>(), L<warn>(), L<throw>()

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

See also   : L<err_string>(), L<strict>()

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

See also   : L<print_err>(), B<Bio::Root::Err.pm>::string

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

See also   : L<verbose>()

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

See Also   : B<Bio::Root::Err.pm>set

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

See also   : B<Bio::Root::IOManager.pm>::set_read

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

See also   : B<Bio::Root::IOManager.pm>::set_log_err

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
    descision, see comments under _set_err().
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


