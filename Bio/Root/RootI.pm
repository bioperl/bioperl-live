# $Id$
#
# BioPerl module for Bio::Root::RootI
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Root::RootI - Abstract interface to root object code

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

=head1 DESCRIPTION

This is just a set of methods which do not assumme B<anything> about the object
they are on. The methods provide the ability to throw exceptions with nice
stack traces.

This is what should be inherieted by all bioperl compliant interfaces, even
if they are exotic XS/CORBA/Other perl systems.

=head1 CONTACT

Functions originally from Steve Chervitz. Refactored by Ewan Birney.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Root::RootI;
use vars qw(@ISA $DEBUG $ID $Revision $VERSION $VERBOSITY 
	    $TEMPMODLOADED $TEMPDIR $TEMPCOUNTER $OPENFLAGS);
use strict;
use Bio::Root::Err;
# determine tempfile
$TEMPCOUNTER = 0;
eval { require 'File/Temp.pm'; $TEMPMODLOADED = 1; };
if( $@ ) { 
    use Fcntl;
    use Symbol;
    $OPENFLAGS = O_CREAT | O_EXCL | O_RDWR;
    for my $oflag (qw/FOLLOW BINARY LARGEFILE EXLOCK NOINHERIT TEMPORARY/) {
	my ($bit, $func) = (0, "Fcntl::O_" . $oflag);
	no strict 'refs';
	$OPENFLAGS |= $bit if eval { $bit = &$func(); 1 };
    }
    eval { require 'File/Spec.pm'; $TEMPDIR = File::Spec->tmpdir(); };
    if( $@ ) {
	if (defined $ENV{'TEMPDIR'} && -d $ENV{'TEMPDIR'} ) {
	    $TEMPDIR = $ENV{'TEMPDIR'};
	} elsif( defined $ENV{'TMPDIR'} && -d $ENV{'TMPDIR'} ) {
	    $TEMPDIR = $ENV{'TMPDIR'};
	} elsif ( -d '/tmp' && -w '/tmp' ) {
	    $TEMPDIR = '/tmp';
	} else {
	    $TEMPDIR = '.';
	}
    }
}

BEGIN { 

    $ID        = 'Bio::Root::RootI';
    $VERSION   = 0.7;
    $Revision  = '$Id$ ';
    $DEBUG     = 0;
    $VERBOSITY = 0;
}


=head2 new

 Purpose   : generic intantiation function can be overridden if 
             special needs of a module cannot be done in _initialize
 
=cut

sub new {
    my ($caller, @args) = @_;
    
    my $caller_is_obj = ref($caller);      #Dave Block
    my $class = $caller_is_obj || $caller; #copied from Conway, OOPerl

    my $self = bless {}, $class;
    eval { 
	$self->_initialize(@args);   
    };
    if( $@ ) {
	&throw(new Bio::Root::RootI, $@);
    }
    return $self;
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

    my $verbosity = 0;


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

    my $verbosity;

    if( $self->can('verbose') ) {
	$verbosity = $self->verbose();
    } else {
	$verbosity = 0;
    }

    if($verbosity < 0 ) {
	# Low verbosity or script is a cgi: don't print anything but set warning.	
	$self->_set_warning(@param);
    } elsif($verbosity > 0) {
	# Extra verbosity: print all data and beep
	print STDERR $self->_set_warning(@param)->string(-BEEP=>1, -CURRENT=>1);
    } else {
	# Default: message and notes only. No beep.
	print STDERR $self->_set_warning(@param)->string(-SHOW=>'msgnotech', -CURRENT=>1);
    }
    0;
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
    
    my $err = $self->_set_err(@data, -STACK_NUM=>4);
    $err->last->set('type','WARNING');
    $self->_set_err_state($err);
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
    my $name;

    # comment from EB. I think this can be slimmed/de-obfustecated
    
    local($^W) = 0;  
    my %err_fields = %Bio::Root::Err::ERR_FIELDS;  # prevents warnings
    my $constructor = 'custom';
    if( !grep exists $Bio::Root::Err::ERR_FIELDS{uc($_)}, @param) {
	$constructor = scalar @param;
    }
    
    if( $self->can('name') ) {
	$name = $self->name;
    } else {
	$name = "Anonymous";
    }

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
		$err = $param[0]->clone(); ## Cloning error object.
	    } else {
		$err = new Bio::Root::Err(-MSG     =>$param[0], 
					  -STACK   =>scalar($self->stack_trace($stackNum)),
					  );
	    }
			last SWITCH; };
	    
	    ## Two arguments: Message and note data.
	    /2/ && do{ $err = new Bio::Root::Err(-MSG     =>$param[0], 
						 -NOTE    =>$param[1], 
						 -STACK   =>scalar($self->stack_trace($stackNum)),
						 );
		       last SWITCH; };
	    
	    ## Three arguments: Message, note, and tech data.
	    /3/ && do{ $err = new Bio::Root::Err(-MSG     =>$param[0], 
						 -NOTE    =>$param[1],  
						 -TECH    =>$param[2], 
						 -STACK   =>scalar($self->stack_trace($stackNum)),
						 );
		       last SWITCH; };
	    
	    ## Default. Pass arguments to Err object for custom construction.
	    ## Note: Stack data is not added. Should be provided in @_.
	    $err = new Bio::Root::Err( @param, 
				       -STACK   =>scalar($self->stack_trace($stackNum)),
				       );
	}
    };
    if($@) {
	
	printf STDERR "%s \"%s\": Failed to create Err object: \n $@", ref($self),$name;<STDIN>;
	print STDERR "\nReturning $self->{'_err'}:";<STDIN>;
	return $self->{'_err'}; 
    }  

    ## Err construction will fail if the err is a duplicate.
    ## In any event, the Err object creation error is
    ## simply reported to STDERR and 
    ## the current 'Err' member is returned since 
    ## there is no call to _set_err_state().
    

#   EB - removed this.
#    $self->_set_err_state($err);

    return $err;
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
    #$DEBUG and do{ print STDERR "$ID: Setting state to $self->{'_errState'} (arg=$data)\n"; <STDIN>; };

#    $self->{'_err'}->last;
    return $err;
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

#------------
sub verbose { 
#------------
    my ($self,$value) = @_; 

    # Using global verbosity
    if( defined $value ) {
	$VERBOSITY = $value;
    }
    return $VERBOSITY;

    # Object-specific verbosity (not used unless above code is commented out)
    if(@_) { $self->{'_verbose'} = shift; }
    defined($self->{'_verbose'}) 
	? return $self->{'_verbose'}
	: (ref $self->{'_parent'} ? $self->{'_parent'}->verbose : 0);
}

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
	$self->verbose($verbose) || undef;
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
    
    # JGRG -- This is wrong, because we don't want
    # to assign empty string to anything, and this
    # code is actually returning an array 1 less
    # than the length of @param:

    ## If there are no parameters, we simply wish to return
    ## an empty array which is the size of the @{$order} array.
    #return ('') x $#{$order} unless @param;
    
    # ...all we need to do is return an empty array:
    return unless @param;
    
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

# from (old) Bio::Root::Object

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
    return defined $self->{'_name'} ? $self->{'_name'} : 
	'anonymous '.ref($self);
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

=head2  tempfile

 Title   : tempfile
 Usage   : my ($handle,$tempfile) = $self->tempfile($dir); 
 Function: returns a temporary file for writing
 Example : my ($handle,$tempfile) = $self->tempfile($dir);
 Returns : temporary handle and temporary file name
 Args    : $dir - directory in which to create new file 

=cut

sub tempfile {
    my ($self, @args) = @_;
    if ( $TEMPMODLOADED ) {
	my ($tfh, $file) = &File::Temp::tempfile(@args);
	push @{$self->{'_rooti_tempfiles'}}, $file;
	return ($tfh,$file);
    }
    my %hash = @args;
    my $dir;
    $dir = ( $hash{DIR} ) ? $hash{DIR} : $dir = $TEMPDIR;

    my $tfilename = sprintf("%s/%s-%s-%s", $dir, 
			    $ENV{USER} || 'unknown', $$, 
			    $TEMPCOUNTER++);
    push @{$self->{'_rooti_tempfiles'}}, $tfilename;

    # taken from File::Temp;
    my $fh;
    if ($] < 5.006) {
        $fh = &Symbol::gensym;
    }
    
    # Try to make sure this will be marked close-on-exec
    # XXX: Win32 doesn't respect this, nor the proper fcntl,
    #      but may have O_NOINHERIT. This may or may not be in Fcntl.
    local $^F = 2; 
    
    # Store callers umask
    my $umask = umask();
    
    # Set a known umaskr
    umask(066);
    
    # Attempt to open the file
    if ( sysopen($fh, $tfilename, $OPENFLAGS, 0600) ) {
	
        # Reset umask
        umask($umask); 
        # Opened successfully - return file handle and name
        return ($fh, $tfilename);
    } else { 
	$self->throw("Could not open tempfile $tfilename $!\n");
    }
}

=head2  tempdir

 Title   : tempdir
 Usage   : my ($tempdir) = $self->tempdir(CLEANUP=>1); 
 Function: returns a temporary directory
 Example : my ($tempdir) = $self->tempdir(CLEANUP=>1); 
 Returns : a temporary directory
 Args    : hash - ( key CLEANUP ) indicates whether or not to cleanup 
           dir on object destruction
=cut

sub tempdir {
    my ( $self, %hash ) = @_;

    if( $TEMPMODLOADED ) {
	my $dir = &File::Temp::tempdir(%hash);
	push @{$self->{'_rooti_tempdirs'}},$dir;
	return $dir;
    }
    # we are planning to cleanup temp files no matter what
    if( $hash{CLEANUP} == 1 ) {	$self->{'_cleanuptempdir'} = 1;
    }
    
    my $tdir = sprintf("%s/%s-%s-%s", $TEMPDIR, 
		    "dir_". $ENV{USER} || 'unknown', $$, 
		    $TEMPCOUNTER++);
    mkdir($tdir, 0755);
    push @{$self->{'_rooti_tempdirs'}}, $tdir; 
    return $tdir;
}


sub DESTROY {
    my ($self) = @_;
    # we are planning to cleanup temp files no matter what
    unlink @{$self->{'_rooti_tempfiles'}} 
    if( defined $self->{'_rooti_tempfiles'} 
	&& ref($self->{'_rooti_tempfiles'}) =~ /array/i );
    foreach ( @{$self->{'_rooti_tempdirs'}} ) {
#	rmdir($_);
    }
}
1;
