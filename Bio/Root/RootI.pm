
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
use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::Root::Err;

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

    my $verbosity;

    if( $self->can('verbose') ) {
	$verbosity = $self->verbose;
    } else {
	$verbosity = 0;
    }


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
	$verbosity = $self->verbose;
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

1;





