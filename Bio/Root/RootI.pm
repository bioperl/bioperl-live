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
# 
# This was refactored to have chained calls to new instead
# of chained calls to _initialize
#
# added debug and deprecated methods --Jason Stajich 2001-10-12
# 
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

use vars qw(@ISA $DEBUG $ID $Revision $VERSION $VERBOSITY);
use strict;
#use Bio::Root::Err; # we don't use that any longer, right?

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
    local($^W) = 0;
    my ($caller, @args) = @_;
    
    my $class = ref($caller) || $caller; #copied from Conway, OOPerl
    my $self = bless({}, $class);

    my %param = @args;
    my($verbose) =  ( $param{'-VERBOSE'} || $param{'-verbose'} );

    ## See "Comments" above regarding use of _rearrange().
    $self->verbose($verbose);

    return $self;
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

   my $out = "-------------------- EXCEPTION --------------------\n".
       "MSG: ".$string."\n".$std."-------------------------------------------\n";
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
 Example :
 Returns : 
 Args    :

=cut

sub warn{
    my ($self,$string) = @_;

    my $verbose = $self->verbose;

    if( $verbose == 2 ) {
	$self->throw($string);
    } elsif( $verbose == -1 ) {
	return;
    } elsif( $verbose == 1 ) {
	my $out = "-------------------- WARNING ---------------------\n".
		"MSG: ".$string."\n";
	$out .= $self->stack_trace_dump;
	
	print STDERR $out;
	return;
    }    

    my $out = "-------------------- WARNING ---------------------\n".
       "MSG: ".$string."\n".
	   "---------------------------------------------------\n";
    print STDERR $out;
}

=head2 debug

 Title   : debug
 Usage   : $obj->debug("This is debugging output");
 Function: Prints a debugging message when verbose is > 0
 Returns : none
 Args    : message string to print to STDERR

=cut

sub debug{
   my ($self,$msg) = @_;
   
   if( $self->verbose > 0 ) { 
       print STDERR $msg;
   }   
}

=head2 deprecated

 Title   : deprecated
 Usage   : $obj->deprecated("Method X is deprecated");
 Function: Prints a message about deprecation 
           unless verbose is < 0 (which means be quiet)
 Returns : none
 Args    : Message string to print to STDERR

=cut

sub deprecated{
   my ($self,$msg) = @_;
   if( $self->verbose >= 0 ) { 
       print STDERR $msg, "\n", $self->stack_trace_dump;
   }
}

		     
=head2 verbose

 Title   : verbose
 Usage   : $self->verbose(1)
 Function: Sets verbose level for how ->warn behaves
           -1 = no warning
            0 = standard, small warning
            1 = warning with stack trace
            2 = warning becomes throw
 Returns : nothing
 Args    : -1,0,1 or 2


=cut

sub verbose{
   my ($self,$value) = @_;

   if(ref($self) && (defined $value || ! defined $self->{'_rootI_verbose'}) ) {
       $value = 0 unless defined $value;
       $self->{'_rootI_verbose'} = $value;
   }
   return (ref($self) ? $self->{'_rootI_verbose'} : $VERBOSITY);
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
   my @out;
   my $prev;
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

=head2 _register_for_cleanup

 Title   : _register_for_cleanup
 Usage   : -- internal --
 Function: Register a method to be called at DESTROY time. This is useful
           and sometimes essential in the case of multiple inheritance for
           classes coming second in the sequence of inheritance.
 Returns : 
 Args    : a reference to a method


=cut

sub _register_for_cleanup {
    my ($self,$method) = @_;
    if($method) {
	if(! exists($self->{'_cleanup_methods'})) {
	    $self->{'_cleanup_methods'} = [];
	}
	push(@{$self->{'_cleanup_methods'}},$method);
    }
}

sub DESTROY {
    my ($self) = shift;
    if(exists($self->{'_cleanup_methods'})) {
	foreach my $method (@{$self->{'_cleanup_methods'}}) {
	    &$method($self);
	}
    }
}

1;
