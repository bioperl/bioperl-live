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
	    $TEMPMODLOADED $FILESPECLOADED $TEMPDIR $TEMPCOUNTER $OPENFLAGS);
use strict;
use Bio::Root::Err;
# determine tempfile
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
    eval { require 'File/Spec.pm'; 
	   $FILESPECLOADED = 1;
	   $TEMPDIR = File::Spec->tmpdir(); 
       };
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
    $TEMPCOUNTER = 0;
}


=head2 new

 Purpose   : generic intantiation function can be overridden if 
             special needs of a module cannot be done in _initialize
 
=cut

sub new {
    local($^W) = 0;
    my ($caller, @args) = @_;
    
    my $caller_is_obj = ref($caller);      #Dave Block
    my $class = $caller_is_obj || $caller; #copied from Conway, OOPerl

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

   if( defined $value || ! defined $self->{'_rootI_verbose'} ) {
       $value = 0 unless defined $value;
       $self->{'_rootI_verbose'} = $value;
   } 
   return $self->{'_rootI_verbose'};
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
    my $tfilename;

    if( $FILESPECLOADED ) {
	$tfilename = File::Spec->catfile($dir, sprintf( "%s-%s-%s",  
					$ENV{USER} || 'unknown', $$, 
							$TEMPCOUNTER++));
    } else {
	$tfilename = join("/", ($dir, sprintf( "%s-%s-%s",  
					       $ENV{USER} || 'unknown', $$, 
					       $TEMPCOUNTER++)));
    }
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
 Args    : args - ( key CLEANUP ) indicates whether or not to cleanup 
           dir on object destruction, other keys as specified by File::Temp

=cut

sub tempdir {
    my ( $self, %args ) = @_;

    if( $TEMPMODLOADED ) {
	my $dir = &File::Temp::tempdir(%args);
	return $dir;
    }
    # we are planning to cleanup temp files no matter what
    $self->{'_cleanuptempdir'} = $args{CLEANUP} == 1;

    my $tdir;
    if( $FILESPECLOADED ) {
	$tdir = File::Spec->catdir($TEMPDIR,
				   sprintf("dir_%s-%s-%s", 
					   $ENV{USER} || 'unknown', $$, 
					   $TEMPCOUNTER++));
    } else { 
	$tdir = join('/', ($TEMPDIR,
			   sprintf("dir_%s-%s-%s", 
				   $ENV{USER} || 'unknown', $$, 
				   $TEMPCOUNTER++)
			   ));
    }
    mkdir($tdir, 0755);
    push @{$self->{'_rooti_tempdirs'}}, $tdir; 
    return $tdir;
}


sub DESTROY {
    my ($self) = @_;
    # we are planning to cleanup temp files no matter what     
    if( defined $self->{'_rooti_tempfiles'} 
	&& ref($self->{'_rooti_tempfiles'}) =~ /array/i) { 
	unlink @{$self->{'_rooti_tempfiles'}}; 
    }
    # cleanup if we are not using File::Temp
    if( $self->{'_cleanuptempdir'} ) {
	foreach ( @{$self->{'_rooti_tempdirs'}} ) {
	    rmdir($_); 
	}
    }
}
1;



