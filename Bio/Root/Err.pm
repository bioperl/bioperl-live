#-----------------------------------------------------------------------------
# PACKAGE : Bio::Root::Err.pm
# AUTHOR  : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED : 22 July 1996
# REVISION: $Id$
# STATUS  : Alpha
# 
# For documentation, run this module through pod2html 
# (preferably from Perl v5.004 or better).
#
# Copyright (c) 1996-8 Steve A. Chervitz. All Rights Reserved.
#           This module is free software; you can redistribute it and/or 
#           modify it under the same terms as Perl itself.
#           Retain this notice and note any modifications made.
#-----------------------------------------------------------------------------

package Bio::Root::Err;
use strict;

use Bio::Root::Global  qw(:devel $CGI);
use Bio::Root::Vector  ();
use Bio::Root::Object  qw();#:std);
use Exporter           ();

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
@ISA         = qw( Bio::Root::Object Bio::Root::Vector Exporter );
@EXPORT      = qw();
@EXPORT_OK   = qw( %ERR_FIELDS @ERR_TYPES &format_stack_entry &throw &warning);
%EXPORT_TAGS = ( 
		 data => [qw(%ERR_FIELDS @ERR_TYPES)],
		 std  => [qw(&throw &warning)]
		);

use vars qw($ID $VERSION);
$ID = 'Bio::Root::Err';
$VERSION = 0.041;

%Bio::Root::Err::ERR_FIELDS = (TYPE=>1, MSG=>1, NOTE=>1, CONTEXT=>1,
			       TECH=>1, STACK=>1 );

@Bio::Root::Err::ERR_TYPES = qw(WARNING EXCEPTION FATAL);


## MAIN POD DOCUMENTATION:

=head1 NAME

Bio::Root::Err.pm -  Exception class for Perl 5 objects

=head1 SYNOPSIS

=head2 Object Creation

B<Bio::Root::Object.pm> is a wrapper for Bio::Root::Err.pm objects so clients 
do not have to create these objects directly. Please see 
B<Bio::Root::Object::throw()> as well as L<_initialize>()
for a more complete treatment 
of how to create Bio::Root::Err.pm objects.

  use Bio::Root::Err;
 
  $err = Bio::Root::Err->new(-MSG     =>"Bad data: $data", 
			     -STACK   =>[\caller(0), \caller(1), ...],
			     );
  

To use the L<throw>() method directly:

  use Bio::Root::Err (:std);

  throw( $object_ref, 'Error message', 'additional note', 'technical note');

The C<$object_ref> argument should be a reference to a Bio::Root::Object.pm.

See also L<USAGE>.

=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.

=head1 DESCRIPTION

A Bio::Root::Err.pm object encapsulates data and methods that facilitate 
working with errors and exceptional conditions that arise in Perl objects.
There are no biological semantics in this module, as one may suspect from its
location in the Bio:: hierarchy. The location of this module serves to
separate it from the namespaces of other Perl Error modules. It also makes it convenient
for use by Bio:: objects.

The motivation for having an error object is to allow 
Perl 5 objects to deal with errors or exceptional conditions that
can arise during their construction or manipulation. For example:

 (1) A complex object can break in many ways.
 (2) Tracking errors within a set of nested objects can be difficult. 
 (3) The way an error is reported should be context-sensitive:
     a web-user needs different information than does the
     software engineer. 

Bio::Root::Err.pm, along with B<Bio::Root::Object.pm>, attempt to make such problems 
tractable. Please see the B<Bio::Root::Object.pm> documentation for more about 
my error handling philosophy.

A B<Bio::Root::Err.pm> object is an example of a Vector-Object: This module inherits
both from B<Bio::Root::Object.pm> and B<Bio::Root::Vector.pm>. This permits a single Err
object to exist within a linked list of Err objects OR alone.
See the B<Bio::Root::Vector.pm> documentation for more about Vector-Objects.

B<The API for this module is not complete since the module is under development.>

=head2 Other Exception Strategies

Exception handling with Perl 5 objects is currently not as evolved as one
would like. The error handling used by B<Bio::Root::Object.pm> and Bio::Root::Err.pm
relies on Perl's built-in error/exception handling with eval/die, 
which is not very object-aware. What I've attempted to do with these
modules is to make eval/die more object-savvy, as well as make Perl 5 
objects more eval/die-savvy (but the current strategy is basically a hack).

It would be great if Perl could throw an object reference with die().
This would permit more intelligent and easy to write exception handlers. 
For now the Err.pm object is reconstructed from the output of L<string>().

There are some other third-party Exception classes such as 
Torsten Ekedahl's B<Experimental::Exception.pm> or
Ken Steven's Throwable.pm or
Graham Barr's Error.pm (see L<Other Exception Modules>). These modules 
attempt to introduce a traditional "try-catch-throw" 
exception handling mechanism into Perl. Future version of my modules  
(and perhaps erl itself) may utilize one of these.

=head1 USAGE

A demo script that illustrates working with Bio::Root::Err objects is available at:

    http://bio.perl.org/Core/Examples/Root_object/error.pl


=head1 DEPENDENCIES

Bio::Root::Err.pm inherits from B<Bio::Root::Object.pm> and B<Bio::Root::Vector.pm>.


=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other Bioperl modules.
Send your comments and suggestions preferably to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org             - General discussion
  bioperl-guts-l@bioperl.org        - Technically-oriented discussion
  http://bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track the bugs and 
their resolution. Bug reports can be submitted via email or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           

=head1 AUTHOR

Steve A. Chervitz, sac@genome.stanford.edu

See the L<FEEDBACK> section for where to send bug reports and comments.

=head1 VERSION

Bio::Root::Err.pm, 0.041


=head1 SEE ALSO

  Bio::Root::Object.pm    - Core object
  Bio::Root::Vector.pm    - Vector object
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

Copyright (c) 1996-8 Steve A. Chervitz. All Rights Reserved.
This module is free software; you can redistribute it and/or 
modify it under the same terms as Perl itself.

=head1 TODO

=over 2

=item * Improve documentation.

=item * Experiment with other Exception modules.

=back

=cut

## END MAIN POD DOCUMENTATION'

=head1 APPENDIX

Methods beginning with a leading underscore are considered private
and are intended for internal use by this module. They are
B<not> considered part of the public interface and are described here
for documentation purposes only.

=cut


########################################################
#                CONSTRUCTOR                           #
########################################################


=head2 _initialize

 Usage     : n/a; automatically called by Bio::Root::Object::new()
 Purpose   : Initializes key Bio::Root::Err.pm data.
 Returns   : String (the -MAKE constructor option.)
 Argument  : Named parameters passed from new()
           : (PARAMETER TAGS CAN BE UPPER OR LOWER CASE).
           :   -MSG     => basic description of the exception.
           :   -NOTE    => additional note to indicate cause of exception
           :               or provide information about how to fix/report it
           :   -TECH    => addition note with technical information
           :               of interest to developer.
           :   -STACK   => array reference containing caller() data
           :   -TYPE    => string, one of @Bio::Root::Err::ERR_TYPES
           :               (default = exception).
           :   -CONTEXT => array reference
           :   -OBJ     => Err object to be cloned.

See Also   : B<Bio::Root::Object::_set_err()>

=cut

#----------------
sub _initialize {
#----------------
    my( $self, @param ) = @_;
    
    my $make = $self->Bio::Root::Object::_initialize( @param );
    
    my( $msg, $note, $tech, $stack, $type, $context, $obj) = 
	$self->_rearrange([qw(MSG NOTE TECH STACK TYPE CONTEXT OBJ)], @param);

    ## NOTE: Don't eval {} the construction process for Err objects.

    if($make =~ /clone/i) { 
	$self->_set_clone($obj);
    } else {
	if(!$self->_build_from_string($msg, $note, $tech)) {
#	    print "Failed to rebuild: msg = $msg";<STDIN>;
	    $self->set('msg', $msg );
	    $self->_set_type( $type ); 
	    $self->_set_context($context);
	    $self->_set_list_data('note', $note );
	    $self->_set_list_data('tech', $tech );
	    $self->_set_list_data('stack', $stack );
	}
	$self->set_display();  
    }

    $DEBUG and do{ print STDERR "---> Initialized Err (${\ref($self)}).\n\n";
		   # $self->print(); 
	       };
    $make;
}

##
## Destructor: Not needed currently. Perhaps if and when Vector is used by delegation.
##

#####################################################################################
##                                  ACCESSORS                                      ##
#####################################################################################



=head2 _set_clone

 Usage     : n/a; internal method used by _initialize()
 Purpose   : Copy all Bio::Root::Err.pm data members into a new object reference.
 Argument  : object ref for object to be cloned.
 Comments  : Does not cloning the vector since this method is
           : typically used to extract a single Err object from its vector.

=cut

#---------------
sub _set_clone {
#---------------
    my($self, $obj) = @_;

    ref($obj) || throw($self, "Can't clone $ID object: Not an object ref ($obj)");

    $self->{'_type'}  = $obj->{'_type'};
    $self->{'_msg'}   = $obj->{'_msg'};
    $self->{'_note'}  = $obj->{'_note'};
    $self->{'_tech'}  = $obj->{'_tech'};
    $self->{'_stack'} = $obj->{'_stack'};
    $self->{'_context'} = $obj->{'_context'};
#    $self->clone_vector($obj);
}


=head2 _build_from_string

 Usage     : n/a; called by _initialize()
 Purpose   : Re-create an Err.pm object from a string containing Err data.
 Returns   : boolean, (was the Err.pm object rebuilt?)
 Argument  : message, note, tech passed from _initialize()
           : The message is examined to see if it contains a stringified error.

See Also   : L<_initialize>(), L<string>(), L<_has_err>()

=cut

#----------------------
sub _build_from_string {
#----------------------
    my ($self, $msg, $note, $tech) = @_;
    my @list = split "\n", $msg;
    my ($mode,$line);
    my $rebuilt = 0;

    # print "$ID: Attempting to build from string: $msg";<STDIN>;

    MEMBER:
    foreach $line (@list) {
	if($line =~ /^-+$/) { last MEMBER; }
	if($line =~ /^-+ (\w+) -+$/) { $self->{'_type'} = $1; $rebuilt = 1; next MEMBER; }
	if($line =~ /^MSG: *(\w.*)/) { my $msg = $1; 
				       if($self->_has_err($msg)) {
					   die "Duplicate error.";
				       }
				       $self->{'_msg'} = $msg; 
				       $mode = 'msg';
				       next MEMBER; }
	if($line =~ /^CONTEXT: *(\w.*)/) { push @{$self->{'_context'}}, $1; $mode = 'context'; next MEMBER; }
	if($line =~ /^NOTE: *(\w.*)/) { push @{$self->{'_note'}}, $1; $mode = 'note'; next MEMBER; }
	if($line =~ /^TECH: *(\w.*)/) { push @{$self->{'_tech'}}, $1; $mode = 'tech'; next MEMBER; }
	if($line =~ /^STACK:/) { $mode = 'stack'; next MEMBER; }
	next MEMBER if !$mode;
	SWITCH: {
	    local $_ = $mode;
	    m/msg/ && do{ $self->{'_msg'} .= "$line\n"; last SWITCH; };
	    m/note/ && do{ push @{$self->{'_note'}}, $line; last SWITCH; };
	    m/context/ && do{ push @{$self->{'_context'}}, $line; last SWITCH; };
	    m/tech/ && do{ push @{$self->{'_tech'}}, $line; last SWITCH; };
	    m/stack/ && do{ push @{$self->{'_stack'}}, $line; last SWITCH; };
	    next MEMBER;
	}
    }
    if($rebuilt) {
	## Optionally add additional notes.
	$self->_set_list_data('note', $note) if defined $note;
	$self->_set_list_data('tech', $tech) if defined $tech;
    }

    $rebuilt;
}


=head2 _has_err

 Usage     : n/a; internal method called by _build_from_string()
 Purpose   : Deterimine if an Err has already been set to prevent duplicate Errs.
 Returns   : boolean

See Also   : L<_build_from_string>()

=cut

#-------------
sub _has_err {
#-------------
    my ($self, $msg) = @_;

    $msg =~ s/^\s+//;
    $msg =~ s/\s+$//;

    my $err = $self->first;
    my ($existing_msg);
    do { 
#	print "checking err object $self\n";
	$existing_msg = $err->msg;
	$existing_msg =~ s/^\s+//;
	$existing_msg =~ s/\s+$//;
#	print "  msg: $existing_msg";<STDIN>;
	return 1 if $existing_msg eq $msg;

    } while($err = $err->next);
    
    0;
}


=head2 _set_type

 Usage     : n/a; internal method
 Purpose   : Sets the type of Err (warning, exception, fatal)
           : Called by _initialize()
 Argument  : string

=cut

#----------------
sub _set_type {
#----------------
    my( $self, $data ) = @_;
    $data ||= 'EXCEPTION';
    
#    printf "\n$ID: Setting type (%s) for err = %s\n", $data, $self->msg;<STDIN>;

    my (@type);
    if( @type = grep /$data/i, @Bio::Root::Err::ERR_TYPES ) {
	$self->{'_type'} = $type[0];
    } else {
	$self->{'_type'} = 'EXCEPTION';
    }

#    print "type = $self->{'_type'} for $self";<STDIN>;
}



=head2 _set_list_data

 Usage     : n/a; internal method used by set().
           : $err->_set_list_data( $member, $data);
 Purpose   : For data members which are anonymous arrays: note, tech, stack,
           : adds the given data to the list.
 Arguments : $member = any of qw(note tech stack)
           : $data   = string
 Comments  : Splits $data on tab. Each item 
           : of the split is a new entry.
           : To clobber the current data (unusual situation), you must first 
           : call set() with no data then call again with desired data.

=cut

#-------------------
sub _set_list_data { 
#-------------------
    my( $self, $member, $data ) = @_;
	
    # Sensitive to data member name changes.
    $member = "_\l$member";

#    $DEBUG && do {printf STDERR "\n$ID: Setting \"%s\" list data (%s)\n", $member, $data;<STDIN>; };

    defined $self->{$member} and return $self->_add_list_data( $member, $data );

    if( $data ) {
	$self->{$member} = [];
	if( $member =~ /stack/) {
	    foreach (@$data) { 
		push @{ $self->{$member}}, format_stack_entry(@$_);
	    }
	} else {
	    my @entries = split "\t", $data;
	    foreach (@entries) {
		next if /^$/;
#		$DEBUG && do {print STDERR "adding $member: $_";<STDIN>;};
		push @{ $self->{$member}}, $_;
	    }
	}
    } else {
	$self->{$member} = undef;
    }
}


=head2 _set_context

 Usage     : n/a; internal method used by set().
 Purpose   : Sets the object containment context for the exception.
           : (this is the hierarchy of objects in which the 
           :  exception occurred.)

=cut

#------------------
sub _set_context {
#------------------
    my($self, $aref) = @_;

    eval {
	if (!ref $aref) { 
#	    push @{$aref}, sprintf "object %s \"%s\"",ref($self->parent), $self->parent->name; 
	    push @{$aref}, "UNKNOWN CONTEXT";
	}
    };
    if($@) { push @{$aref}, 'undefined object'; }

    if($self->type eq 'EXCEPTION') {
	$aref->[0] = "Exception thrown by \l$aref->[0]";
    } else {
	$aref->[0] = "Error in \l$aref->[0]";
    }

    my $script = ($0 =~ /([\w\/\.]+)/, $1);
    push @$aref, "SCRIPT: $script";

    $self->{'_context'} = $aref;

#    print "$ID: _set_context():\n";
#    foreach(@$aref) { print "  $_\n"; }
#    <STDIN>;
}



=head2 set

 Usage     : $err->set( $member, $data );
 Purpose   : General accessor for setting any Err.pm data member. 
 Example   : $err->set('note', 'this is an additional note.');
 Returns   : n/a
 Argument  : $member = string, any of qw(msg type note tech stack)
           : $data   = string
 Throws    : n/a
 Comments  : Note, tech, and stack items are appended to any existing
           : notes, tech notes, and stack.
           : There should be no need to mess with the stack.

=cut

#---------
sub set {
#---------
    my( $self, $member, $data ) = @_;
    
    local $_ = "\u$member";
    SWITCH: {
	/msg/i && do{ $self->{'_msg'} = (defined $data ? $data : 'Unknown error'); last SWITCH; };
	/type/i && do{ $self->_set_type( $data ); last SWITCH; };
	/note|tech|stack/i && do{ $self->_set_list_data( $member, $data); last SWITCH};
	warn "\n*** Invalid or unspecified Err data member: $member\n\n";
    }
}


=head2 msg

 Usage     : $message = $err->msg;
 Purpose   : Get the main message associated with the exception.
 Returns   : String
 Argument  : optional string to be used as a delimiter.

See Also   : L<get>(), L<string>()

=cut


#-------
sub msg   { my($self,$delimiter) = @_; $self->get('msg',$delimiter); }
#-------


=head2 type

 Usage     : $type = $err->type;
 Purpose   : Get the type of Err (warning, exception, fatal)
 Returns   : String
 Argument  : optional string to be used as a delimiter.

See Also   : L<get>(), L<string>()

=cut

#--------
sub type { my($self,$delimiter) = @_; $self->get('type',$delimiter); }
#--------


=head2 note

 Usage     : $note = $err->note;
           : $note = $err->note('<P>');
 Purpose   : Get any general note associated with the exception.
 Returns   : String
 Argument  : optional string to be used as a delimiter.

See Also   : L<get>(), L<string>()

=cut

#---------
sub note  { my($self,$delimiter) = @_; $self->get('note',$delimiter); }
#---------


=head2 tech

 Usage     : $tech = $err->tech;
           : $tech = $err->tech('<P>');
 Purpose   : Get any technical note associate with the exception.
 Returns   : String
 Argument  : optional string to be used as a delimiter.

See Also   : L<get>(), L<string>()

=cut

#------------
sub tech  { my($self,$delimiter) = @_; $self->get('tech',$delimiter); }
#------------



=head2 stack

 Usage     : $stack = $err->stack;
           : $stack = $err->stack('<P>');
 Purpose   : Get the call stack for the exception. 
 Returns   : String
 Argument  : optional string to be used as a delimiter.

See Also   : L<get>(), L<string>()

=cut

#----------
sub stack { my($self,$delimiter) = @_; $self->get('stack',$delimiter); }
#----------



=head2 context

 Usage     : $context = $err->context;
           : $context = $err->context('<P>');
 Purpose   : Get the containment context of the object which generated the exception.
 Returns   : String
 Argument  : optional string to be used as a delimiter.

See Also   : L<get>(), L<string>()

=cut

#------------
sub context { my($self,$delimiter) = @_; $self->get('context',$delimiter); }
#------------



=head2 get

 Usage     : $err->get($member, $delimiter);
 Purpose   : Get specific data from the Err.pm object.
 Returns   : String in scalar context.
           : Array in list context.
 Argument  : $member = any of qw(msg type note tech stack context) or combination.
           : $delimiter = optional string to be used as a delimiter 
           : between member data.

See Also   : L<string>(), L<msg>(), L<note>(), L<tech>(), L<type>(), L<context>(), L<stack>()

=cut

#---------
sub get {
#---------
    my( $self, $member, $delimiter ) = @_;
    
    my $outer_delim = $delimiter || "\n";
#   my $outer_delim = ($CGI ? "\n<P>" : $delimiter);  ## Subtle bug here. 

    my (@out);
    local $_ = $member;
    SWITCH: {
	/type/i  && do{ push (@out, $self->{'_type'},$outer_delim) };
#	/msg/i   && do{ print "getting msg";<STDIN>; push (@out, (defined $self->{'_msg'} ? $self->{'_msg'} : ''),$outer_delim); print "msg: @out<---";<STDIN>; };
	/msg/i   && do{ push (@out, (defined $self->{'_msg'} ? $self->{'_msg'} : ''),$outer_delim); };
	/note/i  && do{ push (@out, $self->_get_list_data('note', $delimiter ),$outer_delim) };
	/tech/i  && do{ push (@out, $self->_get_list_data('tech', $delimiter ),$outer_delim) };
	/stack/i && do{ push (@out, $self->_get_list_data('stack', $delimiter ),$outer_delim) };
	/context/i && do{ push (@out, $self->_get_list_data('context', $delimiter ),$outer_delim) };

	## CAN'T USE THE FOLLOWING FORM SINCE IT FAILS WHEN $member EQUALS 'msgnote'.
#	/note|tech|stack/ && do{ push @out, $self->_get_list_data( $_, $delimiter ); };

	 last SWITCH;
	$self->warn("Invalid or undefined Err data member ($member).");
    }
#    $DEBUG && do{ print STDERR "OUTER DELIM = $outer_delim \nOUT: \n  @out <---";<STDIN>;};
    wantarray ? @out : join('',@out);
}



=head2 _get_list_data

 Usage     : n/a; internal method used by get()
 Purpose   : Gets data for members which are list refs (note, tech, stack, context)
 Returns   : Array
 Argument  : ($member, $delimiter)

See Also   : L<get>()

=cut

#-------------------
sub _get_list_data {
#-------------------
    my( $self, $member, $delimiter ) = @_;
    $delimiter ||= "\t";
    # Sensitive to data member name changes.
    $member = "_\l$member";
    return if !defined $self->{$member};
    join( $delimiter, @{$self->{$member}} );
}



=head2 get_all

 Usage     : (same as get())
 Purpose   : Get specific data from all errors in an Err.pm object.
 Returns   : Array in list context.
           : String in scalar context.
 Argument  : (same as get())

See Also   : L<get>()

=cut

#------------
sub get_all {
#------------
    my( $self, $member, $delimiter ) = @_;
    
    if( $self->size() == 1) {
	return $self->get( $member, $delimiter);
    } else {
	my $err = $self;

	### Return data from multiple errors in a list.
	if(wantarray) {
	    my @out;
	    do{ push @out, $err->get( $member);
	    } while($err = $err->prev());
	    return @out;

	} else {
	    ### Return data from multiple errors in a string with each error's data
	    ### bracketed by a "Error #n\n" line and two delimiters.
	    my $out = '';
	    if($err->size() == 1) {
		$out = $err->get( $member, $delimiter);
	    } else {
		do{ #$out .= "Error #${\$err->rank()}$delimiter";
		    $out .= $err->get( $member, $delimiter);
		    $out .= $delimiter.$delimiter;
		} while($err = $err->prev());
	    }
	    return $out; 
	}
    }  
}

#####################################################################################
##                              INSTANCE METHODS                                   ##
#####################################################################################


=head2 _add_note

 Usage     : n/a; internal method called by _add_list_data()
 Purpose   : adds a new note.

See Also   : L<_add_list_data>()

=cut

#---------------
sub _add_note {
#---------------
    my( $self, $data ) = @_;
    
    if( defined $self->{'_note'} ) {
	push @{ $self->{'_note'}}, $data;
    } else {
	$self->_set_list_data('note', $data );
    }
}

#----------------------------------------------------------------------
=head2 _add_tech()

 Usage     : n/a; internal method called by _add_list_data()
 Purpose   : adds a new technical note.

See Also   : L<_add_list_data>()

=cut

#-------------
sub _add_tech {
#-------------
    my( $self, $data ) = @_;
    
    if( defined $self->{'_tech'} ) {
	push @{ $self->{'_tech'}}, $data;
    } else {
	$self->_set_list_data('Tech', $data );
    }
}


=head2 _add_list_data

 Usage     : n/a; called by _set_list_data()
 Purpose   : adds a new note or tech note.

See Also   : L<_set_list_data>()

=cut

#--------------------
sub _add_list_data {
#--------------------
    my( $self, $member, $data ) = @_;
    
    local $_ = $member;
    SWITCH: {
	/note/i  && do{ $self->_add_note( $data ); };
	/tech/i  && do{ $self->_add_tech( $data ); };
    }
}



=head2 print

 Usage     : $err->print;
 Purpose   : Prints Err data to STDOUT or a FileHandle.
 Returns   : Call to print
 Argument  : Named parameters for string()
 Comments  : Uses string() to get data.

See Also   : L<string>()

=cut

#-----------
sub print { 
#-----------
    my( $self, %param ) = @_;  
#    my $OUT = $self->parent->fh(); 
#    print $OUT $self->string(%param); 
    print $self->string(%param); 
}


=head2 string

 Usage     : $err->string( %named_parameters);
 Purpose   : Stringify the data contained in the Err.pm object.
 Example   : print STDERR $err->string;
 Returns   : String
 Argument  : Named parameters (optional) passed to
           : Bio::Root::IOManager::set_display().

See Also   : L<print>(), L<_build_from_string>(), B<Bio::Root::IOManager::set_display()>

=cut

#-----------
sub string {
#-----------
    my( $self, @param ) = @_;

    my %param = @param;
    $self->set_display( @param );
    my $show = $self->show;
    my $out  = $param{-BEEP} ? "\a" : '';

    my $err = $param{-CURRENT} ? $self->last : $self->first;

#    my $err1 = $err;
#    my $errL = $self->last;
#    print "\n\nERR 1: ${\$err1->msg}";
#    print "\nERR L: ${\$errL->msg}";<STDIN>;

    my $numerate = $err->size() >1;
    my $count = 0;
    my ($title);
    my $hasnote = defined $self->{'_note'};
    my $hastech = defined $self->{'_tech'};

    while (ref $err)  {
	$count++;
#	$out .= sprintf "\nERROR #%d:", $count;

	if(not $title = $err->{'_type'}) {
	    $err = $err->next();
	    next;
	}
	if( $numerate) {
	    ## The rank data is a bit screwy at present.
	    $out .= sprintf "\n%s %s %s\n", '-'x 20, $title,'-'x 20;
	} else {
	    $out .= sprintf "\n%s %s %s\n", '-'x20, $title,'-'x20;
	}
	$show =~ /msg|default/i   and $out .= "MSG: " . $err->msg("\n"); 
	$show =~ /note|default/i  and $hasnote and $out .= "NOTE: ".$err->note("\n"); 
	$show =~ /tech|default/i  and $hastech and $out .= "TECH: ".$err->tech("\n"); 
	$show =~ /context|default/i  and $out .= "CONTEXT: ".$err->context("\n");
	$show =~ /stack|default/i and $out .= "STACK: \n".$err->stack("\n"); 
	$out .= sprintf "%s%s%s\n",'-'x 20, '-'x (length($title)+2), '-'x 20;

#	print "$ID: string: cumulative err:\n$out\n";<STDIN>;

	$err = $err->next();
    }

    $out;
}



=head2 is_fatal

 Usage     : $err->is_fatal;
 Purpose   : Determine if the error is of type 'FATAL'
 Returns   : Boolean
 Status    : Experimental, Deprecated

=cut

#--------------
sub is_fatal { my $self = shift; $self->{'_type'} eq 'FATAL'; }
#--------------

#####################################################################################
##                                CLASS METHODS                                    ##
#####################################################################################


=head2 throw

 Usage     : throw($object, [message], [note], [technical note]);
           : This method is exported.
 Purpose   : Class method version of Bio::Root::Object::throw().
 Returns   : die()s with the contents of the Err object in a string.
           : If the global strictness is less than -1, die is not called and 
           : the error is printed to STDERR.
 Argument  : [0] = object throwing the error.
           : [1] = optional message about the error.
           : [2] = optional note about the error.
           : [3] = optional technical note about the error.
 Comments  : The glogal verbosity level is not used. For verbosity-sensitive
           : behavior, use Bio::Root::Object::throw().
 Status    : Experimental
           : This method is an alternative to Bio::Root::Object::throw()
           : and is not as well developed or documented as that method.

See Also   : L<warning>(), B<Bio::Root::Object::throw()>, B<Bio::Root::Global::strictness()

=cut

#----------
sub throw {
#----------
    my($obj, @param) = @_;

#    print "Throwing exception for object ${\ref $self} \"${\$self->name}\"\n";
    my $err = new Bio::Root::Err(
			   -MSG     =>$param[0], 
			   -NOTE    =>$param[1],  
			   -TECH    =>$param[2], 
			   -STACK   =>scalar(Bio::Root::Object::stack_trace($obj,2)),
			   -CONTEXT =>Bio::Root::Object::containment($obj),
			   -TYPE    =>'EXCEPTION',
		#	    -PARENT =>$obj,
			   );

    if(strictness() < -1) {
	print STDERR $err->string(-BEEP=>1) unless verbosity() < 0;
    } else {
	die $err->string;
    }

    0;
}


=head2 warning

 Usage     : warning($object, [message], [note], [technical note]);
           : This method is exported.
 Purpose   : Class method version of Bio::Root::Object::warn().
 Returns   : Prints the contents of the error to STDERR and returns false (0).
           : If the global strictness() is > 1, warn() calls are converted 
           : into throw() calls.
 Argument  : [0] = object producing the warning.
           : [1] = optional message about the error.
           : [2] = optional note about the error.
           : [3] = optional technical note about the error.
           :
 Comments  : The glogal verbosity level is not used. For verbosity-sensitive
           : behavior, use Bio::Root::Object::warn().
 Status    : Experimental
           : This method is an alternative to Bio::Root::Object::warn()
           : and is not as well developed or documented as that method.

See Also   : L<throw>, B<Bio::Root::Object::warn()>, B<Bio::Root::Global::strictness()>

=cut

#-----------
sub warning {
#-----------
    my($obj, @param) = @_;

#    print "Throwing exception for object ${\ref $self} \"${\$self->name}\"\n";
    my $err = new Bio::Root::Err(
			   -MSG     =>$param[0], 
			   -NOTE    =>$param[1],  
			   -TECH    =>$param[2], 
			   -STACK   =>scalar(Bio::Root::Object::stack_trace($obj,2)),
			   -CONTEXT =>Bio::Root::Object::containment($obj),
			   -TYPE    =>'WARNING',
			#   -PARENT =>$obj,
			   );

    if(strictness() > 1) {
	die $err->string;

    } else {
	print STDERR $err->string(-BEEP=>1) unless $DONT_WARN;
    }

    0;
}


=head2 format_stack_entry

 Usage     : &format_stack_entry(<class>,<file>,<line>,<class_method>,<has_args>,<wantarray>)
           : This function is exported.
 Purpose   : Creates a single stack entry given a caller() list.
 Argument  : List of scalars (output of the caller() method).
 Returns   : String = class_method($line)
           : e.g., Bio::Root::Object::name(1234)

=cut

#------------------------
sub format_stack_entry {  
#------------------------
    my( $class, $file, $line, $classmethod, $hasargs, $wantarray) = @_;

#    if($DEBUG) {
#	print STDERR "format_stack_entry data:\n";
#	foreach(@_) {print STDERR "$_\n"; } <STDIN>;
#    }

    $classmethod ||= 'unknown class/method';
    $line        ||= 'unknown line';
    return "$classmethod($line)";
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

An instance of Bio::Root::Err.pm is a blessed reference to a hash containing
all or some of the following fields:

 FIELD     VALUE
 ------------------------------------------------------------------------
 _type     fatal | warning | exception (one of @Bio::Root::Err::ERR_TYPES).
        
 _msg      Terse description: Main cause of error. 
        
 _note     List reference. Verbose description: probable cause & troubleshooting for user.
        
 _tech     List reference. Technical notes of interest to programmer.
        
 _stack    List reference. Stack trace: list of "class::method(line number)" strings.
        


=cut

1;
