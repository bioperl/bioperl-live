#-----------------------------------------------------------------------------
# PACKAGE : Bio::Root::IOManager.pm
# AUTHOR  : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED : 26 Mar 1997
# REVISION: $Id$
# STATUS  : Alpha
#            
# For documentation, run this module through pod2html 
# (preferably from Perl v5.004 or better).
#
# MODIFIED: 
#   3 Feb 1999, sac:
#      * Added timeout support to read().
#      * Moved the FileHandle creation code out of read() and into 
#        Bio::Root::Utilties since it's of more general use.
#
#    24 Nov 1998, sac:
#      * Modified read(), compress(), and uncompress() to properly
#        deal with file ownership issues.
#
#    19 Aug 1998, sac:
#      * Fixed bug in display(), which wasn't returning true (1).
#
#    0.023, 20 Jul 1998, sac:
#      * read() can now use a supplied FileHandle or GLOB ref (\*IN).
#      * A few other touch-ups in read().
#
#    0.022, 16 Jun 1998, sac:
#      * read() now terminates reading when a supplied &$func_ref
#        returns false.
#
#    0.021, May 1998, sac:
#      * Refined documentation to use 5.004 pod2html.
#      * Properly using typglob refs as necessary 
#        (e.g., set_display(), set_fh()).
#
#    0.031, 2 Sep 1998, sac:
#      * Doc changes only
# 
# Copyright (c) 1997-8 Steve A. Chervitz. All Rights Reserved.
#           This module is free software; you can redistribute it and/or 
#           modify it under the same terms as Perl itself.
#-----------------------------------------------------------------------------

package Bio::Root::IOManager;

use Bio::Root::Global     qw(:devel $CGI);
use Bio::Root::Object     ();
use Bio::Root::Utilities  qw(:obj); 
use FileHandle            ();

@ISA   = qw(Bio::Root::Object);

use strict;
use vars qw($ID $VERSION $revision $Timeout);
$ID = 'Bio::Root::IOManager';
$VERSION = 0.042;
$Timeout = 20;   # Number of seconds to wait for input in read().

## POD Documentation:

=head1 NAME

Bio::Root::IOManager - Input and output manager for Perl5 objects.

=head1 SYNOPSIS

=head2 Object Creation

The creation of Bio::Root::IOManager.pm objects is handled by Bio::Root::Object.pm
which delegates various I/O tasks to this module.

    use Bio::Root::IOManager;
 
    $myIO = new Bio::Root::IOManager(-WHERE   =>'/usr/tmp/data.out', 
   				     -PARENT =>$self);


=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.


=head1 DESCRIPTION

This module encapsulates the data and methods necessary for regulating 
input/output (I/O) of data from Perl objects. 
It is concerned with "where" to get input or send output as opposed to "what" to get.
IOManager.pm is intended to consolidate various I/O issues for
Perl objects and provide an object-oriented way to do I/O things such as:

=over 4

=item * passing filehandles between objects,

=item * opening and reading input from files or STDIN,

=item * routine file management (compressing, uncompressing, and deleting).

=back

Subclasses of B<Bio::Root::Object.pm> have access to all methods defined in 
IOManager.pm since B<Bio::Root::Object.pm> employs Bio::Root::IOManager.pm 
by a delegation mechanism.

It is not clear yet how much objects really need to do the fancy I/O gymnastics as
supported by IOManager. Most of the time, objects simply send output to STDOUT
which is managed at the script/program level. The fancy I/O manipulations are
considered experimental and have not been adequately tested or utilized. 
I'm not really satisfied with the current L<display>()/L<set_display>() strategy.
The additional functionality is not often utilized in typical
applications. Is the extra complexity worth it?

B<The API for this module is under development.>


=head2 Generic Data Access & Manipulation

The L<read>() method provided permits the following:

=over 4

=item * read from a file or STDIN.

=item * read a single record or a stream containing multiple records.

=item * specify a record separator.

=item * store all input data in memory or process the data stream as it is being read.

=back

=head1 DEPENDENCIES

Bio::Root::IOManager.pm inherits from B<Bio::Root::Object.pm> and uses B<FileHandle.pm>.
B<Bio::Root::Utilities.pm> is also used for routine file manipulations 
compression/uncompression/deletion.

=head1 SEE ALSO

  Bio::Root::Object.pm       - Core object
  Bio::Root::Utilities.pm    - Generic utilty object
  Bio::Root::Global.pm       - Manages global variables/constants

  http://bio.perl.org/Projects/modules.html  - Online module documentation
  http://bio.perl.org/                       - Bioperl Project Homepage 
 
 FileHandle.pm (included in the Perl distribution or CPAN).

=head1 TODO

Experiment with using the newer B<IO.pm> included in the Perl distribution,
instead of FileHandle.pm.

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

Bio::Root::IOManager.pm, 0.041

=head1 ACKNOWLEDGEMENTS

This module was developed under the auspices of the Saccharomyces Genome
Database:
    http://genome-www.stanford.edu/Saccharomyces

=head1 COPYRIGHT

Copyright (c) 1997-98 Steve A. Chervitz. All Rights Reserved.
This module is free software; you can redistribute it and/or 
modify it under the same terms as Perl itself.

=cut

#
##
###
#### END of main POD documentation.
###
##
#'


=head1 APPENDIX

Methods beginning with a leading underscore are considered private
and are intended for internal use by this module. They are
B<not> considered part of the public interface and are described here
for documentation purposes only.

=cut



#####################################################################################
##                                 CONSTRUCTOR                                     ##
#####################################################################################


## Using default constructor and destructor inherited from Bio::Root::Object.pm

## Could perhaps set the file data member.


#####################################################################################
##                                 ACCESSORS                                       ##
#####################################################################################
 

=head2 file

 Usage     : $object->file([filename]);
 Purpose   : Set/Get name of a file associated with an object.
 Example   : $object->file('/usr/home/me/data.txt');
 Returns   : String (full path name)
 Argument  : String (full path name) (argument only required for setting)
 Throws    : Exception if the file appears to be empty or non-existent
 Comments  : File can be text or binary.

See Also   : L<compress_file>(), L<uncompress_file>(), L<delete_file>()

=cut

#--------
sub file { 
#--------
    my $self = shift; 
    if($_[0]) { 
	my $file = $_[0];
	if(not -s $file) {
	    $self->throw("File is empty or non-existent: $file");
	}	    
	$self->{'_file'} = $file;
    }
    $self->{'_file'};
}



=head2 set_fh

 Usage     : $self->set_fh( named_parameters )
 Purpose   : Sets various FileHandle data members ('fh', 'fherr').
           : Provides a public interface for _open_fh().
 Returns   : n/a
 Argument  : Named parameters:  (TAGS CAN BE UPPER OR LOWER CASE)
           :   -PATH  => string (filename) or a FileHandle object ref.
           :   -PRE   => string, prefix for opening (e.g., '>', '>>').
           :   -POST  => string, postfix for opening (e.g., '|'), for commands.
           :   -WHICH => string, 'err' for setting output path for errors.
           :
 Throws    : Exception propagated from _open_fh()
 Examples  : $self->set_fh();                   # Create anonymous FileHandle object
           : $self->set_fh(-PATH =>'fileName',  # Open for writing
	   :		   -PRE =>'>');         
           : $self->set_fh(-PATH =>'fileName',  # Open error log file in append mode.
	   :		   -PRE  =>'>>',
	   :		   -WHICH =>'err');  
           : $self->set_fh(-PATH =>$obj->fh()); # Copy a file handle from another object.
           :
 Comments  : set_read() and set_display() provide
           : interfaces for set_fh().
 Status    : Experimental

See also   : L<_open_fh>(), L<set_read>(), L<set_display>().

=cut

#-----------
sub set_fh { 
#-----------
    my( $self, %param) = @_;  

    no strict 'subs';
    my( $path, $prefix, $postfix, $which) = 
	$self->_rearrange([PATH,PRE,POST,WHICH],%param);  
    use strict 'subs';
    $prefix  ||= '';
    $postfix ||= '';
    $which   ||= '';
    my $fullpath = "$prefix$path$postfix";
    my($fh);

    $DEBUG and print STDERR "set_fh($fullpath) for ${\$self->name()}\n";

    if($which eq 'err') {
	if(ref($path) =~ /FileHandle|GLOB/ ) {
	    $fh = $path;
	} else {
	    if(defined $self->{'_fherr'}) { $self->_close_fh('err');}
	    if( not $fh = $self->_open_fh("$fullpath")) {
		$fh = $self->_open_fh("errors.$$");
		$fh || return;
		$self->warn("Couldn't set error output to $fullpath",
			    "Set to file errors.$$");
	    }
	}
	$self->{'_fherr_name'} = $fullpath; 
	$self->{'_fherr'} = $fh;

    } else {
	if(ref($path) =~ /FileHandle|GLOB/ ) {
	    $fh = $path;
	} else {
	    if(defined $self->{'_fh'}) { $self->_close_fh();}
	    if( not $fh = $self->_open_fh("$fullpath")) {
		$fh = $self->_open_fh("out.$$");
		$fh || return;
		$self->warn("Couldn't set output to $fullpath",
			    "Set to file out.$$");
	    }
	}
	$self->{'_fh_name'} = $fullpath; 
	$self->{'_fh'} = $fh;
	$DEBUG && print STDERR "$ID: set fh to: $fh";
    }
}



=head2 _open_fh

 Purpose   : Creates a new FileHandle object and returns it. 
           : This method can be used when you need to
           : pass FileHandles between objects.
 Returns   : The new FileHandle object.
 Throws    : Exception: if the call to new FileHandle fails.
 Examples  : $self->_open_fh();            # Create anonymous FileHandle object
           : $self->_open_fh('fileName');  # Open for reading
           : $self->_open_fh('>fileName'); # Open for writing
 Status    : Experimental

See also   : L<set_fh>(), L<fh>(), L<set_read>(), L<set_display>()

=cut

#-------------
sub _open_fh {
#-------------
    my( $self, $arg) = @_;
    my( $filehandle);

    $DEBUG and print STDERR "_open_fh() $arg\n";

    $filehandle = new FileHandle $arg;

#    if($arg =~ /STD[IO]/) {
#	$filehandle = new FileHandle;
#	$filehandle = *$arg;
#    } else {
#	 $filehandle = new FileHandle $arg;
#    }

    (ref $filehandle) || $self->throw("Can't create new FileHandle $arg",
				      "Cause: $!");
    return $filehandle;
}



=head2 _close_fh

 Purpose   : Destroy a FileHandle object.
 Returns   : n/a
 Status    : Experimental

See also   : L<_open_fh>(), L<set_fh>()

=cut

#--------------
sub _close_fh { 
#--------------
    my( $self, $arg) = @_; 
    $arg ||= '';
    if($arg eq 'err') { 
	close $self->{'_fherr'};
	undef $self->{'_fherr'}; 
    } else {
	close $self->{'_fh'};
	undef $self->{'_fh'}; 
    }
}	       


=head2 set_display

 Usage     : $self->set_display([-WHERE=>'path'],
	   :			[-SHOW =>'what is to be displayed'],
	   :			[-MODE =>'file open mode'])
 Purpose   : Sets a new FileHandle object for output.
           : - Sets the objects 'show' data member to 'default' if it is not defined.
           : - Is a wrapper for setting an object's STDOUT filehandle:
           :   Checks the -WHERE parameter and the status of the object's current
           :   filehandle {'_fh'} and does one of three things:
           :    1. If $param{-WHERE} is defined and is not 'STDOUT', it is sent to 
           :       set_fh() to open a new fh,
           :    2. else, if 'fh' has already been defined, it is returned,
           :    3. else, if where equals 'STDOUT', \*STDOUT is returned.
           :    4. else, \*STDOUT is returned.
           :
           : Thus, if an object has already set its 'fh' to some location,
           : it can still print to 'STDOUT' by explicitly passing -WHERE='STDOUT'
           : to display().
           :
 Arguments : Named parameters: (TAGS CAN BE UPPER OR LOWER CASE).
           : (all are optional).
           :    -WHERE => full path name of file to write to or 'STDOUT'.
           :    -SHOW  => what data is to be displayed. Becomes $self->{'_show'}
           :                     Default = 'default'. This results in a call to
           :                     _display_stats() method when display() is called
           :    -MODE  => mode for opening file. Default is overwrite '>'.
           :
 Returns   : FileHandle object reference or typglob reference (\*STDOUT).
 Throws    : Exception propagated from set_fh().
 Example   : $self->set_display();
           : $self->set_display(-WHERE=>'./data.out');
           : $self->set_display(-WHERE=>$obj->fh());
 Status    : Experimental
 Comments  : I'm not satisfied with the current display()/set_display() strategy.

See also   : L<display>(), L<set_fh>()

=cut

#----------------'
sub set_display {
#----------------
    my( $self, @param ) = @_;
    my ($show, $where, $mode) = $self->_rearrange([qw(SHOW WHERE MODE)], @param);

    ## Default mode: overwrite any existing file.
    $mode  ||= '>';
    $where ||= 'STDOUT';
    
    $self->{'_show'} = ($show || 'default');

    $DEBUG and print STDERR "$ID set_display() show: $self->{'_show'}\twhere: -->$where<--\n";

    if( defined $where and $where !~ /STDOUT/) {
#	print "setting file handle object\n";
	$self->set_fh(-PATH =>$where, 
		      -PRE  =>$mode);
    } elsif( not defined $self->{'_fh'} or $where =~ /STDOUT/)  {	    
	return \*STDOUT;
    } else  {
#	print STDERR "filehandle already set for this object: ${\$self->fh('name')}\n";
    }

    return $self->{'_fh'};
}



=head2 set_read

 Purpose   : Sets a new FileHandle object for input.
           : Same logic as set_display() but creates filehandle for read only.
 Returns   : The input FileHandle object or \*STDIN.
 Arguments : Named parameters: (TAGS CAN BE UPPER OR LOWER CASE).
           :    $param{-WHERE} = full path name of file to write to.
 Access    : Public
 Status    : Experimental, Deprecated
           :
 WARNING   : THIS METHOD HAS NOT BEEN TESTED AND IS LIKELY UNNECESSARY.
           : USE THE read() METHOD INSTEAD.
           :
           : Note also that set_read() uses the same data member as set_display()
           : so it is currently not possible to simultaneously have
           : different displaying and reading filehandles. This degree of
           : I/O control has not been necessary.

See also   : L<read>(), L<set_display>()

=cut

#-------------
sub set_read {
#-------------
    my( $self, @param ) = @_;
    my ($where, $mode) = $self->_rearrange([qw(WHERE MODE)], @param);

    ## Default mode: read only.
    $mode  ||= '<';
    $where ||= 'STDIN';
    
    if( ref($where) and $where !~ /STDIN/) {
#	print "setting file handle object\n";
	$self->set_fh(-PATH =>$where, 
		      -PRE  =>$mode);
    } elsif( not defined $self->{'_fh'} or $where =~ /STDIN/)  {	    
	return \*STDIN;
    } else  {
#	print STDERR "filehandle already set for this object: ${\$self->fh('name')}\n";
    }

    return $self->{'_fh'};
}



=head2 set_display_err

 Purpose   : Sets a new FileHandle object for outputing error information.
           : Same logic as set_display() but creates a filehandle in 
           : append mode.
 Returns   : The output FileHandle object for saving errors or \*STDERR.
 Status    : Experimental
 WARNING   : NOT TESTED

See also   : L<set_display>(), L<set_read>()

=cut

#--------------------
sub set_display_err {
#--------------------
    my( $self, @param ) = @_;
    my ($where, $mode) = $self->_rearrange([qw(WHERE MODE)], @param);

    ## Default mode: read only.
    $mode  ||= '>>';
    $where ||= 'STDERR';
    
    $DEBUG and print STDERR "set_display_err() object: ${\$self->name()}\n";

    if( ref($where) and $where !~ /STDERR/) {
#	print "setting file handle object\n";
	$self->set_fh(-PATH =>$where, 
		      -PRE  =>$mode);
    } elsif( not defined $self->{'_fherr'} or $where =~ /STDERR/)  {	    
	return \*STDERR;
    } else  {
#	print STDERR "filehandle already set for this object: ${\$self->fh('name')}\n";
    }

    return $self->{'_fherr'};
}


#####################################
#    GET ACCESSORS
#####################################


=head2 show

 Usage     : $self->show()
 Purpose   : Get the string used to specify what to display
           : using the display() method.
 Returns   : String or undef if no show data member is defined.
 Arguments : n/a

See also   : L<set_display>()

=cut

#----------
sub show { my $self= shift; $self->{'_show'}; }
#----------



=head2 fh

 Usage     : $object->fh(['name'])
 Purpose   : Accessor for an object's FileHandle object or the argument used
           : to create that object.
 Returns   : One of the following:
           :   1. The arguments used when the filehandle was created ('fh_name').
           :   2. The FileHandle object reference previously assigned to $self->{'_fh'}.
           :   3. Typeglob reference \*STDIN,  \*STDOUT or \*STDERR.
 Example   : $self->fh();          # returns filehandle for the STDIN/STDOUT path.
           : $self->fh('err');     # returns filehandle for the err file.
           : $self->fh('name');    # returns fh creation arguments.
           : $self->fh('errname'); # returns fh creation arguments for the err file.
 Status    : Experimental

See also   : L<set_display>(), L<set_read>(), L<set_fh>(), L<set_display_err>()

=cut

#--------'
sub fh { 
#--------
    my( $self, $type, $stream) = @_;
    $stream ||= 'out';
    $stream = ($stream eq 'in') ? \*STDIN : \*STDOUT;

    ## Problem: Without named parameters, how do you know if
    ## a single argument is to be assigned to $type or $stream?
    ## Function prototypes could be used, or separate methods:
    ## fh_out(), fh_in(), fh_err().
    $type or return ($self->{'_fh'} || $stream);

    if( $type =~ /name/){ 
	if($type =~ /err/ ) { return $self->{'_fherr_name'}; } 
	else                { return $self->{'_fh_name'}; }

    } else {
	if($type =~ /err/ ) { return ($self->{'_fherr'} || \*STDERR); } 
	else                { return ($self->{'_fh'}    || $stream); } 
    }
}


#####################################################################################
##                             INSTANCE METHODS                                    ##
#####################################################################################


##
##  INPUT METHODS:
##


=head2 read

 Usage     : $object->read(<named parameters>);
 Purpose   : Read raw textual data from a file or STDIN.
           : Optionally process each record it as it is read.
 Example   : $data = $object->read(-FILE    =>'usr/people/me/data.txt',
	   :			   -REC_SEP =>"\n:",
	   :			   -FUNC    =>\&process_rec);
           : $data = $object->read(-FILE  =>\*FILEHANDLE);
           : $data = $object->read(-FILE  =>new FileHandle $file, 'r');
           :
 Argument  : Named parameters: (TAGS CAN BE UPPER OR LOWER CASE)
           :  (all optional)
           :    -FILE    => string (full path to file) or a reference
           :                to a FileHandle object or typeglob. This is an
           :                optional parameter (if not defined, STDIN is used).
           :    -REC_SEP => record separator to be used
           :                when reading in raw data. If none is supplied,
           :                the default record separator is used ($/).
           :                $/ is localized to this method but be careful if
           :                you do any additional file reading in functions
           :                called by this method (see the -FUNC parameter).
           :                Such methods will use the value of $/ set
           :                by read() (if a -RE_SEP is supplied).
           :    -FUNC    => reference to a function to be called for each
           :                record. The return value of this function is now checked:
           :                if false, the reading is terminated.
           :                Typically -FUNC supplies a closure.
           :    -HANDLE  => reference to a FileHandle object or a
           :                typeglob to be use for reading input. 
           :                The FileHandle object should be configured to
           :                read from a desired file before calling this
           :                method. If both -handle and -file are defined,
           :                -handle takes precedence.
           :                (The -HANDLE parameter is no longer necessary
           :                 since -FILE can now contain a FileHandle ref.)
           :    -WAIT    => integer (number of seconds to wait for input
           :                before timing out. Default = 20 seconds).
           :
 Returns   : string, array, or undef depending on the arguments.
           : If a function reference is supplied, this function will be
           : called using the contents of each record as it is read in.
           : If no function reference is supplied, the data are returned as a
           : string in scalar context or as a list in array context.
           : The data are not altered; blank lines are not removed. 
           :
 Throws    : Exception if no input is read from source.
           : Exception if no input is read within WAIT seconds.
           : Exception if FUNC is not a function reference.
           : Propagates any exceptions thrown by create_filehandle()
           :
 Comments  : Gets the file name from the current file data member.
           : If no file has been defined, this method will attempt to
           : read from STDIN.
           :
           : COMPRESSED FILES:
           :    read() will attempt to uncompress the file 
           : if it appears to be compressed (binary file test).
           : If the file is compressed and the user is not the owner of 
           : the file, a temporary file is created which is removed after
           : it is read. Compressed files are always left in their
           : original compressed state.
           :
           : If the raw data is to be returned, wantarray is used to
           : determine how the data are to be returned (list or string).
           :
           : Sets the file data member to be the supplied file name.
           : (if any is supplied).

           : The read() method is a fairly new implementation
           : and uses a different approach than display().
           : For example, set_read() is not used.

 Bugs      : The following error is generated by Perl's FileHandle.pm module
           : when using the -w switch. It can be ignored for now:
  "Close on unopened file <GEN0> at /tools/perl/5.003/lib/FileHandle.pm line 255."

See Also   : L<file>(), L<create_filehandle>()

=cut

#----------'
sub read {
#----------
    my($self, @param) = @_;
    my( $rec_sep, $func_ref, $wait ) =
	$self->_rearrange([qw( REC_SEP FUNC WAIT)], @param);

    my $fmt = (wantarray ? 'list' : 'string');
    $wait ||= $Timeout;  # seconds to wait before timing out.

    my $FH = $Util->create_filehandle( -client => $self, @param);

    # Set the record separator (if necessary) using dynamic scope.
    local $/ = $rec_sep if scalar $rec_sep;

    # Verify that we have a proper reference to a function.
    if($func_ref) {
	if(not ref($func_ref) =~ /CODE/) { 
	    $self->throw("Not a function reference: $func_ref, ${\ref $func_ref}");
	}
    }

    $DEBUG && printf STDERR "$ID: read(): rec_sep = %s; func = %s\n",$/, ($func_ref?'defined':'none');
    
    my($data, $lines);

    $SIG{ALRM} = sub { die "Timed out!"; };
 
    eval {
	 alarm($wait);
      READ_LOOP:
	while(<$FH>) {
	    # Default behavior: read all lines.
	    # If &$func_ref returns false, exit this while loop.
	    # Uncomment to skip lines with only white space or record separators
#	    next if m@^(\s*|$/*)$@; 
	    
	    $lines++;
	    my($result);
	    if($func_ref) {
		$result = &$func_ref($_) or last READ_LOOP;
#		print "$ID read(): RESULT = $result\n"; 
	    } else {
		$data .= $_;
	    }
	}
	alarm(0);
    };
    if($@ =~ /Timed out!/) {
	 $self->throw("Timed out while waiting for input from $self->{'_input_type'}.");
    } elsif($@ =~ /\S/) {
         my $err = $@;
	 $self->throw("Unexpected error during read: $err");
    }

    close ($FH) unless $self->{'_input_type'} eq 'STDIN';
    
    # If the file was compressed and the user is the owner of file,
    #   then leave the file in its original compressed state.
    # If the file was compressed and the user was NOT the owner of file,
    #   then removed the compressed file which is a tmp file.

    $self->{'_compressed_file'} and ($self->{'_file_owner'} ? $self->compress_file() : $self->delete_file());

    if($data) {
	$DEBUG && do{ 
	    print STDERR "$ID: $lines records read.\nReturning $fmt.\n" };

	return ($fmt eq 'list') ? split("$/", $data) : $data;

    } elsif(not $func_ref) {
	$self->throw("No data input from $self->{'_input_type'}");
    }
    delete $self->{'_input_type'};
    delete $self->{'_file_owner'};
    delete $self->{'_compressed_file'};
    undef;
}


##
##  OUTPUT METHODS:
##


=head2 display

 Usage     : $self->set_display(named parameters)
 Purpose   : Provides a default display method which calls set_display()
           : and also invokes methods to display an object's stats
           : if necessary ( _print_stats_header() and _displayStats() ).
 Returns   : True (1).
 Throws    : Propagates any exceptions thrown by set_display().
 Arguments : Named parameters for set_display().
 Comments  : I'm not satisfied with the current display()/set_display() strategy.

See also   : L<set_display>(), L<_print_stats_header>()

=cut

#-------------
sub display {
#-------------
    my( $self, %param ) = @_;
    
    $DEBUG && print STDERR "$ID display for ${\ref($self)}\n";

    my $OUT = $self->set_display(%param);
#    my $OUT = $self->set_display( %param );
#    print "$ID: OUT = $OUT";<STDIN>;
    
    $DEBUG && do{ print STDERR "display(): WHERE = $OUT;\nSHOW = $self->{'_show'}";<STDIN>;};

    if($self->{'_show'} =~ /stats|default/i) {
	if($param{-HEADER}) {
	    $self->_print_stats_header($OUT); 
	}
	$self->parent->_display_stats($OUT); 
    }
    1;
}



=head2 _print_stats_header

 Usage     : n/a; internal method.
           : $obj->_print_stats_header(filehandle);
 Purpose   : Prints a header containing basic info about the object 
           : such as the class and name of the object followed by a 
           : line of hyphens.
 Status    : Experimental

=cut

#------------------------
sub _print_stats_header {
#------------------------
    my($self, $OUT) = @_;

    printf $OUT "\nSTATS FOR %s \"%s\"\n",ref($self->parent),$self->parent->name();
    printf $OUT "%s\n", '-'x60;
}




##
##  FILE MANIPULATION METHODS:
##



=head2 file_date

 Usage     : $object->file_date( %named_parameters);
 Purpose   : Get the last modified date of a file.
 Example   : $object->file_date();
           : $object->file_date(-FMT =>'yyyy-mmm-dd',
				-FILE =>'/usr/people/me/data.txt');
           : $object->file_date(-FMT =>'yyyy-mmm-dd');
 Returns   : String (date)
 Argument  : Named parameters:  (TAGS CAN BE UPPER OR LOWER CASE)
           :   -FILE  => string (filename full path)
           :   -FMT   => string (format for the returned date string)
           :
 Throws    : Exception if no file is specified or the file is non-existent
           : (Propagated from Utilities::file_date())
 Comments  : File can be text or binary.

See Also   : L<file>(), B<Bio::Root::Utilities::file_date()>

=cut

#---------------
sub file_date { 
#---------------
    my ($self, @param) = @_;
    my ($file, $fmt) = $self->_rearrange([qw(FILE FMT)], @param);
    
    if(not $file ||= $self->{'_file'}) {
	$self->throw("Can't get file date: no file specified");
    }
    $fmt ||= '';
    $Util->file_date($file, $fmt);  
}



=head2 compress_file

 Usage     : $object->compress_file([filename]);
 Purpose   : Compresses a file if not already compressed.
           : Compresses to a temorary file if user is not owner of supplied file.
 Example   : $object->file('/usr/home/me/data.txt');
           : $object->compress_file();   
 Argument  : String (full path name) (optional).
           : If no argument is provided, the file data member is used.
 Returns   : String (compressed file name, full path).
           : Sets the file data member to the compressed name
           : when not operating on a file supplied as an argument.
           : Returns false (undef) if the file is already compressed
           : (binary test).
 Throws    : Exception if no file is specified.
           : Propagates any exception thrown by Bio::Root::Utilities::compress()
           : if the file cannot be compressed().
           : Tests if file is already compressed to avoid trivial error due to 
           : the file already being compressed.
           :
 Comments  : Relies on the compress() method of Bio::Root::Utilities.pm
           : to implement the file compression functionality.
           : (Currently, Bio::Root::Utilities::compress() uses gzip.)
           :
           : If the user is not the owner of the file, the file is
           : compressed to a tmp file.
           :
           : All file compressing/uncompressing requests should go through
           : compress_file()/uncompress_file(). This serves to confine the
           : dependency between IOManager.pm module and Utilities.pm
           : which helps maintainability.
           :
 Bugs      : Only compresses text files. This obviates a dependency on 
           : particular file suffixes but is not good if you
           : want to compress a binary file.
           :
           : May not be taint-safe.   

See Also   : L<uncompress_file>(), L<file>(), B<Bio::Root::Utilities::compress()>

=cut

#-----------------
sub compress_file {
#-----------------
    my ($self, $file) = @_;
    my $myfile = 0;

    if(!$file) {
	$file = $self->{'_file'};
	$myfile = 1;
    }
    
    $file or $self->throw("Can't compress data file: no file specified");

    #printf STDERR "$ID: Compressing data file for %s\n  $file\n",$self->name();
    
    my ($newfile);
    if (-T $file) {
	$newfile = -o $file ? $Util->compress($file) : $Util->compress($file, 1);
	# set the current file to the new name.
	$self->file($newfile) if $myfile;
    }
    $newfile;
}



=head2 uncompress_file

 Usage     : $object->uncompress_file([filename]);
 Purpose   : Uncompresses the file containing the raw report.
           : Uncompresses to a temorary file if user is not owner of supplied file.
 Example   : $object->file('/usr/home/me/data.txt.gz');
           : $object->uncompress_file();   
 Argument  : String (full path name) (optional).
           : If no argument is provided, the file data member is used.
 Returns   : String (uncompressed file name, full path).
           : Sets the file data member to the uncompressed name
           : when not operating on a file supplied as an argument.
           : Returns false (undef) if the file is already uncompressed.
           :
 Throws    : Exception if no file is specified.
           : Propagates any exception thrown by Bio::Root::Utilities::compress()
           : if the file cannot be uncompressed().
           : Tests if file is already uncompressed to avoid trivial error due to 
           : the file already being uncompressed.
 Comments  : See comments for compress_file(). They apply here as well.
           : 
 Bugs      : Considers all binary files to be compressed. This obviates 
           : a dependency on particular file suffixes.
           : May not be taint safe.

See Also   : L<compress_file>(), L<file>(), B<Bio::Root::Utilities::uncompress()>

=cut

#--------------------
sub uncompress_file {
#--------------------
    my ($self, $file) = @_;
    my $myfile = 0;

    if(!$file) {
	$file = $self->{'_file'};
	$myfile = 1;
    }
    
    $file or $self->throw("Can't compress file: no file specified");

    #printf STDERR "$ID: Uncompressing data file for %s\n  $file",$self->name();

    my ($newfile);
    if (-B $file) {
	$newfile = -o $file ? $Util->uncompress($file) : $Util->uncompress($file, 1);
	# set the current file to the new name & return it.
	$self->file($newfile) if $myfile; 
    }
    $newfile;
}


=head2 delete_file

 Usage     : $object->delete_file([filename]);
 Purpose   : Delete a file.
 Example   : $object->delete_file('/usr/people/me/data.txt');
 Returns   : String (name of file which was deleted) if successful,
           : undef if file does not exist.
           : Sets the file data member to undef 
           : when not operating on a file supplied as an argument.
 Argument  : String (full path name) (optional).
           : If no argument is provided, the file data member is used.
 Throws    : Exception if the user is not the owner of the file.
           : Propagates any exception thrown by Bio::Root::Utilities::delete().
           : if the file cannot be deleted.
 Comments  : Be careful with this method: there is no undelete().
           : Relies on the delete() method provided by Bio::Root::Utilities.pm
           : to implement the file deletion functionality.
           : This method is not taint-safe.   
           : It is intended for off-line maintenance use only. 

See Also   : L<file>(), B<Bio::Root::Utilities::delete()>

=cut

#-----------------
sub delete_file {
#-----------------
    my ($self, $file) = @_;
    my $myfile = 0;

    if(!$file) {
	$file = $self->{'_file'};
	$myfile = 1;
    }
    return undef unless -e $file;

    -o $file or
	$self->throw("Can't delete file $file: Not owner."); 

#    $DEBUG and print STDERR "$ID: Deleting data file for ",$self->name();
    
    eval{ $Util->delete($file); };

    if(!$@ and $myfile) {
	$self->{'_file'} = undef;
    }
    $file;
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

An instance of Bio::Root::IOManager.pm is a blessed reference to a hash containing
all or some of the following fields:

 FIELD          VALUE
 ------------------------------------------------------------------------
  _show         Selects display options.

  _fh           FileHandle object for redirecting STDIN or STDOUT.

  _fherr        FileHandle object for error data. Append mode.

  _fh_name      The arguments used to create fh.

  _fherr_name   The arguments used to create fherr.
  
  INHERITED DATA MEMBERS
  
  _parent       (From Bio::Root::Object.pm> Object reference for the owner of this IOManager.
 
=cut

1;
