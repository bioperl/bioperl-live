# $Id$
#
# BioPerl module for Bio::Root::IO
#
# Cared for by Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Root::IO - module providing several methods often needed when dealing with file IO

=head1 SYNOPSIS

    # utilize stream I/O in your module
    $self->{'io'} = Bio::Root::IO->new(-file => "myfile");
    $self->{'io'}->_print("some stuff");
    $line = $self->{'io'}->_readline();
    $self->{'io'}->_pushback($line);
    $self->{'io'}->close();

    # obtain platform-compatible filenames
    $path = Bio::Root::IO->catfile($dir, $subdir, $filename);
    # obtain a temporary file (created in $TEMPDIR)
    ($handle) = $io->tempfile();

=head1 DESCRIPTION

This module provides methods that will usually be needed for any sort
of file- or stream-related input/output, e.g., keeping track of a file
handle, transient printing and reading from the file handle, a close
method, automatically closing the handle on garbage collection, etc.

To use this for your own code you will either want to inherit from
this module, or instantiate an object for every file or stream you are
dealing with. In the first case this module will most likely not be
the first class off which your class inherits; therefore you need to
call _initialize_io() with the named parameters in order set file
handle, open file, etc automatically.

Most methods start with an underscore, indicating they are private. In
OO speak, they are not private but protected, that is, use them in
your module code, but a client code of your module will usually not
want to call them (except those not starting with an underscore).

In addition this module contains a couple of convenience methods for
cross-platform safe tempfile creation and similar tasks. There are
some CPAN modules related that may be not be available on all
platforms. At present, File::Spec and File::Temp are attempted. This
module exports $TEMPFILE and $ROOTDIR, which will always be set,
$PATHSEP, which will be set if File::Spec fails, and $OPENFLAGS, which
will be set if either of File::Spec or File::Temp fails.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Hilmar Lapp

Email hlapp@gmx.net

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Root::IO;
use vars qw(@ISA $FILESPECLOADED $FILETEMPLOADED $FILEPATHLOADED
	    $TEMPDIR $PATHSEP $ROOTDIR $OPENFLAGS $VERBOSE);
use strict;

use Symbol;
use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);

my $TEMPCOUNTER;

BEGIN {
    $TEMPCOUNTER = 0;
    $FILESPECLOADED = 0;
    $FILETEMPLOADED = 0;
    $FILEPATHLOADED = 0;
    $VERBOSE = 1;

    # try to load those modules that may cause trouble on some systems
    eval { 
	require File::Path;
	$FILEPATHLOADED = 1;
    }; 
    if( $@ ) {
	print STDERR "Cannot load File::Path: $@" if( $VERBOSE > 0 );
	# do nothing
    } 
    eval {
	require File::Spec;
	$FILESPECLOADED = 1;
	$TEMPDIR = File::Spec->tmpdir();
	$ROOTDIR = File::Spec->rootdir();
	require File::Temp; # tempfile creation
	$FILETEMPLOADED = 1;
    };
    if( $@ ) { 
	if(! defined($TEMPDIR)) { # File::Spec failed
	    # determine tempdir
	    if (defined $ENV{'TEMPDIR'} && -d $ENV{'TEMPDIR'} ) {
		$TEMPDIR = $ENV{'TEMPDIR'};
	    } elsif( defined $ENV{'TMPDIR'} && -d $ENV{'TMPDIR'} ) {
		$TEMPDIR = $ENV{'TMPDIR'};
	    }
	    if($^O =~ /mswin/i) {
		$TEMPDIR = 'C:\TEMP' unless $TEMPDIR;
		$PATHSEP = "\\";
		$ROOTDIR = 'C:';
	    } elsif($^O =~ /macos/i) {
		$TEMPDIR = "" unless $TEMPDIR; # what is a reasonable default on Macs?
		$PATHSEP = ":";
		$ROOTDIR = ""; # what is the reasonable
	    } else { # unix
		$TEMPDIR = "/tmp" unless $TEMPDIR;
		$PATHSEP = "/";
		$ROOTDIR = "/";
	    }
	    if (!( -d $TEMPDIR && -w $TEMPDIR )) {
		$TEMPDIR = '.'; # last resort
	    }
	}
	# File::Temp failed (alone, or File::Spec already failed)
	#
	# determine open flags for tempfile creation -- we'll have to do this
	# ourselves
	use Fcntl;
	use Symbol;
	$OPENFLAGS = O_CREAT | O_EXCL | O_RDWR;
	for my $oflag (qw/FOLLOW BINARY LARGEFILE EXLOCK NOINHERIT TEMPORARY/){
	    my ($bit, $func) = (0, "Fcntl::O_" . $oflag);
	    no strict 'refs';
	    $OPENFLAGS |= $bit if eval { $bit = &$func(); 1 };
	}
    }
}

=head2 new

 Title   : new 
 Usage   : 
 Function: Overridden here to automatically call _initialize_io().
 Example :
 Returns : new instance of this class
 Args    : named parameters


=cut

sub new {
    my ($caller, @args) = @_;
    my $self = $caller->SUPER::new(@args);

    $self->_initialize_io(@args);
    return $self;
}

=head2 _initialize_io

 Title   : initialize_io
 Usage   : $self->_initialize_io(@params);
 Function: Initializes filehandle and other properties.

           Currently recognized the following named parameters: -file, -fh
 Example :
 Returns : TRUE
 Args    : named parameters


=cut

sub _initialize_io {
    my($self, @args) = @_;

    $self->_register_for_cleanup(\&_io_cleanup);

    my ($input, $file, $fh) = $self->_rearrange([qw(INPUT FILE FH)], @args);

    delete $self->{'_readbuffer'};
    delete $self->{'_filehandle'};

    # determine whether the input is a file(name) or a stream
    if($input) {
	if(ref(\$input) eq "SCALAR") {
	    # we assume that a scalar is a filename
	    if($file && ($file ne $input)) {
		$self->throw("input file given twice: $file and $input disagree");
	    }
	    $file = $input;
	} elsif(ref($input) &&
		((ref($input) eq "GLOB") || $input->isa('IO::Handle'))) {
	    # input is a stream
	    $fh = $input;
	} else {
	    # let's be strict for now
	    $self->throw("unable to determine type of input $input: ".
			 "not string and not GLOB");
	}
    }
    if(defined($file) && defined($fh)) {
	$self->throw("Providing both a file and a filehandle for reading - only one please!");
    }

    if(defined($file) && ($file ne '')) {
	$fh = Symbol::gensym();
	open ($fh,$file) ||
	    $self->throw("Could not open $file for reading: $!");
    }
    $self->_fh($fh) if $fh; # if not provided, defaults to STDIN and STDOUT
    return 1;
}

=head2 _fh

 Title   : _fh
 Usage   : $obj->_fh($newval)
 Function:
 Example :
 Returns : value of _filehandle
 Args    : newvalue (optional)


=cut

sub _fh {
    my ($obj, $value) = @_;
    if ( defined $value) {
	$obj->{'_filehandle'} = $value;
    }
    return $obj->{'_filehandle'};
}

=head2 _print

 Title   : _print
 Usage   : $obj->_print(@lines)
 Function:
 Example :
 Returns : writes output

=cut

sub _print {
    my $self = shift;
    my $fh = $self->_fh || \*STDOUT;
    print $fh @_;
}

=head2 _readline

 Title   : _readline
 Usage   : $obj->_readline
 Function: Reads a line of input.

           Note that this method implicitely uses the value of $/ that is
           in effect when called.

           Note also that the current implementation does not handle pushed
           back input correctly unless the pushed back input ends with the
           value of $/.
 Example :
 Returns : 

=cut

sub _readline {
    my $self = shift;
    my $fh = $self->_fh || \*STDIN;
    my $line;
    
    # if the buffer been filled by _pushback then return the buffer
    # contents, rather than read from the filehandle
    if(exists($self->{'_readbuffer'})) {
	$line = $self->{'_readbuffer'};
	delete $self->{'_readbuffer'};	
    } else {
	$line = <$fh>;
    }
    $line =~ s/\r\n/\n/g if (defined $line);
    return $line;
}

=head2 _pushback

 Title   : _pushback
 Usage   : $obj->_pushback($newvalue)
 Function: puts a line previously read with _readline back into a buffer
 Example :
 Returns :
 Args    : newvalue

=cut

sub _pushback {
    my ($obj, $value) = @_;
    $value .= $obj->{'_readbuffer'} if(exists($obj->{'_readbuffer'}));
    $obj->{'_readbuffer'} = $value;
}

=head2 close

 Title   : close
 Usage   : $seqio->close()
 Function: Closes the file handle associated with this seqio system
 Example :
 Returns :
 Args    :

=cut

sub close {
   my ($self, @args) = @_;

   $self->{'_filehandle'} = undef;
   delete $self->{'_readbuffer'};
}

sub _io_cleanup {
    my ($self) = @_;

    $self->close();

    # we are planning to cleanup temp files no matter what    
    if( exists($self->{'_rootio_tempfiles'}) &&
	ref($self->{'_rootio_tempfiles'}) =~ /array/i) { 
	if( $self->verbose > 0 ) {
	    print STDERR "going to remove files ", 
	    join(",",  @{$self->{'_rootio_tempfiles'}}), "\n";
	}
	unlink  (@{$self->{'_rootio_tempfiles'}} );
    }
    # cleanup if we are not using File::Temp
    if( $self->{'_cleanuptempdir'} &&
	exists($self->{'_rootio_tempdirs'}) &&
	ref($self->{'_rootio_tempdirs'}) =~ /array/i) {	

	if( $self->verbose > 0 ) {
	    print STDERR "going to remove dirs ", 
	    join(",",  @{$self->{'_rootio_tempdirs'}}), "\n";
	}
	$self->rmtree( $self->{'_rootio_tempdirs'});
    }
}


=head2 tempfile

 Title   : tempfile
 Usage   : my ($handle,$tempfile) = $io->tempfile(); 
 Function: Returns a temporary filename and a handle opened for writing and
           and reading.

 Caveats : If you do not have File::Temp on your system you should avoid
           specifying TEMPLATE and SUFFIX. (We don't want to recode
           everything, okay?)
 Returns : a 2-element array, consisting of temporary handle and temporary 
           file name
 Args    : named parameters compatible with File::Temp: DIR (defaults to
           $Bio::Root::IO::TEMPDIR), TEMPLATE, SUFFIX.

=cut

#'
sub tempfile {
    my ($self, @args) = @_;
    my ($tfh, $file);
    my %params = @args;

    # map between naming with and without dash
    foreach my $key (grep { $_ =~ /^-/; } keys(%params)) {
	$params{substr($key,1)} = $params{$key};
	delete $params{$key};
    }

    $params{'DIR'} = $TEMPDIR if(! exists($params{'DIR'}));
    if($FILETEMPLOADED) {
	if(exists($params{'TEMPLATE'})) {
	    my $template = $params{'TEMPLATE'};
	    delete $params{'TEMPLATE'};
	    ($tfh, $file) = File::Temp::tempfile($template, %params);
	} else {
	    ($tfh, $file) = File::Temp::tempfile(@args);
	}
    } else {
	my $dir = $params{'DIR'};
	$file = $self->catfile($dir,
			       (exists($params{'TEMPLATE'}) ?
				$params{'TEMPLATE'} :
				sprintf( "%s-%s-%s",  
					 $ENV{USER} || 'unknown', $$, 
					 $TEMPCOUNTER++)));
	# taken from File::Temp
	if ($] < 5.006) {
	    $tfh = &Symbol::gensym;
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
	if ( sysopen($tfh, $file, $OPENFLAGS, 0600) ) {
	    # Reset umask
	    umask($umask); 
	} else { 
	    $self->throw("Could not open tempfile $file: $!\n");
	}
    }
    push @{$self->{'_rootio_tempfiles'}}, $file;
    return ($tfh,$file);
}

=head2  tempdir

 Title   : tempdir
 Usage   : my ($tempdir) = $io->tempdir(CLEANUP=>1); 
 Function: Creates and returns the name of a new temporary directory.

           Note that you should not use this function for obtaining "the"
           temp directory. Use $Bio::Root::IO::TEMPDIR for that. Calling this
           method will in fact create a new directory.

 Returns : The name of a new temporary directory.
 Args    : args - ( key CLEANUP ) indicates whether or not to cleanup 
           dir on object destruction, other keys as specified by File::Temp

=cut

sub tempdir {
    my ( $self, @args ) = @_;
    if($FILETEMPLOADED && File::Temp->can('tempdir') ) {
	return File::Temp::tempdir(@args);
    }

    # we have to do this ourselves, not good
    #
    # we are planning to cleanup temp files no matter what
    my %params = @args;
    $self->{'_cleanuptempdir'} = $params{CLEANUP} == 1;
    my $tdir = $self->catfile($TEMPDIR,
			      sprintf("dir_%s-%s-%s", 
				      $ENV{USER} || 'unknown', $$, 
				      $TEMPCOUNTER++));
    mkdir($tdir, 0755);
    push @{$self->{'_rootio_tempdirs'}}, $tdir; 
    return $tdir;
}

=head2 catfile

 Title   : catfile
 Usage   : $path = Bio::Root::IO->catfile(@dirs,$filename);
 Function: Constructs a full pathname in a cross-platform safe way.

           If File::Spec exists on your system, this routine will merely
           delegate to it. Otherwise it tries to make a good guess.

           You should use this method whenever you construct a path name
           from directory and filename. Otherwise you risk cross-platform
           compatibility of your code.

           You can call this method both as a class and an instance method.

 Returns : a string
 Args    : components of the pathname (directories and filename, NOT an
           extension)

=cut

sub catfile {
    my ($self, @args) = @_;

    return File::Spec->catfile(@args) if($FILESPECLOADED);
    # this is clumsy and not very appealing, but how do we specify the
    # root directory?
    if($args[0] eq '/') {
	$args[0] = $ROOTDIR;
    }
    return join($PATHSEP, @args);
}

=head2 rmtree

 Title   : rmtree
 Usage   : Bio::Root::IO->rmtree($dirname );
 Function: Remove a full directory tree

           If File::Path exists on your system, this routine will merely
           delegate to it. Otherwise it runs a local version of that code.

           You should use this method to remove directories which contain 
           files.

           You can call this method both as a class and an instance method.

 Returns : number of files successfully deleted
 Args    : roots - rootdir to delete or reference to list of dirs

           verbose - a boolean value, which if TRUE will cause
                     C<rmtree> to print a message each time it
                     examines a file, giving the name of the file, and
                     indicating whether it's using C<rmdir> or
                     C<unlink> to remove it, or that it's skipping it.
                     (defaults to FALSE)

           safe - a boolean value, which if TRUE will cause C<rmtree>
                  to skip any files to which you do not have delete
                  access (if running under VMS) or write access (if
                  running under another OS).  This will change in the
                  future when a criterion for 'delete permission'
                  under OSs other than VMS is settled.  (defaults to
                  FALSE)

=cut

# taken straight from File::Path VERSION = "1.0403"
sub rmtree {
    my($self,$roots, $verbose, $safe) = @_;
    if( $FILEPATHLOADED ) { 
	return File::Path::rmtree ($roots, $verbose, $safe); 
    }				
    
    my $force_writeable = ($^O eq 'os2' || $^O eq 'dos' || $^O eq 'MSWin32'
		       || $^O eq 'amigaos');
    my $Is_VMS = $^O eq 'VMS';

    my(@files);
    my($count) = 0;
    $verbose ||= 0;
    $safe ||= 0;
    if ( defined($roots) && length($roots) ) {
	$roots = [$roots] unless ref $roots;
    } else {
	$self->warn("No root path(s) specified\n");
	return 0;
    }

    my($root);
    foreach $root (@{$roots}) {
	$root =~ s#/\z##;
	(undef, undef, my $rp) = lstat $root or next;
	$rp &= 07777;	# don't forget setuid, setgid, sticky bits
	if ( -d _ ) {
	    # notabene: 0777 is for making readable in the first place,
	    # it's also intended to change it to writable in case we have
	    # to recurse in which case we are better than rm -rf for 
	    # subtrees with strange permissions
	    chmod(0777, ($Is_VMS ? VMS::Filespec::fileify($root) : $root))
	      or $self->warn("Can't make directory $root read+writeable: $!")
		unless $safe;
	    if (opendir(DIR, $root) ){
		@files = readdir DIR;
		closedir(DIR);
	    } else {
	        $self->warn( "Can't read $root: $!");
		@files = ();
	    }

	    # Deleting large numbers of files from VMS Files-11 filesystems
	    # is faster if done in reverse ASCIIbetical order 
	    @files = reverse @files if $Is_VMS;
	    ($root = VMS::Filespec::unixify($root)) =~ s#\.dir\z## if $Is_VMS;
	    @files = map("$root/$_", grep $_!~/^\.{1,2}\z/s,@files);
	    $count += $self->rmtree([@files],$verbose,$safe);
	    if ($safe &&
		($Is_VMS ? !&VMS::Filespec::candelete($root) : !-w $root)) {
		print "skipped $root\n" if $verbose;
		next;
	    }
	    chmod 0777, $root
	      or $self->warn( "Can't make directory $root writeable: $!")
		if $force_writeable;
	    print "rmdir $root\n" if $verbose;
	    if (rmdir $root) {
		++$count;
	    }
	    else {
		$self->warn( "Can't remove directory $root: $!");
		chmod($rp, ($Is_VMS ? VMS::Filespec::fileify($root) : $root))
		    or $self->warn("and can't restore permissions to "
		            . sprintf("0%o",$rp) . "\n");
	    }
	}
	else {

	    if ($safe &&
		($Is_VMS ? !&VMS::Filespec::candelete($root)
		         : !(-l $root || -w $root)))
	    {
		print "skipped $root\n" if $verbose;
		next;
	    }
	    chmod 0666, $root
	      or $self->warn( "Can't make file $root writeable: $!")
		if $force_writeable;
	    print "unlink $root\n" if $verbose;
	    # delete all versions under VMS
	    for (;;) {
		unless (unlink $root) {
		    $self->warn( "Can't unlink file $root: $!");
		    if ($force_writeable) {
			chmod $rp, $root
			    or $self->warn("and can't restore permissions to "
			            . sprintf("0%o",$rp) . "\n");
		    }
		    last;
		}
		++$count;
		last unless $Is_VMS && lstat $root;
	    }
	}
    }

    $count;
}

1;
