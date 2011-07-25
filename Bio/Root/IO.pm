#
# BioPerl module for Bio::Root::IO
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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
call _initialize_io() with the named parameters in order to set file
handle, open file, etc automatically.

Most methods start with an underscore, indicating they are private. In
OO speak, they are not private but protected, that is, use them in
your module code, but a client code of your module will usually not
want to call them (except those not starting with an underscore).

In addition this module contains a couple of convenience methods for
cross-platform safe tempfile creation and similar tasks. There are
some CPAN modules related that may not be available on all
platforms. At present, File::Spec and File::Temp are attempted. This
module defines $PATHSEP, $TEMPDIR, and $ROOTDIR, which will always be set, 
and $OPENFLAGS, which will be set if either of File::Spec or File::Temp fails.

The -noclose boolean (accessed via the noclose method) prevents a
filehandle from being closed when the IO object is cleaned up.  This
is special behavior when a object like a parser might share a
filehandle with an object like an indexer where it is not proper to
close the filehandle as it will continue to be reused until the end of the
stream is reached.  In general you won't want to play with this flag.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR - Hilmar Lapp

Email hlapp@gmx.net

=head1 CONTRIBUTORS

Mark A. Jensen ( maj -at- fortinbras -dot- us )

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Root::IO;

our ($FILESPECLOADED,   $FILETEMPLOADED,
    $FILEPATHLOADED,    $TEMPDIR,
    $PATHSEP,           $ROOTDIR,
    $OPENFLAGS,         $VERBOSE,
    $ONMAC,             $HAS_LWP,
    $HAS_EOL);

use strict;

use Symbol;
use POSIX qw(dup);
use IO::Handle;
use Bio::Root::HTTPget;

use base qw(Bio::Root::Root);

my $TEMPCOUNTER;
my $HAS_WIN32 = 0;

BEGIN {
    $TEMPCOUNTER = 0;
    $FILESPECLOADED = 0;
    $FILETEMPLOADED = 0;
    $FILEPATHLOADED = 0;
    $VERBOSE = 0;

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
        require LWP::UserAgent;
    };
    if( $@ ) {
        print STDERR "Cannot load LWP::UserAgent: $@" if( $VERBOSE > 0 );
        $HAS_LWP = 0;
    } else {
        $HAS_LWP = 1;
    }

    # If on Win32, attempt to find Win32 package

    if($^O =~ /mswin/i) {
    eval {
        require Win32;
        $HAS_WIN32 = 1;
    };
    }

    # Try to provide a path separator. Why doesn't File::Spec export this,
    # or did I miss it?
    if($^O =~ /mswin/i) {
        $PATHSEP = "\\";
    } elsif($^O =~ /macos/i) {
        $PATHSEP = ":";
    } else { # unix
        $PATHSEP = "/";
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
            $ROOTDIR = 'C:';
        } elsif($^O =~ /macos/i) {
            $TEMPDIR = "" unless $TEMPDIR; # what is a reasonable default on Macs?
            $ROOTDIR = ""; # what is reasonable??
        } else { # unix
            $TEMPDIR = "/tmp" unless $TEMPDIR;
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
    $ONMAC = "\015" eq "\n";
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
 Function: Initializes filehandle and other properties from the parameters.

           Currently recognizes the following named parameters:
              -file     name of file to open
              -string   a string that is to be converted to a filehandle
              -url      name of URL to open
              -input    name of file, or GLOB, or IO::Handle object
              -fh       file handle (mutually exclusive with -file)
              -flush    boolean flag to autoflush after each write
              -noclose  boolean flag, when set to true will not close a
                        filehandle (must explicitly call close($io->_fh)
              -retries  number of times to try a web fetch before failure
                        
              -ua_parms hashref of key => value parameters to pass 
                        to LWP::UserAgent->new()
                        (only meaningful with -url is set)
                        A useful value might be, for example,
                        { timeout => 60 } (ua default is 180 sec)
 Returns : TRUE
 Args    : named parameters


=cut

sub _initialize_io {
    my($self, @args) = @_;

    $self->_register_for_cleanup(\&_io_cleanup);

    my ($input, $noclose, $file, $fh, $string, $flush, $url,
    $retries, $ua_parms) = 
    $self->_rearrange([qw(INPUT NOCLOSE FILE FH STRING FLUSH URL RETRIES UA_PARMS)],
                      @args);

    if($url){
        $retries ||= 5;
    
        if($HAS_LWP) { #use LWP::UserAgent
            require LWP::UserAgent;
            my $ua = LWP::UserAgent->new(%$ua_parms);
            my $http_result;
            my($handle,$tempfile) = $self->tempfile();
            CORE::close($handle);
          
    
            for(my $try = 1 ; $try <= $retries ; $try++){
                $http_result = $ua->get($url, ':content_file' => $tempfile);
                $self->warn("[$try/$retries] tried to fetch $url, but server ".
                            "threw ". $http_result->code . ".  retrying...")
                            if !$http_result->is_success;
                last if $http_result->is_success;
            }
            $self->throw("failed to fetch $url, server threw ".
                         $http_result->code) if !$http_result->is_success;
    
            $input = $tempfile;
            $file  = $tempfile;
        } else { #use Bio::Root::HTTPget
            #$self->warn("no lwp");
    
            $fh = Bio::Root::HTTPget::getFH($url);
        }
        }

    delete $self->{'_readbuffer'};
    delete $self->{'_filehandle'};
    $self->noclose( $noclose) if defined $noclose;
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
        $self->throw("Providing both a file and a filehandle for reading - ".
                     "only one please!");
    }

    if ($string) {
        if(defined($file) || defined($fh)) {
            $self->throw("File or filehandle provided with -string,".
                         " please unset if you are using -string as a file");
        }
        open($fh, "<", \$string)
    }
    
    if(defined($file) && ($file ne '')) {
        $fh = Symbol::gensym();
        open ($fh,$file) || $self->throw("Could not open $file: $!");
        $self->file($file);
    }

    if (defined $fh) {
        # check filehandle to ensure it's one of:
        # a GLOB reference, as in: open(my $fh, "myfile");
        # an IO::Handle or IO::String object
        # the UNIVERSAL::can added to fix Bug2863
        unless ( ( ref $fh && ( ref $fh eq 'GLOB' ) )
                 || ( ref $fh && ( UNIVERSAL::can( $fh, 'can' ) 
                    && ( $fh->isa('IO::Handle') || $fh->isa('IO::String') ) ) ) 
               ) {
            $self->throw("file handle $fh doesn't appear to be a handle");
        }
    }
    if ($HAS_EOL) {
        binmode $fh, ':raw:eol(LF-Native)';
    }
    $self->_fh($fh) if $fh; # if not provided, defaults to STDIN and STDOUT

    $self->_flush_on_write(defined $flush ? $flush : 1);

    return 1;
}

=head2 _fh

 Title   : _fh
 Usage   : $obj->_fh($newval)
 Function: Get/set the file handle for the stream encapsulated.
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

=head2 mode

 Title   : mode
 Usage   : $obj->mode()
 Function:
 Example :
 Returns : mode of filehandle:
           'r' for readable
           'w' for writable
           '?' if mode could not be determined
 Args    : -force (optional), see notes.
 Notes   : once mode() has been called, the filehandle's mode is cached
           for further calls to mode().  to override this behavior so
           that mode() re-checks the filehandle's mode, call with arg
           -force

=cut

sub mode {
    my ($obj, @arg) = @_;
    my %param = @arg;
    return $obj->{'_mode'} if defined $obj->{'_mode'} and !$param{-force};
    
    # Previous system of:
    #  my $iotest = new IO::Handle;
    #  $iotest->fdopen( dup(fileno($fh)) , 'r' );
    #  if ($iotest->error == 0) { ... }
    # didn't actually seem to work under any platform, since there would no
    # no error if the filehandle had been opened writable only. Couldn't be
    # hacked around when dealing with unseekable (piped) filehandles.
    #
    # Just try and do a simple readline, turning io warnings off, instead:
    
    my $fh = $obj->_fh || return '?';
    
    no warnings "io"; # we expect a warning if this is writable only
    my $line = <$fh>;
    if (defined $line) {
        $obj->_pushback($line);
        $obj->{'_mode'} = 'r';
    }
    else {
        $obj->{'_mode'} = 'w';
    }
    
    return $obj->{'_mode'};
}

=head2 file

 Title   : file
 Usage   : $obj->file($newval)
 Function: Get/set the filename, if one has been designated.
 Example :
 Returns : value of file
 Args    : newvalue (optional)

=cut

sub file {
    my ($obj, $value) = @_;
    if ( defined $value) {
    $obj->{'_file'} = $value;
    }
    return $obj->{'_file'};
}

=head2 _print

 Title   : _print
 Usage   : $obj->_print(@lines)
 Function:
 Example :
 Returns : 1 on success, undef on failure

=cut

sub _print {
    my $self = shift;
    my $fh = $self->_fh() || \*STDOUT;
    my $ret = print $fh @_;
    return $ret;
}


=head2 _insert

    Title   : _insert
    Usage   : $obj->_insert($string,1)
    Function: Insert some text in a file at the given line number (1-based).
    Returns : 1 on success
    Args    : string to write in file
              line number to insert the string at

=cut

sub _insert {
    my ($self, $string, $line_num) = @_;
    # Line number check
    if ($line_num < 1) {
        $self->throw("Cannot insert text at line $line_num because the minimum".
            " line number possible is 1");
    }
    # File check
    my $file = $self->file;
    if (not defined $file) {
        $self->throw('Cannot insert a line in a IO object initialized with ".
            "anything else than a file.');
    }
    $file =~ s/^\+?[><]?//; # transform '+>output.ace' into 'output.ace'
    # Everything that needs to be written is written before we read it
    $self->flush;
    # Edit the file in place, line by line (no slurping)
    {
        local @ARGV = ($file);     # input file
        #local $^I = '~';          # backup file extension, e.g. ~, .bak, .ori
        local $^I = '';            # no backup file
        while (<>) {
            if ($. == $line_num) { # right line for new data
                print $string.$_;
            } else {
                print;
            }
        }
    }
    # Line number check (again)
    if ( $. > 0 && $line_num > $. ) {
        $self->throw("Cannot insert text at line $line_num because there are ".
            "only $. lines in file $file");
    }
    # Re-open the file in append mode to be ready to add text at the end of it
    # when the next _print() statement comes
    open my $new_fh, ">>$file" or $self->throw("Cannot append to file $file: $!");
    $self->_fh($new_fh);
    # If file is empty and we're inserting at line 1, simply append text to file
    if ( $. == 0 && $line_num == 1 ) {
        $self->_print($string);
    }
    return 1;
}


=head2 _readline

 Title   : _readline
 Usage   : $obj->_readline(%args)
 Function: Reads a line of input.

           Note that this method implicitely uses the value of $/ that is
           in effect when called.

           Note also that the current implementation does not handle pushed
           back input correctly unless the pushed back input ends with the
           value of $/.

 Example :
 Args    : Accepts a hash of arguments, currently only -raw is recognized
           passing (-raw => 1) prevents \r\n sequences from being changed
           to \n.  The default value of -raw is undef, allowing \r\n to be
           converted to \n.
 Returns : 

=cut

sub _readline {
    my $self = shift;
    my %param =@_;
    my $fh = $self->_fh or return;
    my $line;

    # if the buffer been filled by _pushback then return the buffer
    # contents, rather than read from the filehandle
    if( @{$self->{'_readbuffer'} || [] } ) {
    $line = shift @{$self->{'_readbuffer'}};
    } else {
    $line = <$fh>;
    }
    
    #don't strip line endings if -raw is specified
    # $line =~ s/\r\n/\n/g if( (!$param{-raw}) && (defined $line) );
    # Dave Howorth's fix
    if( !$HAS_EOL && !$param{-raw} && (defined $line) ) {
        $line =~ s/\015\012/\012/g; # Change all CR/LF pairs to LF
        $line =~ tr/\015/\n/ unless $ONMAC; # Change all single CRs to NEWLINE
    }
    return $line;
}

=head2 _pushback

 Title   : _pushback
 Usage   : $obj->_pushback($newvalue)
 Function: puts a line previously read with _readline back into a buffer.
           buffer can hold as many lines as system memory permits.
 Example : $obj->_pushback($newvalue)
 Returns : none
 Args    : newvalue
 Note    : This is only supported for pushing back data ending with the
           current, localized value of $/. Using this method to push modified
           data back onto the buffer stack is not supported; see bug 843.

=cut

# fix for bug 843, this reveals some unsupported behavior
    
#sub _pushback {
#    my ($obj, $value) = @_;    
#    if (index($value, $/) >= 0) {
#        push @{$obj->{'_readbuffer'}}, $value;
#    } else {
#        $obj->throw("Pushing modifed data back not supported: $value");
#    }
#}

sub _pushback {
    my ($obj, $value) = @_;
    return unless $value;
    push @{$obj->{'_readbuffer'}}, $value;
}

=head2 close

 Title   : close
 Usage   : $io->close()
 Function: Closes the file handle associated with this IO instance.
           Will not close the FH if  -noclose is specified
 Returns : none
 Args    : none

=cut

sub close {
   my ($self) = @_;

   # don't close if we explicitly asked not to
   return if $self->noclose;

   if( defined( my $fh = $self->{'_filehandle'} )) {
       $self->flush;
       return if     ref $fh eq 'GLOB'
         && (    \*STDOUT == $fh
              || \*STDERR == $fh
              || \*STDIN  == $fh
                    );

       # don't close IO::Strings
       close $fh unless ref $fh && $fh->isa('IO::String');
   }
   $self->{'_filehandle'} = undef;
   delete $self->{'_readbuffer'};
}


=head2 flush

 Title   : flush
 Usage   : $io->flush()
 Function: Flushes the filehandle
 Returns : none
 Args    : none

=cut

sub flush {
  my ($self) = shift;
  
  if( !defined $self->{'_filehandle'} ) {
    $self->throw("Attempting to call flush but no filehandle active");
  }

  if( ref($self->{'_filehandle'}) =~ /GLOB/ ) {
    my $oldh = select($self->{'_filehandle'});
    $| = 1;
    select($oldh);
  } else {
    $self->{'_filehandle'}->flush();
  }
}

=head2 noclose

 Title   : noclose
 Usage   : $obj->noclose($newval)
 Function: Get/Set the NOCLOSE flag - setting this to true will
           prevent a filehandle from being closed
           when an object is cleaned up or explicitly closed
           This is a bit of hack 
 Returns : value of noclose (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub noclose{
    my $self = shift;

    return $self->{'_noclose'} = shift if @_;
    return $self->{'_noclose'};
}

sub _io_cleanup {
    my ($self) = @_;
    $self->close();
    my $v = $self->verbose;

    # we are planning to cleanup temp files no matter what    
    if( exists($self->{'_rootio_tempfiles'}) &&
    ref($self->{'_rootio_tempfiles'}) =~ /array/i &&
    !$self->save_tempfiles) { 
    if( $v > 0 ) {
        warn( "going to remove files ", 
          join(",",  @{$self->{'_rootio_tempfiles'}}), "\n");
    }
    unlink  (@{$self->{'_rootio_tempfiles'}} );
    }
    # cleanup if we are not using File::Temp
    if( $self->{'_cleanuptempdir'} &&
    exists($self->{'_rootio_tempdirs'}) &&
    ref($self->{'_rootio_tempdirs'}) =~ /array/i &&
    !$self->save_tempfiles) {   
    if( $v > 0 ) {
        warn( "going to remove dirs ", 
          join(",",  @{$self->{'_rootio_tempdirs'}}), "\n");
    }
    $self->rmtree( $self->{'_rootio_tempdirs'});
    }
}

=head2 exists_exe

 Title   : exists_exe
 Usage   : $exists = $obj->exists_exe('clustalw');
           $exists = Bio::Root::IO->exists_exe('clustalw')
           $exists = Bio::Root::IO::exists_exe('clustalw')
 Function: Determines whether the given executable exists either as file
           or within the path environment. The latter requires File::Spec
           to be installed.
           On Win32-based system, .exe is automatically appended to the program
           name unless the program name already ends in .exe.
 Example :
 Returns : 1 if the given program is callable as an executable, and 0 otherwise
 Args    : the name of the executable

=cut

sub exists_exe {
    my ($self, $exe) = @_;
    $self->throw("Must pass a defined value to exists_exe") unless defined $exe;
    $exe = $self if (!(ref($self) || $exe));
    $exe .= '.exe' if(($^O =~ /mswin/i) && ($exe !~ /\.(exe|com|bat|cmd)$/i));
    return $exe if ( -f $exe && -x $exe ); # full path and exists

    # Ewan's comment. I don't think we need this. People should not be
    # asking for a program with a pathseparator starting it
    # $exe =~ s/^$PATHSEP//;

    # Not a full path, or does not exist. Let's see whether it's in the path.
    if($FILESPECLOADED) {
        foreach my $dir (File::Spec->path()) {
            my $f = Bio::Root::IO->catfile($dir, $exe);
            return $f if( -f $f && -x $f );
        }
    }
    return 0;
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
    foreach my $key (keys(%params)) {
    if( $key =~ /^-/  ) {
        my $v = $params{$key};
        delete $params{$key};
        $params{uc(substr($key,1))} = $v;
    } else { 
        # this is to upper case
        my $v = $params{$key};
        delete $params{$key};       
        $params{uc($key)} = $v;
    }
    }
    $params{'DIR'} = $TEMPDIR if(! exists($params{'DIR'}));
    unless (exists $params{'UNLINK'} && 
        defined $params{'UNLINK'} &&
        ! $params{'UNLINK'} ) {
    $params{'UNLINK'} = 1;
    } else { $params{'UNLINK'} = 0 }
        
    if($FILETEMPLOADED) {
    if(exists($params{'TEMPLATE'})) {
        my $template = $params{'TEMPLATE'};
        delete $params{'TEMPLATE'};
        ($tfh, $file) = File::Temp::tempfile($template, %params);
    } else {
        ($tfh, $file) = File::Temp::tempfile(%params);
    }
    } else {
    my $dir = $params{'DIR'};
    $file = $self->catfile($dir,
                   (exists($params{'TEMPLATE'}) ?
                $params{'TEMPLATE'} :
                sprintf( "%s.%s.%s",  
                     $ENV{USER} || 'unknown', $$, 
                     $TEMPCOUNTER++)));

    # sneakiness for getting around long filenames on Win32?
    if( $HAS_WIN32 ) {
        $file = Win32::GetShortPathName($file);
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

    if(  $params{'UNLINK'} ) {
    push @{$self->{'_rootio_tempfiles'}}, $file;
    } 


    return wantarray ? ($tfh,$file) : $tfh;
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
    $self->{'_cleanuptempdir'} = ( defined $params{CLEANUP} && 
                   $params{CLEANUP} == 1);
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

    my $force_writable = ($^O eq 'os2' || $^O eq 'dos' || $^O eq 'MSWin32'
               || $^O eq 'amigaos' || $^O eq 'cygwin');
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
    $rp &= 07777;   # don't forget setuid, setgid, sticky bits
    if ( -d _ ) {
        # notabene: 0777 is for making readable in the first place,
        # it's also intended to change it to writable in case we have
        # to recurse in which case we are better than rm -rf for 
        # subtrees with strange permissions
        chmod(0777, ($Is_VMS ? VMS::Filespec::fileify($root) : $root))
          or $self->warn("Can't make directory $root read+writable: $!")
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
          or $self->warn( "Can't make directory $root writable: $!")
        if $force_writable;
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
          or $self->warn( "Can't make file $root writable: $!")
        if $force_writable;
        warn "unlink $root\n" if $verbose;
        # delete all versions under VMS
        for (;;) {
        unless (unlink $root) {
            $self->warn( "Can't unlink file $root: $!");
            if ($force_writable) {
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

=head2 _flush_on_write

 Title   : _flush_on_write
 Usage   : $obj->_flush_on_write($newval)
 Function: Boolean flag to indicate whether to flush 
           the filehandle on writing when the end of 
           a component is finished (Sequences,Alignments,etc)
 Returns : value of _flush_on_write
 Args    : newvalue (optional)


=cut

sub _flush_on_write {
    my ($self,$value) = @_;
    if( defined $value) {
    $self->{'_flush_on_write'} = $value;
    }
    return $self->{'_flush_on_write'};
}

=head2 save_tempfiles

 Title   : save_tempfiles
 Usage   : $obj->save_tempfiles(1)
 Function: Boolean flag to indicate whether to retain tempfiles/tempdir
 Returns : Boolean value : 1 = save tempfiles/tempdirs, 0 = remove (default)
 Args    : Value evaluating to TRUE or FALSE

=cut

sub save_tempfiles {
    my $self = shift;
    if (@_) {
        my $value = shift;
        $self->{save_tempfiles} = $value ? 1 : 0;
    }
    return $self->{save_tempfiles} || 0;
}

1;
