package Bio::Root::IO;

use strict;
use Symbol;
use IO::Handle;
use File::Copy;
use Fcntl;
use base qw(Bio::Root::Root);

# ABSTRACT: module providing several methods often needed when dealing with file IO
# AUTHOR:   Hilmar Lapp <hlapp@gmx.net>
# OWNER:    Hilmar Lapp
# LICENSE:  Perl_5

# CONTRIBUTOR: Mark A. Jensen <maj@fortinbras.us>

=head1 SYNOPSIS

    # Use stream I/O in your module
    $self->{'io'} = Bio::Root::IO->new(-file => "myfile");
    $self->{'io'}->_print("some stuff");
    my $line = $self->{'io'}->_readline();
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

=cut

our ($FILESPECLOADED,   $FILETEMPLOADED,
     $FILEPATHLOADED,   $TEMPDIR,
     $PATHSEP,          $ROOTDIR,
     $OPENFLAGS,        $VERBOSE,
     $ONMAC,            $HAS_EOL,       );

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

    # If on Win32, attempt to find Win32 package
    if($^O =~ /mswin/i) {
        eval {
            require Win32;
            $HAS_WIN32 = 1;
        };
    }

    # Try to provide a path separator. Why doesn't File::Spec export this,
    # or did I miss it?
    if ($^O =~ /mswin/i) {
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
        # determine open flags for tempfile creation using Fcntl
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
 Usage   : my $io = Bio::Root::IO->new( -file => 'data.txt' );
 Function: Create new class instance. It automatically calls C<_initialize_io>.
 Args    : Same named parameters as C<_initialize_io>.
 Returns : A Bio::Root::IO object

=cut

sub new {
    my ($caller, @args) = @_;
    my $self = $caller->SUPER::new(@args);
    $self->_initialize_io(@args);
    return $self;
}


=head2 _initialize_io

 Title   : _initialize_io
 Usage   : $io->_initialize_io(@params);
 Function: Initializes filehandle and other properties from the parameters.
 Args    : The following named parameters are currently recognized:
              -file     name of file to read or write to
              -fh       file handle to read or write to (mutually exclusive
                        with -file and -string)
              -input    name of file, or filehandle (GLOB or IO::Handle object)
                        to read of write to
              -string   string to read from (will be converted to filehandle)
              -url      name of URL to open
              -flush    boolean flag to autoflush after each write
              -noclose  boolean flag, when set to true will not close a
                        filehandle (must explicitly call close($io->_fh)
              -retries  number of times to try a web fetch before failure
              -ua_parms when using -url, hashref of key => value parameters
                        to pass to LWP::UserAgent->new(). A useful value might
                        be, for example, {timeout => 60 } (ua defaults to 180s)
 Returns : True

=cut

sub _initialize_io {
    my($self, @args) = @_;

    $self->_register_for_cleanup(\&_io_cleanup);

    my ($input, $noclose, $file, $fh, $string,
        $flush, $url, $retries, $ua_parms) =
        $self->_rearrange([qw(INPUT NOCLOSE FILE FH STRING FLUSH URL RETRIES UA_PARMS)],
                          @args);

    my $mode;

    if ($url) {
        $retries ||= 5;

        require LWP::UserAgent;
        my $ua = LWP::UserAgent->new(%$ua_parms);
        my $http_result;
        my ($handle, $tempfile) = $self->tempfile();
        CORE::close($handle);

        for (my $try = 1 ; $try <= $retries ; $try++) {
            $http_result = $ua->get($url, ':content_file' => $tempfile);
            $self->warn("[$try/$retries] tried to fetch $url, but server ".
                        "threw ". $http_result->code . ".  retrying...")
              if !$http_result->is_success;
            last if $http_result->is_success;
        }
        $self->throw("Failed to fetch $url, server threw ".$http_result->code)
          if !$http_result->is_success;

        $file = $tempfile;
        $mode = '>';
    }

    delete $self->{'_readbuffer'};
    delete $self->{'_filehandle'};
    $self->noclose( $noclose) if defined $noclose;
    # determine whether the input is a file(name) or a stream
    if ($input) {
        if (ref(\$input) eq 'SCALAR') {
            # we assume that a scalar is a filename
            if ($file && ($file ne $input)) {
                $self->throw("Input file given twice: '$file' and '$input' disagree");
            }
            $file = $input;
        } elsif (ref($input) &&
            ((ref($input) eq 'GLOB') || $input->isa('IO::Handle'))) {
            # input is a stream
            $fh = $input;
        } else {
            # let's be strict for now
            $self->throw("Unable to determine type of input $input: ".
                         "not string and not GLOB");
        }
    }

    if (defined($file) && defined($fh)) {
        $self->throw("Providing both a file and a filehandle for reading - ".
                     "only one please!");
    }

    if ($string) {
        if (defined($file) || defined($fh)) {
            $self->throw("File or filehandle provided with -string, ".
                         "please unset if you are using -string as a file");
        }
        open $fh, '<', \$string or $self->throw("Could not read string: $!");
    }

    if (defined($file) && ($file ne '')) {
        $self->file($file);
        ($mode, $file) = $self->cleanfile;
        $mode ||= '<';
        my $action = ($mode =~ m/>/) ? 'write' : 'read';
        $fh = Symbol::gensym();
        open $fh, $mode, $file or $self->throw("Could not $action file '$file': $!");
    }

    if (defined $fh) {
        # check filehandle to ensure it's one of:
        # a GLOB reference, as in: open(my $fh, "myfile");
        # an IO::Handle or IO::String object
        # the UNIVERSAL::can added to fix Bug2863
        unless (   ( ref $fh and ( ref $fh eq 'GLOB' ) )
                or ( ref $fh and ( UNIVERSAL::can( $fh, 'can' ) )
                             and (   $fh->isa('IO::Handle')
                                  or $fh->isa('IO::String') ) )
               ) {
            $self->throw("Object $fh does not appear to be a file handle");
        }
        if ($HAS_EOL) {
            binmode $fh, ':raw:eol(LF-Native)';
        }
        $self->_fh($fh); # if $fh not provided, defaults to STDIN and STDOUT
    }

    $self->_flush_on_write(defined $flush ? $flush : 1);

    return 1;
}


=head2 _fh

 Title   : _fh
 Usage   : $io->_fh($newval);
 Function: Get or set the file handle for the stream encapsulated.
 Args    : Optional filehandle to use
 Returns : Filehandle for the stream

=cut

sub _fh {
    my ($self, $value) = @_;
    if ( defined $value) {
        $self->{'_filehandle'} = $value;
    }
    return $self->{'_filehandle'};
}


=head2 mode

 Title   : mode
 Usage   : $io->mode();
           $io->mode(-force => 1);
 Function: Determine if the object was opened for reading or writing
 Args    : -force: Boolean. Once mode() has been called, the mode is cached for
                   further calls to mode(). Use this argument to override this
                   behavior and re-check the object's mode.
 Returns : Mode of the object:
            'r'  for readable
            'w'  for writable
            'rw' for readable and writable
            '?'  if mode could not be determined (e.g. for a -url)

=cut

sub mode {
    my ($self, %arg) = @_;

    # Method 1: IO::Handle::fdopen
    #    my $iotest = new IO::Handle;
    #    $iotest->fdopen( dup(fileno($fh)) , 'r' );
    #    if ($iotest->error == 0) { ... }
    # It did not actually seem to work under any platform, since there would no
    # error if the filehandle had been opened writable only. It could not be
    # hacked around when dealing with unseekable (piped) filehandles.

    # Method 2: readline, a.k.a. the <> operator
    #    no warnings "io";
    #    my $line = <$fh>;
    #    if (defined $line) {
    #       $self->{'_mode'} = 'r';
    #    ...
    # It did not work well either because <> returns undef, i.e. querying the
    # mode() after having read an entire file returned 'w'.

    if ( $arg{-force} || not exists $self->{'_mode'} ) {
        # Determine stream mode
        my $mode;
        my $fh = $self->_fh;
        if (defined $fh) {
            # Determine read/write status of filehandle
            no warnings 'io';
            if ( defined( read $fh, my $content, 0 ) ) {
                # Successfully read 0 bytes
                $mode = 'r'
            }
            if ( defined( syswrite $fh, '') ) {
                # Successfully wrote 0 bytes
                $mode ||= '';
                $mode  .= 'w';
            }
        } else {
           # Stream does not have a filehandle... cannot determine mode
           $mode = '?';
        }
        # Save mode for future use
        $self->{'_mode'} = $mode;
    }
    return $self->{'_mode'};
}


=head2 file

 Title   : file
 Usage   : $io->file('>'.$file);
           my $file = $io->file;
 Function: Get or set the name of the file to read or write.
 Args    : Optional file name (including its mode, e.g. '<' for reading or '>'
           for writing)
 Returns : A string representing the filename and its mode.

=cut

sub file {
    my ($self, $value) = @_;
    if ( defined $value) {
        $self->{'_file'} = $value;
    }
    return $self->{'_file'};
}


=head2 cleanfile

 Title   : cleanfile
 Usage   : my ($mode, $file) = $io->cleanfile;
 Function: Get the name of the file to read or write, stripped of its mode
           ('>', '<', '+>', '>>', etc).
 Args    : None
 Returns : In array context, an array of the mode and the clean filename.

=cut

sub cleanfile {
    my ($self) = @_;
    return ($self->{'_file'} =~ m/^ (\+?[><]{1,2})?\s*(.*) $/x);
}


=head2 format

 Title   : format
 Usage   : $io->format($newval)
 Function: Get the format of a Bio::Root::IO sequence file or filehandle. Every
           object inheriting Bio::Root::IO is guaranteed to have a format.
 Args    : None
 Returns : Format of the file or filehandle, e.g. fasta, fastq, genbank, embl.

=cut

sub format {
    my ($self) = @_;
    my $format = (split '::', ref($self))[-1];
    return $format;
}


=head2 variant

 Title   : format
 Usage   : $io->format($newval)
 Function: Get the variant of a Bio::Root::IO sequence file or filehandle.
           The format variant depends on the specific format used. Note that
           not all formats have variants. Also, the Bio::Root::IO-implementing
           modules that require access to variants need to define a global hash
           that has the allowed variants as its keys.
 Args    : None
 Returns : Variant of the file or filehandle, e.g. sanger, solexa or illumina for
           the fastq format, or undef for formats that do not have variants.

=cut

sub variant {
    my ($self, $variant) = @_;
    if (defined $variant) {
        $variant = lc $variant;
        my $var_name = '%'.ref($self).'::variant';
        my %ok_variants = eval $var_name; # e.g. %Bio::Assembly::IO::ace::variant
        if (scalar keys %ok_variants == 0) {
            $self->throw("Could not validate variant because global variant ".
                         "$var_name was not set or was empty\n");
        }
        if (not exists $ok_variants{$variant}) {
            $self->throw("$variant is not a valid variant of the " .
                         $self->format . ' format');
        }
        $self->{variant} = $variant;
    }
    return $self->{variant};
}


=head2 _print

 Title   : _print
 Usage   : $io->_print(@lines)
 Function: Print lines of text to the IO stream object.
 Args    : List of strings to print
 Returns : True on success, undef on failure

=cut

sub _print {
    my $self = shift;
    my $fh = $self->_fh() || \*STDOUT;
    my $ret = print $fh @_;
    return $ret;
}


=head2 _insert

 Title   : _insert
 Usage   : $io->_insert($string,1)
 Function: Insert some text in a file at the given line number (1-based).
 Args    : * string to write in file
           * line number to insert the string at
 Returns : True

=cut

sub _insert {
    my ($self, $string, $line_num) = @_;
    # Line number check
    if ($line_num < 1) {
        $self->throw("Could not insert text at line $line_num: the minimum ".
                     "line number possible is 1.");
    }
    # File check
    my ($mode, $file) = $self->cleanfile;
    if (not defined $file) {
        $self->throw('Could not insert a line: IO object was initialized with '.
                     'something else than a file.');
    }
    # Everything that needs to be written is written before we read it
    $self->flush;

    # Edit the file line by line (no slurping)
    $self->close;
    my $temp_file;
    my $number = 0;
    while (-e "$file.$number.temp") {
        $number++;
    }
    $temp_file = "$file.$number.temp";
    copy($file, $temp_file);
    open my $fh1, '<', $temp_file or $self->throw("Could not read temporary file '$temp_file': $!");
    open my $fh2, '>', $file      or $self->throw("Could not write file '$file': $!");
    while (my $line = <$fh1>) {
        if ($. == $line_num) { # right line for new data
            print $fh2 $string . $line;
        }
        else {
            print $fh2 $line;
        }
    }
    CORE::close $fh1;
    CORE::close $fh2;
    unlink $temp_file or $self->throw("Could not delete temporary file '$temp_file': $!");

    # Line number check (again)
    if ( $. > 0 && $line_num > $. ) {
        $self->throw("Could not insert text at line $line_num: there are only ".
                     "$. lines in file '$file'");
    }
    # Re-open the file in append mode to be ready to add text at the end of it
    # when the next _print() statement comes
    open my $new_fh, '>>', $file or $self->throw("Could not append to file '$file': $!");
    $self->_fh($new_fh);
    # If file is empty and we're inserting at line 1, simply append text to file
    if ( $. == 0 && $line_num == 1 ) {
        $self->_print($string);
    }
    return 1;
}


=head2 _readline

 Title   : _readline
 Usage   : local $Bio::Root::IO::HAS_EOL = 1;
           my $io = Bio::Root::IO->new(-file => 'data.txt');
           my $line = $io->_readline();
           $io->close;
 Function: Read a line of input and normalize all end of line characters.

           End of line characters are typically "\n" on Linux platforms, "\r\n"
           on Windows and "\r" on older Mac OS. By default, the _readline()
           method uses the value of $/, Perl's input record separator, to
           detect the end of each line. This means that you will not get the
           expected lines if your input has Mac-formatted end of line characters.
           Also, note that the current implementation does not handle pushed
           back input correctly unless the pushed back input ends with the
           value of $/. For each line parsed, its line ending, e.g. "\r\n" is
           converted to "\n", unless you provide the -raw argument.

           Altogether it is easier to let the PerlIO::eol module automatically
           detect the proper end of line character and normalize it to "\n". Do
           so by setting $Bio::Root::IO::HAS_EOL to 1.

 Args    : -raw : Avoid converting end of line characters to "\n" This option
                  has no effect when using $Bio::Root::IO::HAS_EOL = 1.
 Returns : Line of input, or undef when there is nothing to read anymore

=cut

sub _readline {
    my ($self, %param) = @_;
    my $fh = $self->_fh or return;
    my $line;

    # if the buffer been filled by _pushback then return the buffer
    # contents, rather than read from the filehandle
    if( @{$self->{'_readbuffer'} || [] } ) {
        $line = shift @{$self->{'_readbuffer'}};
    } else {
        $line = <$fh>;
    }

    # Note: In Windows the "-raw" parameter has no effect, because Perl already discards
    # the '\r' from the line when reading in text mode from the filehandle
    # ($line = <$fh>), and put it back automatically when printing
    if( !$HAS_EOL && !$param{-raw} && (defined $line) ) {
        # don't strip line endings if -raw or $HAS_EOL is specified
        $line =~ s/\015\012/\012/g;         # Change all CR/LF pairs to LF
        $line =~ tr/\015/\n/ unless $ONMAC; # Change all single CRs to NEWLINE
    }
    return $line;
}


=head2 _pushback

 Title   : _pushback
 Usage   : $io->_pushback($newvalue)
 Function: Puts a line previously read with _readline back into a buffer.
           buffer can hold as many lines as system memory permits.

           Note that this is only supported for pushing back data ending with
           the current, localized value of $/. Using this method to push
           modified data back onto the buffer stack is not supported; see bug
           843.

 Args    : newvalue
 Returns : True

=cut

# fix for bug 843, this reveals some unsupported behavior

#sub _pushback {
#    my ($self, $value) = @_;
#    if (index($value, $/) >= 0) {
#        push @{$self->{'_readbuffer'}}, $value;
#    } else {
#        $self->throw("Pushing modifed data back not supported: $value");
#    }
#}

sub _pushback {
    my ($self, $value) = @_;
    return unless $value;
    unshift @{$self->{'_readbuffer'}}, $value;
    return 1;
}


=head2 close

 Title   : close
 Usage   : $io->close()
 Function: Closes the file handle associated with this IO instance,
           excepted if -noclose was specified.
 Args    : None
 Returns : True

=cut

sub close {
    my ($self) = @_;

    # do not close if we explicitly asked not to
    return if $self->noclose;

    if( defined( my $fh = $self->{'_filehandle'} )) {
        $self->flush;
        return if ref $fh eq 'GLOB' && (
            \*STDOUT == $fh || \*STDERR == $fh || \*STDIN  == $fh
        );

        # don't close IO::Strings
        CORE::close $fh unless ref $fh && $fh->isa('IO::String');
    }
    $self->{'_filehandle'} = undef;
    delete $self->{'_readbuffer'};
    return 1;
}


=head2 flush

 Title   : flush
 Usage   : $io->flush()
 Function: Flushes the filehandle
 Args    : None
 Returns : True

=cut

sub flush {
    my ($self) = shift;

    if( !defined $self->{'_filehandle'} ) {
        $self->throw("Flush failed: no filehandle was active");
    }

    if( ref($self->{'_filehandle'}) =~ /GLOB/ ) {
        my $oldh = select($self->{'_filehandle'});
        $| = 1;
        select($oldh);
    } else {
        $self->{'_filehandle'}->flush();
    }
    return 1;
}


=head2 noclose

 Title   : noclose
 Usage   : $io->noclose($newval)
 Function: Get or set the NOCLOSE flag - setting this to true will prevent a
           filehandle from being closed when an object is cleaned up or
           explicitly closed.
 Args    : Optional new value (a scalar or undef)
 Returns : Value of noclose (a scalar)

=cut

sub noclose {
    my $self = shift;
    return $self->{'_noclose'} = shift if @_;
    return $self->{'_noclose'};
}


=head2 _io_cleanup

=cut

sub _io_cleanup {
    my ($self) = @_;
    $self->close();
    my $v = $self->verbose;

    # we are planning to cleanup temp files no matter what
    if (    exists($self->{'_rootio_tempfiles'})
        and ref($self->{'_rootio_tempfiles'}) =~ /array/i
        and not $self->save_tempfiles
        ) {
        if( $v > 0 ) {
            warn( "going to remove files ",
                  join(",",  @{$self->{'_rootio_tempfiles'}}),
                  "\n");
        }
        unlink  (@{$self->{'_rootio_tempfiles'}} );
    }
    # cleanup if we are not using File::Temp
    if (    $self->{'_cleanuptempdir'}
        and exists($self->{'_rootio_tempdirs'})
        and ref($self->{'_rootio_tempdirs'}) =~ /array/i
        and not $self->save_tempfiles
        ) {
        if( $v > 0 ) {
            warn( "going to remove dirs ",
                  join(",",  @{$self->{'_rootio_tempdirs'}}),
                  "\n");
        }
        $self->rmtree( $self->{'_rootio_tempdirs'});
    }
}


=head2 exists_exe

 Title   : exists_exe
 Usage   : $exists = $io->exists_exe('clustalw');
           $exists = Bio::Root::IO->exists_exe('clustalw')
           $exists = Bio::Root::IO::exists_exe('clustalw')
 Function: Determines whether the given executable exists either as file
           or within the path environment. The latter requires File::Spec
           to be installed.
           On Win32-based system, .exe is automatically appended to the program
           name unless the program name already ends in .exe.
 Args    : Name of the executable
 Returns : 1 if the given program is callable as an executable, and 0 otherwise

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
        for my $dir (File::Spec->path()) {
            my $f = Bio::Root::IO->catfile($dir, $exe);
            return $f if( -f $f && -x $f );
        }
    }
    return 0;
}


=head2 tempfile

 Title   : tempfile
 Usage   : my ($handle,$tempfile) = $io->tempfile();
 Function: Create a temporary filename and a handle opened for reading and
           writing.
           Caveats: If you do not have File::Temp on your system you should
           avoid specifying TEMPLATE and SUFFIX.
 Args    : Named parameters compatible with File::Temp: DIR (defaults to
           $Bio::Root::IO::TEMPDIR), TEMPLATE, SUFFIX.
 Returns : A 2-element array, consisting of temporary handle and temporary
           file name.

=cut

sub tempfile {
    my ($self, @args) = @_;
    my ($tfh, $file);
    my %params = @args;

    # map between naming with and without dash
    for my $key (keys(%params)) {
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
    } else {
        $params{'UNLINK'} = 0;
    }

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
        $file = $self->catfile(
            $dir,
            (exists($params{'TEMPLATE'}) ?
             $params{'TEMPLATE'} :
             sprintf( "%s.%s.%s", $ENV{USER} || 'unknown', $$, $TEMPCOUNTER++))
        );

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
            $self->throw("Could not write temporary file '$file': $!");
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

 Args    : args - ( key CLEANUP ) indicates whether or not to cleanup
           dir on object destruction, other keys as specified by File::Temp
 Returns : The name of a new temporary directory.

=cut

sub tempdir {
    my ($self, @args) = @_;
    if ($FILETEMPLOADED && File::Temp->can('tempdir')) {
        return File::Temp::tempdir(@args);
    }

    # we have to do this ourselves, not good
    # we are planning to cleanup temp files no matter what
    my %params = @args;
    print "cleanup is " . $params{CLEANUP} . "\n";
    $self->{'_cleanuptempdir'} = ( defined $params{CLEANUP} &&
                                   $params{CLEANUP} == 1);
    my $tdir = $self->catfile( $TEMPDIR,
                               sprintf("dir_%s-%s-%s",
                                       $ENV{USER} || 'unknown',
                                       $$,
                                       $TEMPCOUNTER++));
    mkdir($tdir, 0755);
    push @{$self->{'_rootio_tempdirs'}}, $tdir;
    return $tdir;
}


=head2 catfile

 Title   : catfile
 Usage   : $path = Bio::Root::IO->catfile(@dirs, $filename);
 Function: Constructs a full pathname in a cross-platform safe way.

           If File::Spec exists on your system, this routine will merely
           delegate to it. Otherwise it tries to make a good guess.

           You should use this method whenever you construct a path name
           from directory and filename. Otherwise you risk cross-platform
           compatibility of your code.

           You can call this method both as a class and an instance method.

 Args    : components of the pathname (directories and filename, NOT an
           extension)
 Returns : a string

=cut

sub catfile {
    my ($self, @args) = @_;

    return File::Spec->catfile(@args) if $FILESPECLOADED;
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
 Returns : number of files successfully deleted

=cut

# taken straight from File::Path VERSION = "1.0403"
sub rmtree {
    my ($self, $roots, $verbose, $safe) = @_;
    if ( $FILEPATHLOADED ) {
        return File::Path::rmtree ($roots, $verbose, $safe);
    }

    my $force_writable = ($^O eq 'os2' || $^O eq 'dos' || $^O eq 'MSWin32' ||
                          $^O eq 'amigaos' || $^O eq 'cygwin');
    my $Is_VMS = $^O eq 'VMS';

    my @files;
    my $count = 0;
    $verbose ||= 0;
    $safe    ||= 0;
    if ( defined($roots) && length($roots) ) {
        $roots = [$roots] unless ref $roots;
    } else {
        $self->warn("No root path(s) specified\n");
        return 0;
    }

    my $root;
    for $root (@{$roots}) {
        $root =~ s#/\z##;
        (undef, undef, my $rp) = lstat $root or next;
        $rp &= 07777;   # don't forget setuid, setgid, sticky bits
        if ( -d _ ) {
            # notabene: 0777 is for making readable in the first place,
            # it's also intended to change it to writable in case we have
            # to recurse in which case we are better than rm -rf for
            # subtrees with strange permissions
            chmod(0777, ($Is_VMS ? VMS::Filespec::fileify($root) : $root))
              or $self->warn("Could not make directory '$root' read+writable: $!")
            unless $safe;
            if (opendir DIR, $root){
                @files = readdir DIR;
                closedir DIR;
            } else {
                $self->warn("Could not read directory '$root': $!");
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
                print "skipped '$root'\n" if $verbose;
                next;
            }
            chmod 0777, $root
              or $self->warn("Could not make directory '$root' writable: $!")
              if $force_writable;
            print "rmdir '$root'\n" if $verbose;
            if (rmdir $root) {
                ++$count;
            }
            else {
                $self->warn("Could not remove directory '$root': $!");
                chmod($rp, ($Is_VMS ? VMS::Filespec::fileify($root) : $root))
                  or $self->warn("and can't restore permissions to "
                                 . sprintf("0%o",$rp) . "\n");
            }
        }
        else {
            if (     $safe
                and ($Is_VMS ? !&VMS::Filespec::candelete($root)
                             : !(-l $root || -w $root))
                ) {
                print "skipped '$root'\n" if $verbose;
                next;
            }
            chmod 0666, $root
              or $self->warn( "Could not make file '$root' writable: $!")
              if $force_writable;
            warn "unlink '$root'\n" if $verbose;
            # delete all versions under VMS
            for (;;) {
                unless (unlink $root) {
                    $self->warn("Could not unlink file '$root': $!");
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

    return $count;
}


=head2 _flush_on_write

 Title   : _flush_on_write
 Usage   : $io->_flush_on_write($newval)
 Function: Boolean flag to indicate whether to flush
           the filehandle on writing when the end of
           a component is finished (Sequences, Alignments, etc)
 Args    : Optional new value
 Returns : Value of _flush_on_write

=cut

sub _flush_on_write {
    my ($self, $value) = @_;
    if (defined $value) {
        $self->{'_flush_on_write'} = $value;
    }
    return $self->{'_flush_on_write'};
}


=head2 save_tempfiles

 Title   : save_tempfiles
 Usage   : $io->save_tempfiles(1)
 Function: Boolean flag to indicate whether to retain tempfiles/tempdir
 Args    : Value evaluating to TRUE or FALSE
 Returns : Boolean value : 1 = save tempfiles/tempdirs, 0 = remove (default)

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
