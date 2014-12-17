package Bio::Root::Utilities;
use strict;
use Bio::Root::IO;
use Bio::Root::Exception;
use base qw(Bio::Root::Root Exporter);

# ABSTRACT: general-purpose utility module
# AUTHOR:   Steve Chervitz <sac@bioperl.org>
# OWNER:    1996-2007 Steve Chervitz
# LICENSE:  Perl_5

=head1 SYNOPSIS

=head2 Object Creation

    # Using the supplied singleton object:
    use Bio::Root::Utilities qw(:obj);
    $Util->some_method();

    # Create an object manually:
    use Bio::Root::Utilities;
    my $util = Bio::Root::Utilities->new();
    $util->some_method();

    $date_stamp = $Util->date_format('yyy-mm-dd');

    $clean = $Util->untaint($dirty);

    $compressed = $Util->compress('/home/me/myfile.txt')

    my ($mean, $stdev) = $Util->mean_stdev( @data );

    $Util->authority("me@example.com");
    $Util->mail_authority("Something you should know about...");

    ...and a host of other methods. See below.

=head1 DESCRIPTION

Provides general-purpose utilities of potential interest to any Perl script.

The C<:obj> tag is a convenience that imports a $Util symbol into your
namespace representing a Bio::Root::Utilities object. This saves you
from creating your own Bio::Root::Utilities object via
C<Bio::Root::Utilities-E<gt>new()> or from prefixing all method calls with
C<Bio::Root::Utilities>, though feel free to do these things if desired.
Since there should normally not be a need for a script to have more
than one Bio::Root::Utilities object, this module thus comes with it's
own singleton.

=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://www.bioperl.org/wiki/Getting_BioPerl
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.

=head1 DEPENDENCIES

Inherits from L<Bio::Root::Root>, and uses L<Bio::Root::IO>
and L<Bio::Root::Exception>.

Relies on external executables for file compression/uncompression
and sending mail. No paths to these are hard coded but are located
as needed.

=head1 SEE ALSO

  http://bioperl.org  - Bioperl Project Homepage

=head1 ACKNOWLEDGEMENTS

This module was originally developed under the auspices of the
Saccharomyces Genome Database: http://www.yeastgenome.org/

=cut

use vars qw(@EXPORT_OK %EXPORT_TAGS);
@EXPORT_OK   = qw($Util);
%EXPORT_TAGS = ( obj => [qw($Util)],
                 std => [qw($Util)],);

use vars qw($ID $Util $GNU_PATH $TIMEOUT_SECS
            @COMPRESSION_UTILS @UNCOMPRESSION_UTILS
            $DEFAULT_NEWLINE $NEWLINE $AUTHORITY
            @MONTHS @DAYS $BASE_YEAR $DEFAULT_CENTURY
            );

$ID = 'Bio::Root::Utilities';
# Number of seconds to wait before timing out when reading input (taste_file())
$TIMEOUT_SECS  = 30;
$NEWLINE = $ENV{'NEWLINE'} || undef;
$BASE_YEAR = 1900;  # perl's localtime() assumes this for it's year data.
# TODO: update this every hundred years. Y2K-sensitive code.
$DEFAULT_CENTURY = $BASE_YEAR + 100;
@MONTHS = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
@DAYS   = qw(Sun Mon Tue Wed Thu Fri Sat);
# Sets the preference for compression utilities to be used by compress().
# The first executable in this list to be found in the current PATH will be used,
# unless overridden in the call to that function. See docs for details.
@COMPRESSION_UTILS = qw(gzip bzip2 zip compress);
@UNCOMPRESSION_UTILS = qw(gunzip gzip bunzip2 unzip uncompress);

# Default person to receive feedback from users and possibly automatic error messages.
$AUTHORITY = '';

# Note: $GNU_PATH is now deprecated, shouldn't be needed since now this module
# will automatically locate the compression utility in the current PATH.
# Retaining $GNU_PATH for backward compatibility.
#
# $GNU_PATH points to the directory containing the gzip and gunzip
# executables. It may be required for executing gzip/gunzip
# in some situations (e.g., when $ENV{PATH} doesn't contain this dir.
# Customize $GNU_PATH for your site if the compress() or
# uncompress() functions are generating exceptions.
$GNU_PATH  = '';
#$GNU_PATH  = '/tools/gnu/bin/';

$DEFAULT_NEWLINE = "\012";  # \n  (used if get_newline() fails for some reason)

## Static UTIL object.
$Util = Bio::Root::Root->new();


=head2 date_format

 Title     : date_format
 Usage     : $Util->date_format( [FMT], [DATE])
 Purpose   : -- Get a string containing the formated date or time
           :    taken when this routine is invoked.
           : -- Provides a way to avoid using `date`.
           : -- Provides an interface to localtime().
           : -- Interconverts some date formats.
           :
           : (For additional functionality, use Date::Manip or
           :  Date::DateCalc available from CPAN).
 Example   : $Util->date_format();
           : $date = $Util->date_format('yyyy-mmm-dd', '11/22/92');
 Returns   : String (unless 'list' is provided as argument, see below)
           :
           :   'yyyy-mm-dd'  = 1996-05-03    # default format.
           :   'yyyy-dd-mm'  = 1996-03-05
           :   'yyyy-mmm-dd' = 1996-May-03
           :   'd-m-y'       = 3-May-1996
           :   'd m y'       = 3 May 1996
           :   'dmy'         = 3may96
           :   'mdy'         = May 3, 1996
           :   'ymd'         = 96may3
           :   'md'          = may3
           :   'year'        = 1996
           :   'hms'         = 23:01:59  # when not converting a format, 'hms' can be
           :                             # tacked on to any of the above options
           :                             # to add the time stamp: eg 'dmyhms'
           :   'full' | 'unix' = UNIX-style date: Tue May  5 22:00:00 1998
           :   'list'          = the contents of localtime(time) in an array.
 Argument  : (all are optional)
           : FMT  = yyyy-mm-dd | yyyy-dd-mm | yyyy-mmm-dd |
           :        mdy | ymd | md | d-m-y | hms | hm
           :        ('hms' may be appended to any of these to
           :        add a time stamp)
           :
           : DATE = String containing date to be converted.
           :        Acceptable input formats:
           :           12/1/97 (for 1 December 1997)
           :           1997-12-01
           :           1997-Dec-01
 Throws    :
 Comments  : If you don't care about formatting or using backticks, you can
           : always use: $date = `date`;
           :
           : For more features, use Date::Manip.pm, (which I should
           : probably switch to...)

See Also   : L<file_date()|file_date>, L<month2num()|month2num>

=cut

#---------------'
sub date_format {
#---------------
    my $self   = shift;
    my $option = shift;
    my $date   = shift;  # optional date to be converted.

    my($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst);

    $option ||= 'yyyy-mm-dd';

    my ($month_txt, $day_txt, $month_num, $fullYear);
    my ($converting, @date);

    # Load a supplied date for conversion:
    if(defined($date) && ($date =~ /[\D-]+/)) {
        $converting = 1;
        if( $date =~ m{/}) {
            ($mon,$mday,$year) = split(m{/}, $date);
        } elsif($date =~ /(\d{4})-(\d{1,2})-(\d{1,2})/) {
            ($year,$mon,$mday) = ($1, $2, $3);
        } elsif($date =~ /(\d{4})-(\w{3,})-(\d{1,2})/) {
            ($year,$mon,$mday) = ($1, $2, $3);
            $mon = $self->month2num($2);
        } else {
            print STDERR "\n*** Unsupported input date format: $date\n";
        }
        if(length($year) == 4) {
            $fullYear = $year;
            $year = substr $year, 2;
        } else {
            # Heuristics to guess what century was intended when a 2-digit year is given
            # If number is over 50, assume it's for prev century; under 50 = default century.
            # TODO: keep an eye on this Y2K-sensitive code
            if ($year > 50) {
                $fullYear = $DEFAULT_CENTURY + $year - 100;
            } else {
                $fullYear = $DEFAULT_CENTURY + $year;
            }
        }
        $mon -= 1;
    } else {
        ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = @date =
            localtime(($date ? $date : time()));
        return @date if $option =~ /list/i;
        $fullYear = $BASE_YEAR+$year;
    }
    $month_txt = $MONTHS[$mon];
    $day_txt   = $DAYS[$wday] if defined $wday;
    $month_num = $mon+1;

#    print "sec: $sec, min: $min, hour: $hour, month: $mon, m-day: $mday, year: $year\nwday: $wday, yday: $yday, dst: $isdst";<STDIN>;

    if( $option =~ /yyyy-mm-dd/i ) {
        $date = sprintf "%4d-%02d-%02d",$fullYear,$month_num,$mday;
    } elsif( $option =~ /yyyy-dd-mm/i ) {
        $date = sprintf "%4d-%02d-%02d",$fullYear,$mday,$month_num;
    } elsif( $option =~ /yyyy-mmm-dd/i ) {
        $date = sprintf "%4d-%3s-%02d",$fullYear,$month_txt,$mday;
    } elsif( $option =~ /full|unix/i ) {
        $date = sprintf "%3s %3s %2d %02d:%02d:%02d %d",$day_txt, $month_txt, $mday, $hour, $min, $sec, $fullYear;
    } elsif( $option =~ /mdy/i ) {
        $date = "$month_txt $mday, $fullYear";
    } elsif( $option =~ /ymd/i ) {
        $date = $year."\l$month_txt$mday";
    } elsif( $option =~ /dmy/i ) {
        $date = $mday."\l$month_txt$year";
    } elsif( $option =~ /md/i ) {
        $date = "\l$month_txt$mday";
    } elsif( $option =~ /d-m-y/i ) {
        $date = "$mday-$month_txt-$fullYear";
    } elsif( $option =~ /d m y/i ) {
        $date = "$mday $month_txt $fullYear";
    } elsif( $option =~ /year/i ) {
        $date = $fullYear;
    } elsif( $option =~ /dmy/i ) {
        $date = $mday.'-'.$month_txt.'-'.$fullYear;
    } elsif($option and $option !~ /hms/i) {
        print STDERR "\n*** Unrecognized date format request: $option\n";
    }

    if( $option =~ /hms/i and not $converting) {
        $date .= " $hour:$min:$sec" if $date;
        $date ||= "$hour:$min:$sec";
    }

    return $date || join(" ", @date);
}


=head2 month2num

 Title      : month2num
 Purpose    : Converts a string containing a name of a month to integer
            : representing the number of the month in the year.
 Example    : $Util->month2num("march");  # returns 3
 Argument   : The string argument must contain at least the first
            : three characters of the month's name. Case insensitive.
 Throws     : Exception if the conversion fails.

=cut

#--------------'
sub month2num {
#--------------
    my ($self, $str) = @_;

    # Get string in proper format for conversion.
    $str = substr($str, 0, 3);
    for my $month (0..$#MONTHS) {
        return $month+1 if $str =~ /$MONTHS[$month]/i;
    }
    $self->throw("Invalid month name: $str");
}

=head2 num2month

 Title   : num2month
 Purpose : Does the opposite of month2num.
         : Converts a number into a string containing a name of a month.
 Example : $Util->num2month(3);  # returns 'Mar'
 Throws  : Exception if supplied number is out of range.

=cut

#-------------
sub num2month {
#-------------
    my ($self, $num) = @_;

    $self->throw("Month out of range: $num") if $num < 1 or $num > 12;
    return $MONTHS[$num-1];
}

=head2 compress

 Title     : compress
 Usage     : $Util->compress(full-path-filename);
           : $Util->compress(<named parameters>);
 Purpose   : Compress a file.
 Example   : $Util->compress("/usr/people/me/data.txt");
           : $Util->compress(-file=>"/usr/people/me/data.txt",
           :                 -tmp=>1,
           :                 -outfile=>"/usr/people/share/data.txt.gz",
           :                 -exe=>"/usr/local/bin/fancyzip");
 Returns   : String containing full, absolute path to compressed file
 Argument  : Named parameters (case-insensitive):
           :   -FILE => String (name of file to be compressed, full path).
           :            If the supplied filename ends with '.gz' or '.Z',
           :            that extension will be removed before attempting to compress.
           : Optional:
           :   -TMP  => boolean. If true, (or if user is not the owner of the file)
           :            the file is compressed to a temp file. If false, file may be
           :            clobbered with the compressed version (if using a utility like
           :            gzip, which is the default)
           :   -OUTFILE => String (name of the output compressed file, full path).
           :   -EXE  => Name of executable for compression utility to use.
           :            Will supercede those in @COMPRESSION_UTILS defined by
           :            this module. If the absolute path to the executable is not provided,
           :            it will be searched in the PATH env variable.
 Throws    : Exception if file cannot be compressed.
           : If user is not owner of the file, generates a warning and compresses to
           : a tmp file. To avoid this warning, use the -o file test operator
           : and call this function with -TMP=>1.
 Comments  : Attempts to compress using utilities defined in the @COMPRESSION_UTILS
           : defined by this module, in the order defined. The first utility that is
           : found to be executable will be used. Any utility defined in optional -EXE param
           : will be tested for executability first.
           : To minimize security risks, the -EXE parameter value is untained using
           : the untaint() method of this module (in 'relaxed' mode to permit path separators).

See Also   : L<uncompress()|uncompress>

=cut

#------------'
sub compress {
#------------
    my ($self, @args) = @_;
    # This method formerly didn't use named params and expected fileName, tmp
    # in that order. This should be backward compatibile.
    my ($fileName, $tmp, $outfile, $exe) = $self->_rearrange([qw(FILE TMP OUTFILE EXE)], @args);
    my ($file, $get, $fmt);

    # in case the supplied name already has a compressed extension
    if($fileName =~ /(\.gz|\.Z|\.bz2|\.zip)$/) { $fileName =~ s/$1$//; };
    $self->debug("compressing file $fileName");

    my @util_to_use = @COMPRESSION_UTILS;

    if (defined $exe){
        $exe = $self->untaint($exe, 1);
        unshift @util_to_use, $exe;
    }

    my @checked = @util_to_use;
    $exe ||= '';
    while (not -x $exe and scalar(@util_to_use)) {
        $exe = $self->find_exe(shift @util_to_use);
    }

    unless (-x $exe) {
        $self->throw("Can't find compression utility. Looked for @checked");
    }

    my ($compressed, @cmd, $handle);

    if(defined($outfile) or $tmp or not -o $fileName) {
        if (defined $outfile) {
            $compressed = $outfile;
        } else {
            # obtain a temporary file name (not using the handle)
            # and insert some special text to flag it as a bioperl-based temp file
            my $io = Bio::Root::IO->new();
            ($handle, $compressed) = $io->tempfile();
            $compressed .= '.tmp.bioperl.gz';
        }

        # Use double quotes if executable path have empty spaces
        if ($exe =~ m/ /) {
            $exe = "\"$exe\"";
        }

        if ($exe =~ /gzip|bzip2|compress/) {
            @cmd = ("$exe -f < \"$fileName\" > \"$compressed\"");
        } elsif ($exe eq 'zip') {
            @cmd = ("$exe -r \"$fileName.zip\" \"$fileName\"");
        }
        not $tmp and
            $self->warn("Not owner of file $fileName. Compressing to temp file $compressed.");
        $tmp = 1;
    } else {
        # Need to compute the compressed name based on exe since we're returning it.
        $compressed = $fileName;
        if ($exe =~ /gzip/) {
            $compressed .= '.gz';
        } elsif ($exe =~ /bzip2/) {
            $compressed .= '.bz2';
        } elsif ($exe =~ /zip/) {
            $compressed .= '.zip';
        } elsif ($exe =~ /compress/) {
            $compressed .= '.Z';
        }
        if ($exe =~ /gzip|bzip2|compress/) {
            @cmd = ($exe, '-f', $fileName);
        } elsif ($exe eq 'zip') {
            @cmd = ($exe, '-r', "$compressed", $fileName);
        }
    }

    if(system(@cmd) != 0) {
        $self->throw( -class => 'Bio::Root::SystemException',
                      -text => "Failed to compress file $fileName using $exe: $!");
    }

    return $compressed;
}

=head2 uncompress

 Title     : uncompress
 Usage     : $Util->uncompress(full-path-filename);
           : $Util->uncompress(<named parameters>);
 Purpose   : Uncompress a file.
 Example   : $Util->uncompress("/usr/people/me/data.txt");
           : $Util->uncompress(-file=>"/usr/people/me/data.txt.gz",
           :                   -tmp=>1,
           :                   -outfile=>"/usr/people/share/data.txt",
           :                   -exe=>"/usr/local/bin/fancyzip");
 Returns   : String containing full, absolute path to uncompressed file
 Argument  : Named parameters (case-insensitive):
           :   -FILE => String (name of file to be uncompressed, full path).
           :            If the supplied filename ends with '.gz' or '.Z',
           :            that extension will be removed before attempting to uncompress.
           : Optional:
           :   -TMP  => boolean. If true, (or if user is not the owner of the file)
           :            the file is uncompressed to a temp file. If false, file may be
           :            clobbered with the uncompressed version (if using a utility like
           :            gzip, which is the default)
           :   -OUTFILE => String (name of the output uncompressed file, full path).
           :   -EXE  => Name of executable for uncompression utility to use.
           :            Will supercede those in @UNCOMPRESSION_UTILS defined by
           :            this module. If the absolute path to the executable is not provided,
           :            it will be searched in the PATH env variable.
 Throws    : Exception if file cannot be uncompressed.
           : If user is not owner of the file, generates a warning and uncompresses to
           : a tmp file. To avoid this warning, use the -o file test operator
           : and call this function with -TMP=>1.
 Comments  : Attempts to uncompress using utilities defined in the @UNCOMPRESSION_UTILS
           : defined by this module, in the order defined. The first utility that is
           : found to be executable will be used. Any utility defined in optional -EXE param
           : will be tested for executability first.
           : To minimize security risks, the -EXE parameter value is untained using
           : the untaint() method of this module (in 'relaxed' mode to permit path separators).

See Also   : L<compress()|compress>

=cut

#------------'
sub uncompress {
#------------
    my ($self, @args) = @_;
    # This method formerly didn't use named params and expected fileName, tmp
    # in that order. This should be backward compatibile.
    my ($fileName, $tmp, $outfile, $exe) = $self->_rearrange([qw(FILE TMP OUTFILE EXE)], @args);
    my ($file, $get, $fmt);

    # in case the supplied name lacks a compressed extension
    if(not $fileName =~ /(\.gz|\.Z|\.bz2|\.zip)$/) { $fileName .= $1; };
    $self->debug("uncompressing file $fileName");

    my @util_to_use = @UNCOMPRESSION_UTILS;

    if (defined $exe){
        $exe = $self->untaint($exe, 1);
        unshift @util_to_use, $exe;
    }

    $exe ||= '';
    while (not -x $exe and scalar(@util_to_use)) {
        $exe = $self->find_exe(shift @util_to_use);
    }

    unless (-x $exe) {
        $self->throw("Can't find compression utility. Looked for @util_to_use");
    }

    my ($uncompressed, @cmd, $handle);

    $uncompressed = $fileName;
    $uncompressed =~ s/\.\w+$//;

    if(defined($outfile) or $tmp or not -o $fileName) {
        if (defined $outfile) {
            $uncompressed = $outfile;
        } else {
            # obtain a temporary file name (not using the handle)
            my $io = Bio::Root::IO->new();
            ($handle, $uncompressed) = $io->tempfile();
            # insert some special text to flag it as a bioperl-based temp file
            $uncompressed .= '.tmp.bioperl';
        }

        # Use double quotes if executable path have empty spaces
        if ($exe =~ m/ /) {
            $exe = "\"$exe\"";
        }

        if ($exe =~ /gunzip|bunzip2|uncompress/) {
            @cmd = ("$exe -f < \"$fileName\" > \"$uncompressed\"");
        } elsif ($exe =~ /gzip/) {
            @cmd = ("$exe -df < \"$fileName\" > \"$uncompressed\"");
        } elsif ($exe eq 'unzip') {
            @cmd = ("$exe -p \"$fileName\" > \"$uncompressed\"");
        }
        not $tmp and
            $self->warn("Not owner of file $fileName. Uncompressing to temp file $uncompressed.");
        $tmp = 1;
    } else {
        if ($exe =~ /gunzip|bunzip2|uncompress/) {
            @cmd = ($exe, '-f', $fileName);
        } elsif ($exe =~ /gzip/) {
            @cmd = ($exe, '-df', $fileName);
        } elsif ($exe eq 'zip') {
            @cmd = ($exe, $fileName);
        }
    }

    if(system(@cmd) != 0) {
        $self->throw( -class => 'Bio::Root::SystemException',
                      -text => "Failed to uncompress file $fileName using $exe: $!");
    }

    return $uncompressed;
}


=head2 file_date

 Title    : file_date
 Usage    : $Util->file_date( filename [,date_format])
 Purpose  : Obtains the date of a given file.
          : Provides flexible formatting via date_format().
 Returns  : String = date of the file as: yyyy-mm-dd (e.g., 1997-10-15)
 Argument : filename = string, full path name for file
          : date_format = string, desired format for date (see date_format()).
          :               Default = yyyy-mm-dd
 Thows    : Exception if no file is provided or does not exist.
 Comments : Uses the mtime field as obtained by stat().

=cut

#--------------
sub file_date {
#--------------
    my ($self, $file, $fmt) = @_;

    $self->throw("No such file: $file") if not $file or not -e $file;

    $fmt ||= 'yyyy-mm-dd';

    my @file_data = stat($file);
    return $self->date_format($fmt, $file_data[9]); # mtime field
}


=head2 untaint

 Title   : untaint
 Purpose : To remove nasty shell characters from untrusted data
         : and allow a script to run with the -T switch.
         : Potentially dangerous shell meta characters:  &;`'\"|*?!~<>^()[]{}$\n\r
         : Accept only the first block of contiguous characters:
         :  Default allowed chars = "-\w.', ()"
         :  If $relax is true  = "-\w.', ()\/=%:^<>*"
 Usage   : $Util->untaint($value, $relax)
 Returns : String containing the untained data.
 Argument: $value = string
         : $relax = boolean
 Comments:
     This general untaint() function may not be appropriate for every situation.
     To allow only a more restricted subset of special characters
     (for example, untainting a regular expression), then using a custom
     untainting mechanism would permit more control.

     Note that special trusted vars (like $0) require untainting.

=cut

#------------`
sub untaint {
#------------
    my($self,$value,$relax) = @_;
    $relax ||= 0;
    my $untainted;

    $self->debug("\nUNTAINT: $value\n");

    unless (defined $value and $value ne '') {
        return $value;
    }

    if( $relax ) {
        $value =~ /([-\w.\', ()\/=%:^<>*]+)/;
        $untainted = $1
#    } elsif( $relax == 2 ) {  # Could have several degrees of relax.
#        $value =~ /([-\w.\', ()\/=%:^<>*]+)/;
#        $untainted = $1
    } else {
        $value =~ /([-\w.\', ()]+)/;
        $untainted = $1
    }

    $self->debug("UNTAINTED: $untainted\n");

    $untainted;
}


=head2 mean_stdev

 Title    : mean_stdev
 Usage    : ($mean, $stdev) = $Util->mean_stdev( @data )
 Purpose  : Calculates the mean and standard deviation given a list of numbers.
 Returns  : 2-element list (mean, stdev)
 Argument : list of numbers (ints or floats)
 Thows    : n/a

=cut

#---------------
sub mean_stdev {
#---------------
    my ($self, @data) = @_;
    return (undef, undef) if not @data; # case of empty @data list
    my $mean = 0;
    my $N = 0;
    foreach my $num (@data) {
        $mean += $num;
        $N++
    }
    $mean /= $N;
    my $sum_diff_sqd = 0;
    foreach my $num (@data) {
        $sum_diff_sqd += ($mean - $num) * ($mean - $num);
    }
    # if only one element in @data list, unbiased stdev is undefined
    my $stdev = $N <= 1 ? undef : sqrt( $sum_diff_sqd / ($N-1) );
    return ($mean, $stdev);
}


=head2 count_files

 Title    : count_files
 Purpose  : Counts the number of files/directories within a given directory.
          : Also reports the number of text and binary files in the dir
          : as well as names of these files and directories.
 Usage    : count_files(\%data)
          :   $data{-DIR} is the directory to be analyzed. Default is ./
          :   $data{-PRINT} = 0|1; if 1, prints results to STDOUT, (default=0).
 Argument : Hash reference (empty)
 Returns  : n/a;
          : Modifies the hash ref passed in as the sole argument.
          :  $$href{-TOTAL}            scalar
          :  $$href{-NUM_TEXT_FILES}   scalar
          :  $$href{-NUM_BINARY_FILES} scalar
          :  $$href{-NUM_DIRS}         scalar
          :  $$href{-T_FILE_NAMES}     array ref
          :  $$href{-B_FILE_NAMES}     array ref
          :  $$href{-DIRNAMES}         array ref

=cut

#----------------
sub count_files {
#----------------
    my $self = shift;
    my $href = shift;   # Reference to an empty hash.
    my( $name, @fileLine);
    my $dir = $$href{-DIR} || './'; # THIS IS UNIX SPECIFIC? FIXME/TODO
    my $print = $$href{-PRINT} || 0;

    ### Make sure $dir ends with /
    $dir !~ m{/$} and do{ $dir .=  '/'; $$href{-DIR} = $dir; };

    open ( my $PIPE, "ls -1 $dir |" ) || $self->throw("Can't open input pipe: $!");

    ### Initialize the hash data.
    $$href{-TOTAL} = 0;
    $$href{-NUM_TEXT_FILES} = $$href{-NUM_BINARY_FILES} = $$href{-NUM_DIRS} = 0;
    $$href{-T_FILE_NAMES} = [];
    $$href{-B_FILE_NAMES} = [];
    $$href{-DIR_NAMES} = [];
    while( my $line = <$PIPE> ) {
        chomp();
        $$href{-TOTAL}++;
        if( -T $dir.$line ) {
            $$href{-NUM_TEXT_FILES}++;
            push @{$$href{-T_FILE_NAMES}}, $line; }
        if( -B $dir.$line and not -d $dir.$line) {
            $$href{-NUM_BINARY_FILES}++;
            push @{$$href{-B_FILE_NAMES}}, $line; }
        if( -d $dir.$line ) {
            $$href{-NUM_DIRS}++;
            push @{$$href{-DIR_NAMES}}, $line; }
    }
    close $PIPE;

    if( $print) {
        printf( "\n%4d %s\n", $$href{-TOTAL},            "total files+dirs in $dir");
        printf( "%4d %s\n",   $$href{-NUM_TEXT_FILES},   "text files");
        printf( "%4d %s\n",   $$href{-NUM_BINARY_FILES}, "binary files");
        printf( "%4d %s\n",   $$href{-NUM_DIRS},         "directories");
    }
}


=head2 file_info

 Title   : file_info
 Purpose : Obtains a variety of date for a given file.
         : Provides an interface to Perl's stat().
 Status  : Under development. Not ready. Don't use!

=cut

#--------------
sub file_info {
#--------------
    my ($self, %param) = @_;
    my ($file, $get, $fmt) = $self->_rearrange([qw(FILE GET FMT)], %param);
    $get ||= 'all';
    $fmt ||= 'yyyy-mm-dd';

    my($dev, $ino, $mode, $nlink, $uid, $gid, $rdev, $size,
       $atime, $mtime, $ctime, $blksize, $blocks) = stat $file;

    if($get =~ /date/i) {
        ## I can  get the elapsed time since the file was modified but
        ## it's not so straightforward to get the date in a nice format...
        ## Think about using a standard CPAN module for this, like
        ## Date::Manip or Date::DateCalc.

        my $date = $mtime;
        my $elsec = time - $mtime;
        printf "\nFile age: %.0f sec %.0f hrs %.0f days", $elsec, $elsec/3600, $elsec/(3600*24);<STDIN>;
        my $days = sprintf "%.0f", $elsec/(3600*24);
    } elsif($get eq 'all') {
        return stat $file;
    }
}

=head2 delete

 Title   : delete
 Purpose :

=cut

#------------
sub delete {
#------------
    my $self = shift;
    my $fileName = shift;
    if(not -e $fileName) {
        $self->throw("Could not delete file '$fileName': Does not exist.");
    } elsif(not -o $fileName) {
        $self->throw("Could not delete file '$fileName': Not owner.");
    }
    my $ulval = unlink($fileName) > 0
        or $self->throw("Failed to delete file '$fileName': $!");
}


=head2 create_filehandle

 Usage     : $object->create_filehandle(<named parameters>);
 Purpose   : Create a FileHandle object from a file or STDIN.
           : Mainly used as a helper method by read() and get_newline().
 Example   : $data = $object->create_filehandle(-FILE =>'usr/people/me/data.txt')
 Argument  : Named parameters (case-insensitive):
           :  (all optional)
           :    -CLIENT  => object reference for the object submitting
           :                the request. Default = $Util.
           :    -FILE    => string (full path to file) or a reference
           :                to a FileHandle object or typeglob. This is an
           :                optional parameter (if not defined, STDIN is used).
 Returns   : Reference to a FileHandle object.
 Throws    : Exception if cannot open a supplied file or if supplied with a
           : reference that is not a FileHandle ref.
 Comments  : If given a FileHandle reference, this method simply returns it.
           : This method assumes the user wants to read ascii data. So, if
           : the file is binary, it will be treated as a compressed (gzipped)
           : file and access it using gzip -ce. The problem here is that not
           : all binary files are necessarily compressed. Therefore,
           : this method should probably have a -mode parameter to
           : specify ascii or binary.

See Also :  L<get_newline()|get_newline>

=cut

#---------------------
sub create_filehandle {
#---------------------
    my($self, @param) = @_;
    my($client, $file, $handle) =
        $self->_rearrange([qw( CLIENT FILE HANDLE )], @param);

    if(not ref $client) {  $client = $self; }
    $file ||= $handle;
    if( $client->can('file')) {
        $file = $client->file($file);
    }

    my $FH;
    my ($handle_ref);

    if($handle_ref = ref($file)) {
        if($handle_ref eq 'FileHandle') {
            $FH = $file;
            $client->{'_input_type'} = "FileHandle";
        } elsif($handle_ref eq 'GLOB') {
            $FH = $file;
            $client->{'_input_type'} = "Glob";
        } else {
            $self->throw(-class => 'Bio::Root::IOException',
                         -text  => "Could not read file '$file': Not a FileHandle or GLOB ref.");
        }
        $self->verbose > 0 and printf STDERR "$ID: reading data from FileHandle\n";

    } elsif($file) {
        $client->{'_input_type'} = "FileHandle for $file";

        # Use gzip -cd to access compressed data.
        if( -B $file ) {
            $client->{'_input_type'} .= " (compressed)";
            my $gzip = $self->find_exe('gzip');
            $file = "$gzip -cd $file |"
        }

        require FileHandle;
        $FH = FileHandle->new();
        open ($FH, $file) || $self->throw(-class=>'Bio::Root::FileOpenException',
                                          -text =>"Could not access data file '$file': $!");
        $self->verbose > 0 and printf STDERR "$ID: reading data from file '$file'\n";

    } else {
        # Read from STDIN.
        $FH = \*STDIN;
        $self->verbose > 0 and printf STDERR "$ID: reading data from STDIN\n";
        $client->{'_input_type'} = "STDIN";
    }

    return $FH;
}

=head2 get_newline

 Usage     : $object->get_newline(<named parameters>);
 Purpose   : Determine the character(s) used for newlines in a given file or
           : input stream. Delegates to Bio::Root::Utilities::get_newline()
 Example   : $data = $object->get_newline(-CLIENT => $anObj,
           :                                   -FILE =>'usr/people/me/data.txt')
 Argument  : Same arguemnts as for create_filehandle().
 Returns   : Reference to a FileHandle object.
 Throws    : Propogates any exceptions thrown by Bio::Root::Utilities::get_newline().

See Also : L<taste_file()|taste_file>, L<create_filehandle()|create_filehandle>

=cut

#-----------------
sub get_newline {
#-----------------
    my($self, @param) = @_;

    return $NEWLINE if defined $NEWLINE;

    my($client ) =
        $self->_rearrange([qw( CLIENT )], @param);

    my $FH = $self->create_filehandle(@param);

    if(not ref $client) {  $client = $self;   }

    if($client->{'_input_type'} =~ /STDIN|Glob|compressed/) {
        # Can't taste from STDIN since we can't seek 0 on it.
        # Are other non special Glob refs seek-able?
        # Attempt to guess newline based on platform.
        # Not robust since we could be reading Unix files on a Mac, e.g.
        if(defined $ENV{'MACPERL'}) {
            $NEWLINE = "\015";  # \r
        } else {
            $NEWLINE = "\012";  # \n
        }
    } else {
        $NEWLINE = $self->taste_file($FH);
    }

    close ($FH) unless ($client->{'_input_type'} eq 'STDIN' ||
                        $client->{'_input_type'} eq 'FileHandle' ||
                        $client->{'_input_type'} eq 'Glob' );

    delete $client->{'_input_type'};

    return $NEWLINE || $DEFAULT_NEWLINE;
}


=head2 taste_file

 Usage     : $object->taste_file( <FileHandle> );
           : Mainly a utility method for get_newline().
 Purpose   : Sample a filehandle to determine the character(s) used for a newline.
 Example   : $char = $Util->taste_file($FH)
 Argument  : Reference to a FileHandle object.
 Returns   : String containing an octal represenation of the newline character string.
           :   Unix = "\012"  ("\n")
           :   Win32 = "\012\015" ("\r\n")
           :   Mac = "\015"  ("\r")
 Throws    : Exception if no input is read within $TIMEOUT_SECS seconds.
           : Exception if argument is not FileHandle object reference.
           : Warning if cannot determine neewline char(s).
 Comments  : Based on code submitted by Vicki Brown (vlb@deltagen.com).

See Also : L<get_newline()|get_newline>

=cut

#---------------
sub taste_file {
#---------------
    my ($self, $FH) = @_;
    my $BUFSIZ = 256;   # Number of bytes read from the file handle.
    my ($buffer, $octal, $str, $irs, $i);

    ref($FH) eq 'FileHandle' or $self->throw("Can't taste file: not a FileHandle ref");

    $buffer = '';

    # this is a quick hack to check for availability of alarm(); just copied
    # from Bio/Root/IOManager.pm HL 02/19/01
    my $alarm_available = 1;
    eval {
        alarm(0);
    };
    if($@) {
        # alarm() not available (ActiveState perl for win32 doesn't have it.
        # See jitterbug PR#98)
        $alarm_available = 0;
    }
    $SIG{ALRM} = sub { die "Timed out!"; };
    my $result;
    eval {
        $alarm_available && alarm( $TIMEOUT_SECS );
        $result = read($FH, $buffer, $BUFSIZ); # read the $BUFSIZ characters of file
        $alarm_available && alarm(0);
    };
    if($@ =~ /Timed out!/) {
        $self->throw( "Timed out while waiting for input.",
                      "Timeout period = $TIMEOUT_SECS seconds.\n"
                     ."For longer time before timing out, edit \$TIMEOUT_SECS in Bio::Root::Utilities.pm.");

    } elsif(not $result) {
        my $err = $@;
        $self->throw("read taste failed to read from FileHandle.", $err);

    } elsif($@ =~ /\S/) {
        my $err = $@;
        $self->throw("Unexpected error during read: $err");
    }

    seek($FH, 0, 0) or $self->throw("seek failed to seek 0 on FileHandle.");

    my @chars = split(//, $buffer);
    my $flavor;

    for ($i = 0; $i <$BUFSIZ; $i++) {
        if (($chars[$i] eq "\012")) {
            unless ($chars[$i-1] eq "\015") {
                $flavor='Unix';
                $octal = "\012";
                $str = '\n';
                $irs = "^J";
                last;
            }
        } elsif (($chars[$i] eq "\015") && ($chars[$i+1] eq "\012")) {
            $flavor='DOS';
            $octal = "\015\012";
            $str = '\r\n';
            $irs = "^M^J";
            last;
        } elsif (($chars[$i] eq "\015")) {
            $flavor='Mac';
            $octal = "\015";
            $str = '\r';
            $irs = "^M";
            last;
        }
    }
    if (not $octal) {
        $self->warn("Could not determine newline char. Using '\012'");
        $octal = "\012";
    } else {
        #print STDERR "FLAVOR=$flavor, NEWLINE CHAR = $irs\n";
    }
    return($octal);
}

=head2 file_flavor

 Usage     : $object->file_flavor( <filename> );
 Purpose   : Returns the 'flavor' of a given file (unix, dos, mac)
 Example   : print "$file has flavor: ", $Util->file_flavor($file);
 Argument  : filename = string, full path name for file
 Returns   : String describing flavor of file and handy info about line endings.
           : One of these is returned:
           :   unix (\n or 012 or ^J)
           :   dos (\r\n or 015,012 or ^M^J)
           :   mac (\r or 015 or ^M)
           :   unknown
 Throws    : Exception if argument is not a file
           : Propogates any exceptions thrown by Bio::Root::Utilities::get_newline().

See Also : L<get_newline()|get_newline>,  L<taste_file()|taste_file>

=cut

#---------------
sub file_flavor {
#---------------
    my ($self, $file) = @_;
    my %flavors=("\012"    =>'unix (\n or 012 or ^J)',
                 "\015\012" =>'dos (\r\n or 015,012 or ^M^J)',
                 "\015"     =>'mac (\r or 015 or ^M)'
                );

    -f $file or $self->throw("Could not determine flavor: arg '$file' is either non existant or is not a file.\n");
    my $octal = $self->get_newline($file);
    my $flavor = $flavors{$octal} || "unknown";
    return $flavor;
}

######################################
#####     Mail Functions      ########
######################################

=head2 mail_authority

 Title    : mail_authority
 Usage    : $Util->mail_authority( $message )
 Purpose  : Syntactic sugar to send email to $Bio::Root::Global::AUTHORITY

See Also  : L<send_mail()|send_mail>

=cut

#---------------
sub mail_authority {
#---------------
    my( $self, $message ) = @_;
    my $script = $self->untaint($0,1);

    my $email = $self->{'_auth_email'} || $AUTHORITY;
    if (defined $email) {
        $self->send_mail( -TO=>$AUTHORITY, -SUBJ=>$script, -MSG=>$message);
    } else {
        $self->throw("Can't email authority. No email defined.");
    }
}

=head2 authority

 Title    : authority
 Usage    : $Util->authority('admin@example.com');
 Purpose  : Set/get the email address that should be notified by mail_authority()

See Also  : L<mail_authority()|mail_authority>

=cut

#-------------
sub authority {
#-------------
    my( $self, $email ) = @_;
    $self->{'_auth_email'} = $email if defined $email;
    return $self->{'_auth_email'};
}


=head2 send_mail

 Title    : send_mail
 Usage    : $Util->send_mail( named_parameters )
 Purpose  : Provides an interface to mail or sendmail, if available
 Returns  : n/a
 Argument : Named parameters:  (case-insensitive)
          :  -TO   => e-mail address to send to
          :  -SUBJ => subject for message  (optional)
          :  -MSG  => message to be sent   (optional)
          :  -CC   => cc: e-mail address   (optional)
 Thows    : Exception if TO: address appears bad or is missing.
          : Exception if mail cannot be sent.
 Comments : Based on  TomC's tip at:
          :   http://www.perl.com/CPAN/doc/FMTEYEWTK/safe_shellings
          :
          : Using default 'From:' information.
          :   sendmail options used:
          :      -t: ignore the address given on the command line and
          :          get To:address from the e-mail header.
          :     -oi: prevents send_mail from ending the message if it
          :          finds a period at the start of a line.

See Also  : L<mail_authority()|mail_authority>

=cut


#-------------
sub send_mail {
#-------------
    my( $self, @param) = @_;
    my($recipient,$subj,$message,$cc) = $self->_rearrange([qw(TO SUBJ MSG CC)],@param);

    $self->throw("Invalid or missing e-mail address: $recipient")
        if not $recipient =~ /\S+\@\S+/;

    $subj ||= 'empty subject'; $message ||= '';

    # Best to use mail rather than sendmail. Permissions on sendmail in
    # linux distros have been significantly locked down in recent years,
    # due to the perception that it is insecure.
    my ($exe, $ccinfo);
    if ($exe = $self->find_exe('mail')) {
        if (defined $cc) {
            $ccinfo = "-c $cc";
        }
        $self->debug("send_mail: $exe -s '$subj' $ccinfo $recipient\n");
        open (MAIL, "| $exe -s '$subj' $ccinfo $recipient") ||
            $self->throw("Can't send email: mail cannot fork: $!");
        print MAIL <<QQ_EOFM_QQ;
$message
QQ_EOFM_QQ
        $? and $self->warn("mail didn't exit nicely: $?");
        close(MAIL);
    } elsif ($exe = $self->find_exe('sendmail')) {
        open (SENDMAIL, "| $exe -oi -t") ||
            $self->throw("Can't send email: sendmail cannot fork: $!");
        print SENDMAIL <<QQ_EOFSM_QQ;
To: $recipient
Subject: $subj
Cc: $cc

$message

QQ_EOFSM_QQ
        $? and $self->warn("sendmail didn't exit nicely: $?");

        close(SENDMAIL);
    } else {
        $self->throw("Can't find executable for mail or sendmail.");
    }
}


=head2 find_exe

 Title     : find_exe
 Usage     : $Util->find_exe(name);
 Purpose   : Locate an executable (for use in a system() call, e.g.))
 Example   : $Util->find_exe("gzip");
 Returns   : String containing executable that passes the -x test.
             Returns undef if an executable of the supplied name cannot be found.
 Argument  : Name of executable to be found.
           : Can be a full path. If supplied name is not executable, an executable
           : of that name will be searched in all directories in the currently
           : defined PATH environment variable.
 Throws    : No exceptions, but issues a warning if multiple paths are found
           : for a given name. The first one is used.
 Comments  : TODO: Confirm functionality on all bioperl-supported platforms.
             May get tripped up by variation in path separator character used
             for splitting ENV{PATH}.
See Also   :

=cut

#------------
sub find_exe {
#------------
    my ($self, $name) = @_;
    my @bindirs;
    if ($^O =~ m/mswin/i) {
        @bindirs = split ';', $ENV{'PATH'};
        # Add usual executable extension if missing or -x won't work
        $name.= '.exe' if ($name !~ m/\.exe$/i);
    }
    else {
        @bindirs = split ':', $ENV{'PATH'};
    }
    my $exe = $name;
    unless (-x $exe) {
        undef $exe;
        my @exes;
        foreach my $d (@bindirs) {
            # Note: Windows also understand '/' as folder separator,
            # so there is no need to use a conditional with '\'
            push(@exes, "$d/$name") if -x "$d/$name";
        }
        if (scalar @exes) {
            $exe = $exes[0];
            if (defined $exes[1]) {
                $self->warn("find_exe: Multiple paths to '$name' found. Using $exe.");
            }
        }
    }
    return $exe;
}


######################################
###   Interactive Functions      #####
######################################


=head2 yes_reply

 Title   : yes_reply()
 Usage   : $Util->yes_reply( [query_string]);
 Purpose : To test an STDIN input value for affirmation.
 Example : print +( $Util->yes_reply('Are you ok') ? "great!\n" : "sorry.\n" );
         : $Util->yes_reply('Continue') || die;
 Returns : Boolean, true (1) if input string begins with 'y' or 'Y'
 Argument: query_string = string to be used to prompt user (optional)
         : If not provided, 'Yes or no' will be used.
         : Question mark is automatically appended.

=cut

#-------------
sub yes_reply {
#-------------
    my $self = shift;
    my $query = shift;
    my $reply;
    $query ||= 'Yes or no';
    print "\n$query? (y/n) [n] ";
    chomp( $reply = <STDIN> );
    $reply =~ /^y/i;
}



=head2 request_data

 Title   : request_data()
 Usage   : $Util->request_data( [value_name]);
 Purpose : To request data from a user to be entered via keyboard (STDIN).
 Example : $name = $Util->request_data('Name');
         : # User will see: % Enter Name:
 Returns : String, (data entered from keyboard, sans terminal newline.)
 Argument: value_name = string to be used to prompt user.
         : If not provided, 'data' will be used, (not very helpful).
         : Question mark is automatically appended.

=cut

#----------------
sub request_data {
#----------------
    my $self = shift;
    my $data = shift || 'data';
    print "Enter $data: ";
    # Remove the terminal newline char.
    chomp($data = <STDIN>);
    $data;
}

=head2 quit_reply

 Title   : quit_reply
 Usage   :
 Purpose :

=cut

sub quit_reply {
# Not much used since you can use request_data()
# and test for an empty string.
    my $self = shift;
    my $reply;
    chop( $reply = <STDIN> );
    $reply =~ /^q.*/i;
}


=head2 verify_version

 Purpose : Checks the version of Perl used to invoke the script.
         : Aborts program if version is less than the given argument.
 Usage   : verify_version('5.000')

=cut

#------------------
sub verify_version {
#------------------
    my $self = shift;
    my $reqVersion  = shift;

    $] < $reqVersion and do {
        printf STDERR ( "\a\n%s %0.3f.\n", "** Sorry. This Perl script requires at least version", $reqVersion);
        printf STDERR ( "%s %0.3f %s\n\n", "You are running Perl version", $], "Please update your Perl!\n\n" );
        exit(1);
    }
}

1;

__END__
