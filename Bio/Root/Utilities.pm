#-----------------------------------------------------------------------------
# PACKAGE : Bio::Root::Utilities.pm
# PURPOSE : Provides general-purpose utilities of potential interest to any Perl script.
# AUTHOR  : Steve Chervitz (sac@bioperl.org)
# CREATED : Feb 1996
# REVISION: $Id$
# STATUS  : Alpha
#
# This module manages file compression and uncompression using gzip or
# the UNIX compress programs (see the compress() and uncompress() methods).
# Also, it can create filehandles from gzipped files. If you want to use a
# different compression utility (such as zip, pkzip, stuffit, etc.) you
# are on your own.
#
# If you manage to incorporate an alternate compression utility into this
# module, please post a note to the bio.perl.org mailing list
# bioperl-l@bioperl.org
#
# TODO    : Configure $GNU_PATH during installation.
#           Improve documentation (POD).
#           Make use of Date::Manip and/or Date::DateCalc as appropriate.
#
# MODIFICATIONS: See bottom of file.
#
# Copyright (c) 1996-2000 Steve Chervitz. All Rights Reserved.
#          This module is free software; you can redistribute it and/or
#          modify it under the same terms as Perl itself.
#
#-----------------------------------------------------------------------------

package	Bio::Root::Utilities;
use strict;

BEGIN {
    use vars qw($Loaded_POSIX $Loaded_IOScalar);
    $Loaded_POSIX = 1;
    unless( eval "require POSIX" ) {
	$Loaded_POSIX = 0;
    }
}

use Bio::Root::Global  qw(:data :std $TIMEOUT_SECS);
use Bio::Root::Object  ();
#use AutoLoader;
#*AUTOLOAD = \&AutoLoader::AUTOLOAD;

use vars qw(@EXPORT_OK %EXPORT_TAGS);
use base qw(Bio::Root::Root Exporter);
@EXPORT_OK   = qw($Util);
%EXPORT_TAGS = ( obj => [qw($Util)],
		 std => [qw($Util)],);

use vars qw($ID $Util $GNU_PATH $DEFAULT_NEWLINE);

$ID        = 'Bio::Root::Utilities';

# $GNU_PATH points to the directory containing the gzip and gunzip
# executables. It may be required for executing gzip/gunzip
# in some situations (e.g., when $ENV{PATH} doesn't contain this dir.
# Customize $GNU_PATH for your site if the compress() or
# uncompress() functions are generating exceptions.
$GNU_PATH  = '';
#$GNU_PATH  = '/tools/gnu/bin/';

$DEFAULT_NEWLINE = "\012";  # \n  (used if get_newline() fails for some reason)

## Static UTIL object.
$Util = {};
bless $Util, $ID;
$Util->{'_name'} = 'Static Utilities object';

## POD Documentation:

=head1 NAME

Bio::Root::Utilities - General-purpose utility module

=head1 SYNOPSIS

=head2 Object Creation

    use Bio::Root::Utilities qw(:obj);

There is no need to create a new Bio::Root::Utilities.pm object when
the C<:obj> tag is used. This tag will import the static $Util
object created by Bio::Root::Utilities.pm into your name space. This
saves you from having to call C<new Bio::Root::Utilities>.

You are free to not use the :obj tag and create the object as you
like, but a Bio::Root::Utilities object is not configurable; any given
script only needs a single copy.

    $date_stamp = $Util->date_format('yyy-mm-dd');

    $clean = $Util->untaint($dirty);

    $Util->mail_authority("Something you should know about...");

    ...and other methods. See below.

=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.

=head1 DESCRIPTION

Provides general-purpose utilities of potential interest to any Perl script.
Scripts and modules are expected to use the static $Util object exported by
this package with the C<:obj> tag.

=head1 DEPENDENCIES

L<Bio::Root::Utilities> inherits from L<Bio::Root::Object>.
It also relies on the GNU gzip program for file compression/uncompression.

=head1 SEE ALSO

  Bio::Root::Object.pm       - Core object
  Bio::Root::Global.pm       - Manages global variables/constants

  http://bio.perl.org/                       - Bioperl Project Homepage

  FileHandle.pm (included in the Perl distribution or CPAN).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other Bioperl modules.
Send your comments and suggestions preferably to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR

Steve Chervitz E<lt>sac@bioperl.orgE<gt>

See L<the FEEDBACK section | FEEDBACK> for where to send bug reports and comments.

=head1 VERSION

Bio::Root::Utilities.pm, 0.042

=head1 ACKNOWLEDGEMENTS

This module was developed under the auspices of the Saccharomyces Genome
Database:
    http://genome-www.stanford.edu/Saccharomyces

=head1 COPYRIGHT

Copyright (c) 1997-98 Steve Chervitz. All Rights Reserved.
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


############################################################################
##                 INSTANCE METHODS                                       ##
############################################################################

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
           :   'hms'         = 23:01:59  # 'hms' can be tacked on to any of the above options
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
 Comments  : Relies on the $BASE_YEAR constant exported by Bio:Root::Global.pm.
           :
           : If you don't care about formatting or using backticks, you can
           : always use: $date = `date`;
           :
           : For more features, use Date::Manip.pm, (which I should
           : probably switch to...)

See Also   : L<file_date>(), L<month2num>()

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
    my (@date);

    # Load a supplied date for conversion:
    if(defined($date) && ($date =~ /[\D-]+/)) {
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
	if(length($year) == 4) { $year = substr $year, 2; }
	$mon -= 1;
    } else {
	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = @date =
	    localtime(($date ? $date : time()));
	return @date if $option =~ /list/i;
    }
    $month_txt = $MONTHS[$mon];
    $day_txt   = $DAYS[$wday] if defined $wday;
    $month_num = $mon+1;
    $fullYear = $BASE_YEAR+$year;

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

    if( $option =~ /hms/i) {
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
    for(0..$#MONTHS) {
	return $_+1 if $str =~ /$MONTHS[$_]/i;
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
 Usage     : $Util->compress(filename, [tmp]);
 Purpose   : Compress a file to conserve disk space.
 Example   : $Util->compress("/usr/people/me/data.txt");
 Returns   : String (name of compressed file, full path).
 Argument  : filename = String (name of file to be compressed, full path).
           :            If the supplied filename ends with '.gz' or '.Z',
           :            that extension will be removed before attempting to compress.
           : tmp = boolean,
           :    If true, (or if user is not the owner of the file)
           :         the file is compressed to a tmp file
           :    If false, file is clobbered with the compressed version.
 Throws    : Exception if file cannot be compressed
           : If user is not owner of the file, generates a warning
           :   and compresses to a tmp file.
           :   To avoid this warning, use the -o file test operator
           :   and call this function with a true second argument.
 Comments  : Attempts to compress using gzip (default compression level).
           : If that fails, will attempt to use compress.
           : In some situations, the full path to the gzip executable
           : may be required. This can be specified with the $GNU_PATH
           : package global variable. When installed, $GNU_PATH is an
           : empty string.

See Also   : L<uncompress>()

=cut

#------------'
sub compress {
#------------
    my $self = shift;
    my $fileName = shift;
    my $tmp = shift || 0;

    if($fileName =~ /(\.gz|\.Z)$/) { $fileName =~ s/$1$//; };
    $DEBUG && print STDERR "gzipping file $fileName";

    my ($compressed, @args);

    if($tmp or not -o $fileName) {
	if($Loaded_POSIX) {
	    $compressed = POSIX::tmpnam;
	} else {
	    $compressed = _get_pseudo_tmpnam();
	}
	$compressed .= ".tmp.bioperl";
	$compressed .= '.gz';
	@args = ($GNU_PATH."gzip -f < $fileName > $compressed");
	not $tmp and
	    $self->warn("Not owner of file $fileName\nCompressing to tmp file $compressed.");
	$tmp = 1;
    } else {
	$compressed = "$fileName.gz";
	@args = ($GNU_PATH.'gzip', '-f', $fileName);
    }

    if(system(@args) != 0) {
	# gzip may not be present. Try compress.
	$compressed = "$fileName.Z";
	if($tmp) {
	    @args = ("/usr/bin/compress -f < $fileName > $compressed");
	} else {
	    @args = ('/usr/bin/compress', '-f', $fileName);
	}
	system(@args) == 0 or
	    $self->throw("Failed to gzip/compress file $fileName: $!",
			 "Confirm current \$GNU_PATH: $GNU_PATH",
			 "Edit \$GNU_PATH in Bio::Root::Utilities.pm if necessary.");
    }

    return $compressed;
}


=head2 uncompress

 Title     : uncompress
 Usage     : $Util->uncompress(filename, [tmp]);
 Purpose   : Uncompress a file.
 Example   : $Util->uncompress("/usr/people/me/data.txt.gz");
 Returns   : String (name of uncompressed file, full path).
 Argument  : filename = String (name of file to be uncompressed, full path).
           :           If the supplied filename does not end with '.gz' or '.Z'
           :           a '.gz' will be appended before attempting to uncompress.
           : tmp = boolean,
           :    If true, (or if user is not the owner of the file)
           :         the file is uncompressed to a tmp file
           :    If false, file is clobbered with the uncompressed version.
 Throws    : Exception if file cannot be uncompressed
           : If user is not owner of the file, generates a warning
           :   and uncompresses to a tmp file.
           :   To avoid this warning, use the -o file test operator
           :   and call this function with a true second argument.
 Comments  : Attempts to uncompress using gunzip.
           : If that fails, will use uncompress.
           : In some situations, the full path to the gzip executable
           : may be required. This can be specified with the $GNU_PATH
           : package global variable. When installed, $GNU_PATH is an
           : empty string.

See Also   : L<compress>()

=cut

#---------------
sub uncompress {
#---------------
    my $self = shift;
    my $fileName = shift;
    my $tmp = shift || 0;

    if(not $fileName =~ /(\.gz|\.Z)$/) { $fileName .= '.gz'; }
    $DEBUG && print STDERR "gunzipping file $fileName";

    my($uncompressed, @args);

    if($tmp or not -o $fileName) {
	if($Loaded_POSIX) {
	    $uncompressed = POSIX::tmpnam;
	} else {
	    $uncompressed = _get_pseudo_tmpnam();
	}
	$uncompressed .= ".tmp.bioperl";
	@args = ($GNU_PATH."gunzip -f < $fileName > $uncompressed");
	not $tmp and $self->verbose > 0 and
	    $self->warn("Not owner of file $fileName\nUncompressing to tmp file $uncompressed.");
	$tmp = 1;
    } else {
	@args = ($GNU_PATH.'gunzip', '-f', $fileName);
	($uncompressed = $fileName) =~ s/(\.gz|\.Z)$//;
    }

#    $ENV{'PATH'} = '/tools/gnu/bin';

    if(system(@args) != 0) {
	# gunzip may not be present. Try uncompress.
	($uncompressed = $fileName) =~ s/(\.gz|\.Z)$//;
	if($tmp) {
	    @args = ("/usr/bin/uncompress -f < $fileName > $uncompressed");
	} else {
	    @args = ('/usr/bin/uncompress', '-f', $fileName);
	}
	system(@args) == 0 or
	    $self->throw("Failed to gunzip/uncompress file $fileName: $!",
			 "Confirm current \$GNU_PATH: $GNU_PATH",
			 "Edit \$GNU_PATH in Bio::Root::Utilities.pm if necessary.");
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

    $DEBUG and print STDERR "\nUNTAINT: $value\n";

    defined $value || return;

    if( $relax ) {
	$value =~ /([-\w.\', ()\/=%:^<>*]+)/;
	$untainted = $1
#    } elsif( $relax == 2 ) {  # Could have several degrees of relax.
#	$value =~ /([-\w.\', ()\/=%:^<>*]+)/;
#	$untainted = $1
    } else {
	$value =~ /([-\w.\', ()]+)/;
	$untainted = $1
    }

    $DEBUG and print STDERR "UNTAINTED: $untainted\n";

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
    return (undef,undef) if not @data; # case of empty @data list
    my $mean = 0;
    my $N = 0;
    foreach (@data) { $mean += $_; $N++ }
    $mean /= $N;
    my $sum_diff_sqd = 0;
    foreach (@data) { $sum_diff_sqd += ($mean - $_) * ($mean - $_); }
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
    while( <$PIPE> ) {
	chomp();
	$$href{-TOTAL}++;
	if( -T $dir.$_ ) {
	    $$href{-NUM_TEXT_FILES}++; push @{$$href{-T_FILE_NAMES}}, $_; }
	if( -B $dir.$_ and not -d $dir.$_) {
	    $$href{-NUM_BINARY_FILES}++; push @{$$href{-B_FILE_NAMES}}, $_; }
	if( -d $dir.$_ ) {
	    $$href{-NUM_DIRS}++; push @{$$href{-DIR_NAMES}}, $_; }
    }
    close $PIPE;

    if( $print) {
	printf( "\n%4d %s\n", $$href{-TOTAL}, "total files+dirs in $dir");
	printf( "%4d %s\n", $$href{-NUM_TEXT_FILES}, "text files");
	printf( "%4d %s\n", $$href{-NUM_BINARY_FILES}, "binary files");
	printf( "%4d %s\n", $$href{-NUM_DIRS}, "directories");
    }
}


#=head2 file_info
#
# Title   : file_info
# Purpose : Obtains a variety of date for a given file.
#	  : Provides an interface to Perl's stat().
# Status  : Under development. Not ready. Don't use!
#
#=cut

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


#------------
sub delete {
#------------
  my $self = shift;
  my $fileName = shift;
  if(not -e $fileName) {
    $self->throw("Can't delete file $fileName: Does not exist.");
  } elsif(not -o $fileName) {
    $self->throw("Can't delete file $fileName: Not owner.");
  }
  my $ulval = unlink($fileName) > 0 or
    $self->throw("Failed to delete file $fileName: $!");
}

=head2 create_filehandle

 Usage     : $object->create_filehandle(<named parameters>);
 Purpose   : Create a FileHandle object from a file or STDIN.
           : Mainly used as a helper method by read() and get_newline().
 Example   : $data = $object->create_filehandle(-FILE =>'usr/people/me/data.txt')
 Argument  : Named parameters (case-insensitive):
           :  (all optional)
           :    -CLIENT  => object reference for the object submitting
           :                the request. This facilitates use by
           :                Bio::Root::IOManager::read(). Default = $Util.
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

See Also :  L<get_newline>(), L<Bio::Root::IOManager::read>(),

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
	$self->throw("Can't read from $file: Not a FileHandle or GLOB ref.");
      }
      $self->verbose > 0 and printf STDERR "$ID: reading data from FileHandle\n";

    } elsif($file) {
      $client->{'_input_type'} = "FileHandle for $file";

      # Use gzip -cd to access compressed data.
      if( -B $file ) {
	$client->{'_input_type'} .= " (compressed)";
	$file = "${GNU_PATH}gzip -cd $file |"
      }

      require FileHandle;
      $FH = new FileHandle;
      open ($FH, $file) || $self->throw("Can't access data file: $file",
					"$!");
      $self->verbose > 0 and printf STDERR "$ID: reading data from file $file\n";

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

See Also : L<taste_file>(), L<create_filehandle>()

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

See Also : L<get_newline>()

=cut

#---------------
sub taste_file {
#---------------
  my ($self, $FH) = @_;
  my $BUFSIZ = 256;   # Number of bytes read from the file handle.
  my ($buffer, $octal, $str, $irs, $i);
  my $wait = $TIMEOUT_SECS;

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
    $alarm_available && alarm( $wait );
    $result = read($FH, $buffer, $BUFSIZ); # read the $BUFSIZ characters of file
    $alarm_available && alarm(0);
  };
  if($@ =~ /Timed out!/) {
    $self->throw("Timed out while waiting for input.",
		 "Timeout period = $wait seconds.\nFor longer time before timing out, edit \$TIMEOUT_SECS in Bio::Root::Global.pm.");

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
#    print STDERR "FLAVOR=$flavor, NEWLINE CHAR = $irs\n";
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

See Also : L<get_newline>(),  L<taste_file>()

=cut

#---------------
sub file_flavor {
#---------------
    my ($self, $file) = @_;
    my %flavors=("\012"    =>'unix (\n or 012 or ^J)',
                 "\015\012" =>'dos (\r\n or 015,012 or ^M^J)',
                 "\015"     =>'mac (\r or 015 or ^M)'
                );

    -f $file or $self->throw("Can't determine flavor: arg '$file' is either non existant or is not a file.\n");
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

See Also  : L<send_mail>()

=cut

sub mail_authority {

    my( $self, $message ) = @_;
    my $script = $self->untaint($0,1);

    send_mail( -TO=>$AUTHORITY, -SUBJ=>$script, -MSG=>$message);

}


=head2 send_mail

 Title    : send_mail
 Usage    : $Util->send_mail( named_parameters )
 Purpose  : Provides an interface to /usr/lib/sendmail
 Returns  : n/a
 Argument : Named parameters:  (case-insensitive)
          :  -TO   => e-mail address to send to
          :  -SUBJ => subject for message  (optional)
          :  -MSG  => message to be sent   (optional)
          :  -CC   => cc: e-mail address   (optional)
 Thows    : Exception if TO: address appears bad or is missing
 Comments : Based on  TomC's tip at:
          :   http://www.perl.com/CPAN/doc/FMTEYEWTK/safe_shellings
          :
          : Using default 'From:' information.
          :   sendmail options used:
          :      -t: ignore the address given on the command line and
          :          get To:address from the e-mail header.
          :     -oi: prevents send_mail from ending the message if it
          :          finds a period at the start of a line.

See Also  : L<mail_authority>()

=cut


#-------------'
sub send_mail {
#-------------
    my( $self, @param) = @_;
    my($recipient,$subj,$message,$cc) = $self->_rearrange([qw(TO SUBJ MSG CC)],@param);

    $self->throw("Invalid or missing e-mail address: $recipient")
	if not $recipient =~ /\S+\@\S+/;

    $cc ||= ''; $subj ||= ''; $message ||= '';

    open (SENDMAIL, "|/usr/lib/sendmail -oi -t") ||
	$self->throw("Can't send mail: sendmail cannot fork: $!");

print SENDMAIL <<QQ_EOF_QQ;
To: $recipient
Subject: $subj
Cc: $cc

$message

QQ_EOF_QQ

    close(SENDMAIL);
    if ($?) { warn "sendmail didn't exit nicely: $?" }
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

# Purpose : Returns a string that can be used as a temporary file name.
#           Based on localtime.
#           This is used if POSIX is not available.

sub _get_pseudo_tmpnam {

    my $date = localtime(time());

    my $tmpnam = 'tmpnam';

    if( $date =~ /([\d:]+)\s+(\d+)\s*$/ ) {
    	$tmpnam = $2. '_' . $1;
    	$tmpnam =~ s/:/_/g;
    }
    return $tmpnam;
}


1;
__END__

MODIFICATION NOTES:
---------------------

17 Feb 1999, sac:
  * Using global $TIMEOUT_SECS in taste_file().

13 Feb 1999, sac:
  * Renamed get_newline_char() to get_newline() since it could be >1 char.

3 Feb 1999, sac:
  * Added three new methods: create_filehandle, get_newline_char, taste_file.
    create_filehandle represents functionality that was formerly buried
    within Bio::Root::IOManager::read().

2 Dec 1998, sac:
  * Removed autoloading code.
  * Modified compress(), uncompress(), and delete() to properly
    deal with file ownership issues.

3 Jun 1998, sac:
    * Improved file_date() to be less reliant on the output of ls.
      (Note the word 'less'; it still relies on ls).

5 Jul 1998, sac:
    * compress() & uncompress() will write files to a temporary location
      if the first attempt to compress/uncompress fails.
      This allows users to access compressed files in directories in which they
      lack write permission.



