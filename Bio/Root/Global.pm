#--------------------------------------------------------------------------------
# PACKAGE : Bio::Root::Global.pm
# PURPOSE : Provides global data, objects, and methods potentially useful to 
#           many different modules and scripts.
# AUTHOR  : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED : 3 Sep 1996
# REVISION: $Id$
#
# INSTALLATION:
#   This module is included with the central Bioperl distribution:
#   http://bio.perl.org/Core/Latest
#   ftp://bio.perl.org/pub/DIST
#   Follow the installation instructions included in the README file.
#
# COMMENTS: Edit the $AUTHORITY string to a desired e-mail address.
#
#           STRICTNESS, VERBOSITY, and variables containing the words WARN and FATAL
#           are considered experimental. The purpose & usage of these is explained 
#           in Bio::Root::Object.pm.
#             
# MODIFIED: 
#    sac --- Fri Jan  8 00:04:28 1999
#      * Added BEGIN block to set $CGI if script is running as a cgi.
#    sac --- Tue Dec 1 1998
#      * Added $STRICTNESS and $VERBOSITY.
#      * Deprecated WARN_ON_FATAL, FATAL_ON_WARN, DONT_WARN and related methods.
#        These will eventually be removed.
#    sac --- Fri 5 Jun 1998: Added @DAYS.
#    sac --- Sun Aug 16 1998: Added $RECORD_ERR and &record_err().
#--------------------------------------------------------------------------------
package	 Bio::Root::Global;

BEGIN {
    use vars qw($CGI);

    # $CGI is a boolean to indicate if the script is running as a CGI.
    # Useful for conditionally producing HTML-formatted messages
    # or suppressing messages appropriate only for interactive sessions.

    $CGI = 1 if $ENV{REMOTE_ADDR} || $ENV{REMOTE_HOST};
}

use Exporter ();
use vars qw($BASE_YEAR @DAYS @MONTHS);

@ISA       = qw( Exporter );
@EXPORT_OK = qw($AUTHORITY $NEWLINE
		$DEBUG $MONITOR $TESTING 
		$DONT_WARN $WARN_ON_FATAL $FATAL_ON_WARN $RECORD_ERR
		$STRICTNESS $VERBOSITY
		$CGI $GLOBAL 
		$BASE_YEAR %ROMAN_NUMS @MONTHS @DAYS 
		&roman2int &debug &monitor &testing &dont_warn &record_err
		&warn_on_fatal &fatal_on_warn &strictness &verbosity
		);

%EXPORT_TAGS = (
		
		std   =>[qw($DEBUG $MONITOR $TESTING $NEWLINE
			    $DONT_WARN $WARN_ON_FATAL $FATAL_ON_WARN $RECORD_ERR
			    $STRICTNESS $VERBOSITY			    
			    &debug &monitor &testing &dont_warn 
			    &warn_on_fatal &fatal_on_warn &record_err
			    &strictness &verbosity
			    &roman2int $AUTHORITY $CGI $GLOBAL)],

		obj   =>[qw($GLOBAL)],

		devel =>[qw($DEBUG $MONITOR $TESTING $DONT_WARN 
			    $WARN_ON_FATAL $FATAL_ON_WARN $RECORD_ERR
			    $STRICTNESS $VERBOSITY $NEWLINE		    
			    &debug &monitor &testing &dont_warn 
			    &strictness &verbosity
			    &warn_on_fatal &fatal_on_warn)], 

		data  =>[qw($BASE_YEAR %ROMAN_NUMS @MONTHS @DAYS)],

		);

# Note: record_err() is not included in the devel tag to allow Bio::Root:Object.pm
#       to define it without a name clash.

######################################
##             Data                 ##
######################################

# Who should receive feedback from users and possibly automatic error messages.
$AUTHORITY     = 'sac@genome.stanford.edu';
 
$DEBUG         = 0;
$MONITOR       = 0;
$TESTING       = 0;
$DONT_WARN     = 0;
$WARN_ON_FATAL = 0; 
$FATAL_ON_WARN = 0; 
$RECORD_ERR    = 0;
$STRICTNESS    = 0;
$VERBOSITY     = 0;

$BASE_YEAR = 1900;
$NEWLINE   = undef;

%ROMAN_NUMS  = ('1'=>'I',    '2'=>'II',    '3'=>'III',    '4'=>'IV',    '5'=>'V',
		'6'=>'VI',   '7'=>'VII',   '8'=>'VIII',   '9'=>'IX',   '10'=>'X',
               '11'=>'XI',  '12'=>'XII',  '13'=>'XIII',  '14'=>'XIV',  '15'=>'XV',
               '16'=>'XVI', '17'=>'XVII', '18'=>'XVIII', '19'=>'XIX',  '20'=>'XX',
               '21'=>'XXI', '22'=>'XXII', 
		);

@MONTHS = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
@DAYS   = qw(Sun Mon Tue Wed Thu Fri Sat);

# The implicit global object. Used for trapping miscellaneous errors/exceptions.
# Created without using or requiring Bio::Root::Object.pm, because Object.pm uses Global.pm.
# Just be sure to use Bio::Root::Object.pm, or a module that uses it.

$GLOBAL = {};
bless $GLOBAL, 'Bio::Root::Object';
$GLOBAL->{'_name'} = 'Global object';


######################################
##         Methods                  ##
######################################

sub roman2int {
    my $roman = uc(shift);
    foreach (keys %ROMAN_NUMS) {
	return $_ if $ROMAN_NUMS{$_} eq $roman;
    }
# Alternatively:
#    my @int = grep $ROMAN_NUMS{$_} eq $roman, keys %ROMAN_NUMS;
#    return $int[0];
    undef;
}

sub debug {
    my $level = shift;
    if( defined $level) { $DEBUG = $level }
    else { $DEBUG = 0 }
#    $MONITOR and do{ print STDERR $DEBUG ? "Debug on ($DEBUG).\n\n" : "Debug off.\n\n"; };
    $MONITOR and do{ print STDERR $DEBUG ? "Debug on ($DEBUG).\n\n" : ""; };
    $DEBUG;
}

sub monitor {
    my $level = shift;
    if( defined $level) { $MONITOR = $level }
    else { $MONITOR = 0 }
    $DEBUG and (print STDERR "Monitor on ($MONITOR).\n\n");
    $MONITOR;
}

sub testing {
    my $level = shift;
    if( defined $level) { $TESTING = $level }
    else { $TESTING = 0 }
    $TESTING ? ($MONITOR && print STDERR "Testing on ($TESTING).\n\n") : ($MONITOR && print STDERR "Testing off.\n\n");
    $TESTING;
}

sub strictness {
# Values can integers from -2 to 2
# See Bio::Root::Object::strict() for more explanation.
    my $arg = shift;
    if( defined $arg) { $STRICTNESS = $arg}
    $DEBUG && print STDERR "\n*** STRICTNESS: $arg ***\n\n";
    $STRICTNESS;
}

sub verbosity {
# Values can integers from -1 to 1
# See Bio::Root::Object::verbose() for more explanation.
    my $arg = shift;
    if( defined $arg) { $VERBOSITY = $arg}
    $DEBUG && print STDERR "\n*** VERBOSITY: $arg ***\n\n";
    $VERBOSITY;
}

sub record_err {
    if( defined shift) { $RECORD_ERR = 1}
    else { $RECORD_ERR = 0 }
    $RECORD_ERR ? ($DEBUG && print STDERR "\n*** RECORD_ERR on. ***\n\n") : ($DEBUG && print STDERR "RECORD_ERR off.\n\n");
    $RECORD_ERR;
}

##
## The following methods are deprecated and will eventually be removed.
##

sub dont_warn {
    my $arg = shift;
    !$CGI and print STDERR "\n$0: Deprecated method dont_warn() called. Use verbosity(-1) instead\n";
    if( $arg) { verbosity(-1)}
    else { verbosity(0); }
}

sub warn_on_fatal {
    my $arg = shift;
    !$CGI and print STDERR "\n$0: Deprecated method warn_on_fatal() called. Use strictness(-2) instead\n";
    if( $arg) { strictness(-2)}
    else { strictness(0); }
}

sub fatal_on_warn {
    my $arg = shift;
    !$CGI and print STDERR "\n$0: Deprecated method fatal_on_warn() called. Use strictness(2) instead\n";
    if( $arg) { strictness(2)}
    else { strictness(0); }
}

#####################################################################################
#                            END OF PACKAGE 
#####################################################################################

1;
