#!/usr/local/bin/perl -w

# This program is originally from www.webtechs.com  --Alex Dong Li
#
# MODIFICATIONS BY STEVE A. CHERVITZ MARKED 'SAC:'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# INSTALLATION:
#
#  1) Change the first line to the location of perl on your system if
#     perl is not located at above path.
#
#  2) Set $SOCK_STREAM to the correct value for your system (line 38)
#     OR, let the script configure it automatically.
#
#     SAC: Automatic configuration of $SOCK_STREAM by
#          trail and error: there are only two possibilities;
#          try one, if it fails try the other. If that fails, die.
#          This is handy if you tend to run the script on different
#          systems. There is probably a better way to do this.
#
#          Setting $SOCK_STREAM correctly for your system is 
#          recommended since it will be more efficient.
#
#          Curious note: Why doesn't this script call "use Socket"?
#

###############################################################
# Determine which OS you are using
# and set the correct SOCK_STREAM value
#
# Sun 4.1.3 and NT 4.0
#    $SOCK_STREAM = 1;
# Solaris 2.x
#    $SOCK_STREAM = 2;
 
$SOCK_STREAM = 2;


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%                                                                     %
#%             PROGRAM  :  postclient.pl                               %
#%             CREATOR  :  Mark Gaither                                %
#%       CREATION DATE  :  12 Feb 1994                                 %
#%         DESCRIPTION  :  Command line Web client which posts data    %
#% to a CGI script which handles form data.                            %
#%                                                                     %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%                                                                     %
#% Copyright (c) 1996 - WebTechs.  All rights reserved.                %
#% UNPUBLISHED -- rights reserved under the Copyright Laws of the      %
#% United States.  Use of a copyright notice is precautionary only and %
#% does not imply publication or disclosure.                           %
#%                                                                     %
#%                        WebTechs                                     %
#%              1809 Azalea, Cedar Park, TX 78613                      %
#%                                                                     %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Example:
# postclient.pl -u "http://www.webtechs.com/cgi-bin/test-cgi"
#    -f One=1 Two=frog Three=dirty+deed+done+dirt+cheap
$usage = 'Usage: postclient.pl {-h} [-u URL]
    [-f FieldName=FieldValue FieldName=FieldValue FieldName=FieldValue ...]
';

# This script implements its own version of Getopt::Std (sub MyGetOpts)
my $opt_h = 0;  # Usage info
my $opt_u = 0;  # URL
my @opt_f = 0;  # Name-value pairs

&MyGetOpts('hu:f@');
$ENV{'PATH'} = "/usr/ucb:" . $ENV{'PATH'};

if(defined($opt_h)) { print $usage; exit; }
if(!defined($opt_u)) { print $usage; exit; }
else { $url = $opt_u; }
if(!defined(@opt_f)) { print $usage; exit; }
else {
    foreach $f (@opt_f) {
	push(@fields,$f);
    }
}

$ENV{'PATH'} = "/usr/ucb:" . $ENV{'PATH'};

$object = '';

# decipher the URL
if($url =~ /^([a-z]+)\:\/\/([^:\/]+)(:\d+)?(.+)$/) {
    $protocol = $1;
    $server = $2;
    $port = $3;
    $object = $4;

    if(!defined($port)) { $port = 80; }
    else { $port =~ s/://; }
}

if($protocol ne 'http') { exit; }

# build the GET URL command
# if the leading slash is omitted, try adding one
if($object eq '') { $object = "/"; }

$name = $aliases = $proto = '';

$AF_INET = 2;

$SIG{'INT'} = 'dokill';

$sockaddr = 'S n a4 x8';

chop($hostname = `hostname`);

($name,$aliases,$proto) = getprotobyname('tcp');
($name,$aliases,$port) = getservbyname($port,'tcp') unless $port =~ /^\d+$/;;
($name,$aliases,$type,$len,$thisaddr) = gethostbyname($hostname);

# SAC: Added $0 to the die calls.

($name,$aliases,$type,$len,$thataddr) = gethostbyname($server);
if($name eq '') { die "\n$0: $!\n"; } # handle gethostbyname failure

$this = pack($sockaddr,$AF_INET,0,$thisaddr);
$that = pack($sockaddr,$AF_INET,$port,$thataddr);

# make the socket filehandle
#
# SAC: Using trial and error to configure SOCK_STREAM.

if (!socket(S,$AF_INET,$SOCK_STREAM,$proto)) {
#    print STDERR "\n\npostclient.pl: open socket($SOCK_STREAM) failed. Trying again.\n\n";
    $SOCK_STREAM = ($SOCK_STREAM == 1 ? 2 : 1);
    if(!socket(S,$AF_INET,$SOCK_STREAM,$proto)) {
	die "\n$0: Can't open socket: $!\n";
    }
}

# give the socket an address
if (!bind(S,$this)) { die "\n$0: $!\n"; }

#call up the server
if (!connect(S,$that)) { die "\n$0: $!\n"; }

# Build the string of field name and values pairs
my $str = '';
for($i = 0; $i < $#fields; $i++) {
    $str .= $fields[$i] . "&";
}
$str .= $fields[$#fields]; # this appends the last pair

$str =~ tr/ /+/; # replace blanks with '+' per URL spec
$str =~ tr/&/&amp;/; # treat URL ampersands
$length = length($str);

# set socket to be command buffered
select(S);
$| = 1;

# Build the HTTP POST object
print <<"DATA";
POST $object HTTP/1.0
Content-Type:application/x-www-form-urlencoded
Content-Length:$length

$str
DATA

#######################################################
# If you want to see the return HTML resulting from
# the POST operation, uncomment the following lines.
# The following code reads from the socket you just
# posted to. The code chunk is like a GET operation.
#######################################################

select(STDOUT);
$| = 1;
$_ = <S>;

if(/^HTTP/) {
    while(<S>) {
	# look for first blank line
	if(/^\w/) { next; }
	else {
	    while(<S>) { print STDOUT $_; }
	}
    }
}

close(S);

exit;

sub MyGetOpts {
    local($argumentative) = @_;
    local(@args,$_,$first,$rest);
    local($errs) = 0;
    local($[) = 0;

    @args = split( / */, $argumentative );
    while(@ARGV && ($_ = $ARGV[0]) =~ /^-(.)(.*)/) {
	($first,$rest) = ($1,$2);
	$pos = index($argumentative,$first);
	if($pos >= $[) {
	    if($args[$pos+1] eq ':') {
		shift(@ARGV);
		if($rest eq '') {
		    $rest = shift(@ARGV);
		}
		eval "\$opt_$first = \$rest;";
	    }
	    elsif($args[$pos+1] eq '@') {
		undef @array;
		shift(@ARGV);
		while(@ARGV && ($_ = $ARGV[0]) !~ /^-/) {
		    shift(@ARGV);
		    push(@array,$_);
		}
		eval "\@opt_$first = \@array;";
	    }
	    else {
		eval "\$opt_$first = 1";
		if($rest eq '') {
		    shift(@ARGV);
		}
		else {
		    $ARGV[0] = "-$rest";
		}
	    }
	}
	else {
	    print STDERR "Unknown option: $first\n";
	    ++$errs;
	    if($rest ne '') {
		$ARGV[0] = "-$rest";
	    }
	    else {
		shift(@ARGV);
	    }
	}
    }
    $errs == 0;
}

