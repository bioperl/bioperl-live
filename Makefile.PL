## --Makefile.PL for Bioperl core module distribution--##
## ----------------------------------------------------##

# $Id$

#
# Basically this is a standard makefile except 
# that we are experimenting with an 'optional autoload'
# on the Seq.pm file. 
#
# At the moment, (oct 1998) this don't work, so it
# is switched off, hence $do_autoload_finesse =0; [EB]
#
# If you want help with this Makefile get in touch with
# the bioperl core by sending a message to bioperl-bugs@bio.perl.org
# for real bugs or vsns-bcd-perl-guts@lists.uni-bielefeld.de
# to send a mail to all the bioperl core developpers.
#

#
# Talk to Ewan (<birney@sanger.ac.uk>) for more info.
#

# Let the code begin...

require 5.004;
use ExtUtils::MakeMaker;

$do_autoload_finesse = 0;
$NAME     = 'Bio';
$DISTNAME = "bioperl";
$VERSION  = "0.05.1";

if( $do_autoload_finesse == 1 ) {
    fatal("Autoload finessing does not work at the moment!");
}

# don't catch signals - it is very boring!

#BEGIN {
##Catch signals and invoke a cleanup function
#$SIG{__WARN__} = $SIG{INT} = $SIG{QUIT} =\&cleanup;
#}

@BACKUP_FILES = qw(./Seq.pm.bak);

#
# This does an inplace edit on Seq.pm
#

if( $do_autoload_finesse ) {

    print "Bioperl module installation\nPreparing modules for AUTOLOADing...";
    `perl -p -i.bak -e 's/#__END__/__END__/' Seq.pm`;
    print "done.\n\n";
}

# Check for LWP and warn that Webblast requires it.

unless( eval "require HTTP::Request::Common" and
	eval "require LWP::UserAgent") {

print STDERR <<QQ_MISSING_LWP_QQ;

======================================================
WARNING: LWP IS NOT INSTALLED.

   Webblast.pm now requires the LWP (libwww-perl) modules.
   To launch remote BLASTs, you will need to install
   these most-excellent modules which can be readily
   obtained from your local CPAN site:
   http://www.perl.com/CPAN/

   But don't worry, the rest of bioperl will work fine  
======================================================

QQ_MISSING_LWP_QQ
}

# Check for Ace and warn that Bio::DB::Ace

unless( eval "require Ace" ) {

print STDERR <<QQ_MISSING_ACE_QQ;

======================================================
WARNING: Ace IS NOT INSTALLED.

  Ace.pm (AcePerl) is not installed. Bio::DB::Ace 
  requires it for access to acedb databases. You can
  download this sweet module from
  http://stein.cshl.org/AcePerl
  
  But don't worry, the rest of bioperl will work fine  
======================================================

QQ_MISSING_ACE_QQ
}

WriteMakefile(
    NAME	=> $NAME,
    DISTNAME    => $DISTNAME,
    VERSION     => $VERSION,
    'dist'      => {  COMPRESS     => 'gzip -9f',
                      SUFFIX       => '.gz', 
                      DIST_DEFAULT => 'all tardist',
                   },

   'realclean'  => {  FILES        => join(' ',@BACKUP_FILES) }
);


sub cleanup {
my $sig = @_;
print "\n** [$0] Caught a SIG$sig -- cleaning up and shutting down.\n";

## Try to keep module pristine & recover from any inplace edits
if(-e "./Seq.pm.bak")    { `cp ./Seq.pm.bak ./Seq.pm.dist`; `mv ./Seq.pm.bak ./Seq.pm`;  }
exit(0);
}



