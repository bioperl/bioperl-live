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



# Novemeber 1998.
#
# A note for CVS maintainers *only*
#
# Currently this makefile is kept *in the wrong place* 
# for the distribution. It should be kept in the directory
# above. This means that all the autoloading finesseing
# needs rewriting once we get around to it. 
#
# Talk to Ewan (<birney@sanger.ac.uk>) for more info.
#

# Let the code begin...

require 5.003;
use ExtUtils::MakeMaker;

$do_autoload_finesse = 0;
$NAME     = 'Bio';
$DISTNAME = "bioperl";
$VERSION  = "0.04";

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



