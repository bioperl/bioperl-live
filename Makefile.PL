## --Makefile.PL for Bioperl core module distribution--##
## ----------------------------------------------------##


#
# Basically this is a standard Makefile with a hash of
# needed packages which are tested at the start
#
# If you want help with this Makefile get in touch with
# the bioperl core by sending a message to bioperl-bugs@bio.perl.org
# for real bugs or vsns-bcd-perl-guts@lists.uni-bielefeld.de
# to send a mail to all the bioperl core developpers.
#


%packages = (
	     'HTTP::Request::Common' => '0.0/Part of the LWP package/remote http Blast jobs/Bio::Tools::Blast',
	     'LWP::UserAgent' => '0.0/Part of the LWP package/remote http Blast jobs/Bio::Tools::Blast',
	     'Ace'            => '0.0/Aceperl/access of ACeDB database/Bio::DB::Ace'
	     );

sub check_package {
    my ($name,$str) = @_;
    my($ver,$desc,$expl,@modules);

    ($ver,$desc,$expl,$module) = split( /\// , $str);
    if( !eval "require $name") {
	print "External Module $name, $desc,\n is not installed on this computer.\n  The $module in Bioperl needs it for $expl\n\n";
	return 1;
    }

    return 0;
}

#
# Generate sub testing package system
#
#

# make mini test targets for playing around with
# in testing...

opendir(DIR,"t") || warn "No test directory. Weird!";
@files = readdir(DIR);
shift @files;
shift @files;
foreach my $file (@files) {
    $file =~ /(\w+)\.t/ || next;
    $file = $1;
    $line = "test_$file :: pure_all\n\tPERL_DL_NONLAZY=1 \$(FULLPERL) -I\$(INST_ARCHLIB) -I\$(INST_LIB) -I\$(PERL_ARCHLIB) -I\$(PERL_LIB) -e \'use Test::Harness qw(&runtests \$\$verbose); \$\$verbose=\$(TEST_VERBOSE); runtests \@ARGV;\' t/$file.t\n";
    push(@tline,$line);
    push(@targets,$file);
    $tset  .= "$file \\\n";
}

$tline = join('',@tline);
$tline = "show_tests :\n\t\@echo \'type make test_<subtest> to run\'\n\t\@echo '$tset'\n$tline\n";

print STDERR "Generated sub tests. go make show_tests to see available subtests\n";

sub MY::postamble {
             $tline;
           }


				           
#
# Talk to Ewan (<birney@sanger.ac.uk>) for more info.
#

# Let the code begin...

require 5.004;
use ExtUtils::MakeMaker;

$do_autoload_finesse = 0;
$NAME     = 'Bio';
$DISTNAME = "bioperl";
$VERSION  = "0.06.0";

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

$error = 0;
foreach $name ( keys %packages ) {
    if( &check_package($name,$packages{$name}) == 1 ) {
	$error = 1;
    }
}

if( $error == 1 ) {
    print <<QQ_ERROR_QQ;

Warning:

   There are some external packages, listed above, which bioperl 
   uses. This only effects the functionality which is listed above:
   the rest of bioperl will work fine.

   The installation of these external packages is very simple. You
   can read more about the external packages at

   http://bio.perl.org/Download/Core/External.html

   The ftp site

   ftp://bio.perl.org/pub/External/Core/

   has all the external packages for easy retrieval, or you can use CPAN

   Enjoy the rest of bioperl, which you can use after going 'make install'

QQ_ERROR_QQ
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



