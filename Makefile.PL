## --Makefile.PL for Bioperl core module distribution--##
## ----------------------------------------------------##
## $Id$
##----------------------------------------------------

#
# Basically this is a standard Makefile with a hash of
# needed packages which are tested at the start
#
# If you want help with this Makefile get in touch with
# the developers by sending a message to bioperl-l@bioperl.org.



# All Bioperl dependencies on external modules are defined here
%packages = (
	     'HTTP::Request::Common'    => '0.0/Part of the LWP package/GenBank+GenPept sequence retrieval, remote http Blast jobs/Bio::DB::*,Bio::Tools::Run::RemoteBlast',
	     'LWP::UserAgent'           => '0.0/Part of the LWP package/GenBank+GenPept sequence retrieval, remote http Blast jobs/Bio::DB::*,Bio::Tools::Run::RemoteBlast',
	     'Ace'                      => '0.0/Aceperl/access of ACeDB database/Bio::DB::Ace',
	     'IO::Scalar'               => '0.0/IO handle to read or write to a scalar/remote http Blast jobs/Bio::Tools::Blast::Run::Webblast',
	     'IO::String'               => '0.0/IO handle to read or write to a string/GenBank+GenPept sequence retrieval, Variation code/Bio::DB::*,Bio::Variation::*,Bio::Tools::Blast::Run::Webblast, Bio::Index::Blast',
	     'XML::Parser'              => '0.0/Parsing of XML documents/Bio::Variation code, GAME parser/Bio::SeqIO::game,Bio::Variation::*',
	     'XML::Writer'              => '0.0/Parsing + writing of XML documents/Bio::Variation code, GAME parser/Bio::SeqIO::game,Bio::Variation::*',
	     'XML::Parser::PerlSAX'              => '0.0/Parsing of XML documents/Bio::Variation code, GAME parser/Bio::SeqIO::game,Bio::Variation::*',
	     'XML::Twig'         	=> '0.0/Available on CPAN/parsing of XML documents/Module Bio::Variation::IO::xml.pm',
	     'File::Temp'               => '0.0/Temporary File creation/Bio::DB::WebDBSeqI, Bio::Seq::LargePrimarySeq',
	     'SOAP::Lite'               => '0.0/SOAP protocol/XEMBL Services (also Bibliographic queries in Biblio::) /Bio::DB::XEMBLService',
	     'HTML::Parser'             => '0.0/HTML parsing of GDB page/Bio::DB::GDB',
	     'DBD::mysql'               => '0.0/Mysql driver/loading and querying of Mysql-based GFF feature databases/Bio::DB::GFF',
	     'GD'                       => '1.3/Graphical Drawing Toolkit/Rendering Sequences and Features/Bio::Graphics',
# rare function let's not scare people
# 	     'srsperl.pm'               => '0.0/Sequence Retrieval System (SRS) perl API/alternative way of retrieving sequences/Bio::LiveSeq::IO::SRS.pm'
	     );

sub check_package {
    my ($name,$str) = @_;
    my ($ver,$desc,$expl,$module) = split( /\// , $str);
    if( !eval "require($name)") {
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
# Prompt user for BioGFFDB stuff
#

if(  eval { require DBI } ) {
    my $proceed = user_prompt("Do you want to run the BioGFFDB tests (requires mysql) y/n",'n');
    if( $proceed eq 'y' ) {
	system("touch t/do_biodbgff.tests");
    } else { 
	if( -e "t/do_biodbgff.tests" ) {
	    unlink "t/do_biodbgff.tests";
	}
    }
}

sub user_prompt ($;$) {
    my($mess,$def)=@_;
    Carp::confess("prompt function called without an argument") unless defined $mess;
    my $dispdef = defined $def ? "[$def] " : " ";
    $def = defined $def ? $def : "";
    my $ans;
    local $|=1;
    print STDERR "$mess $dispdef";
    chomp($ans = <STDIN>);
    return ($ans ne '') ? $ans : $def;
}


				           
#
# Talk to Ewan (<birney@ebi.ac.uk>) or bioperl-l@bioperl.org for more info.
#

# Let the code begin...

require 5.005;
use ExtUtils::MakeMaker;

$do_autoload_finesse = 0;
$NAME     = 'Bio';
$DISTNAME = "bioperl";
$VERSION  = "1.00";

if( $do_autoload_finesse == 1 ) {
    fatal("Autoload finessing does not work at the moment!");
}



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
print "\n";
foreach $name ( keys %packages ) {
    if( &check_package($name,$packages{$name}) == 1 ) {
	$error = 1;
    }
}

if( $error == 1 ) {
    print <<QQ_ERROR_QQ;

Warning:

   There are some external packages and perl modules, listed above, which 
   bioperl uses. This only effects the functionality which is listed above:
   the rest of bioperl will work fine.

   The installation of these external packages is very simple. You
   can read more about bioperl external dependencies at

   http://bioperl.org/Core/external.shtml

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
    ($] ge '5.005')  ? (
	'AUTHOR' => 'Bioperl Team (bioperl-l@bioperl.org)',
	'ABSTRACT' => 'Bioinformatics Toolkit',
	) : (),
     'realclean'  => {  FILES        => join(' ',@BACKUP_FILES) }
);




