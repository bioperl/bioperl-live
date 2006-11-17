## --Makefile.PL for Bioperl core module distribution--##
## ----------------------------------------------------##
## $Id$
##----------------------------------------------------


# This still works, but has been superceded by Build.PL
# Please use Build.PL in preference

warn <<'WARN';
!--
  Bioperl is transitioning from Makefile.PL-based installations to Build.PL

  This Makefile.PL should still work fine, but please try Build.PL instead,
  especially if you run into any non-test-related problems:

  perl Build.PL
  ./Build test
  ./Build install
--!

WARN

#
# Start with object preamble and setup of variables that might need changing
# in the future
#

require 5.006_001;
use strict;
use warnings;
use Config;
use File::Spec::Functions;
use ExtUtils::MakeMaker;
use File::Path 'rmtree';
use Cwd;
use IO::File;

# All Bioperl dependencies on external modules are defined here
# (these are all considered optional)
my %packages = (
  'Ace'                       => '0.0/Aceperl/access of ACeDB database/Bio::DB::Ace,Bio::DB::GFF::Adaptor::ace',
  'Bio::ASN1::EntrezGene'     => '0.0/Entrez Gene Bioperl Parser/parsing entrezgene/Bio::SeqIO::entrezgene',
  # we actually need 1.01 of Class::AutoClass, but unfortunately it is versioned as 1.0
  'Class::AutoClass'          => '1.0/Protocol for object creation/creating objects/Bio::Graph::SimpleGraph,Bio::Graph::SimpleGraph::Traversal,Bio::Graph::ProteinGraph',
  'Clone'                     => '0.0/Recursively copy Perl datatypes/cloning objects/Bio::Graph::ProteinGraph,Bio::Tools::Primer3',
  'Convert::Binary::C'        => '0.0/Binary Data Conversion using C Types/strider functionality/Bio::SeqIO::strider',
  'Data::Stag::XMLWriter'     => '0.0/XML writer/writing choas xml files/Bio::SeqIO::chaosxml',
  'DB_File'                   => '0.0/Persistent object storage and retrieval/persistent database storage/Bio::DB::GFF::Adaptor::berkeleydb,Bio::SeqFeature::Collection,Bio::DB::Taxonomy::flatfile,Bio::DB::FileCache,Bio::DB::Flat::BDB,Bio::Index::Abstract',
  'GD'                        => '1.3/Graphical Drawing Toolkit/rendering Sequences and Features/Bio::Graphics::Glyph*',
  'GD::SVG'                   => '0.0/Generate SVG output/creating SVG images/Bio::Graphics::Panel',
  'Graph::Directed'           => '0.0/Generic Graph data structure and algorithms/ontology engine implementation for the GO parser/Bio::Ontology::SimpleGOEngine::GraphAdaptor',
  'HTML::Entities'            => '0.0/Parse POST request URL forms/remote analysis POST submissions/Bio::SearchIO::blastxml',
  'HTML::Parser'              => '3.0/Parse HTML documents/screen scraping www.gdb.org/Bio::DB::GDB',
  'HTTP::Request::Common'     => '0.0/Part of the LWP package/GenBank+GenPept sequence retrieval, remote http Blast jobs/Bio::DB::*,Bio::Tools::Run::RemoteBlast,Bio::Tools::Analysis::Protein*,Bio::Tools::Analysis::DNA*',
  'IO::String'                => '0.0/IO handle to read or write to a string/a wide variety of functions/Bio::SeqIO::game::gameWriter,Bio::DB::*,Bio::Variation::IO::xml,Bio::Tools::*,Bio::Index::Blast,Bio::Index::Hmmer,Bio::TreeIO::cluster,Bio::TreeIO::nexus',
  'LWP::UserAgent'            => '0.0/Part of the LWP package/remote access/Bio::DB::*,Bio::Tools::Run::RemoteBlast,Bio::Tools::WebBlat,Bio::WebAgent,Bio::Graphics::Glyph::image',
  'Set::Scalar'               => '0.0/Basic set operations/proper operation/Bio::Tree::Compatible',
  'SOAP::Lite'                => '0.0/SOAP protocol/XEMBL Services and Bibliographic queries/Bio::DB::XEMBLService,Bio::DB::Biblio::soap',
  'Spreadsheet::ParseExcel'   => '0.0/Excel parser/parsing Excel files/Bio::SeqIO::excel',
  'Storable'                  => '0.0/Persistent object storage and retrieval/storing sequence objects in local file cache/Bio::DB::FileCache,Bio::SeqFeature::Collection,Bio::PopGen::HtSNP,Bio::PopGen::TagHaplotype,Bio::DB::GFF::Adaptor::berkeleydb',
  'SVG'                       => '2.26/SVG output/Bio::Graphics/Bio::Graphics::Pictogram',
  'SVG::Graph'                => '0.01/SVG output/creating SVG images/Bio::TreeIO::svggraph',
  'Text::Shellwords'          => '0.0/Execute shell commands/test scripts/Bio::Graphics',
  'XML::DOM::XPath'           => '0.0/XML parser/parsing interpro features/Bio::FeatureIO::interpro',
  'XML::Parser'               => '0.0/Parsing of XML documents/parsing xml/Bio::Biblio::IO::medlinexml',
  'XML::Parser::PerlSAX'      => '0.0/Parsing of XML documents/parsing xml/Bio::SeqIO::tinyseq,Bio::SeqIO::game::gameSubs,Bio::OntologyIO::InterProParser,Bio::ClusterIO::dbsnp',
  'XML::SAX'                  => '0.0/New style SAX parser/parsing xml/Bio::SearchIO::blastxml,Bio::SeqIO::tigrxml,Bio::SeqIO::bsml_sax',
  'XML::SAX::Base'            => '0.0/New style SAX parser/parsing xml/Bio::SeqIO::tigrxml,Bio::SeqIO::bsml_sax',
  'XML::SAX::Writer'          => '0.0/New style SAX parser/writing xml/Bio::SeqIO::tigrxml',
  'XML::Simple'               => '0.0/Simple XML parser/parsing xml/Bio::DB::EUtilities',
  'XML::Twig'                 => '0.0/Parsing of XML documents/parsing xml/Bio::Variation::IO::xml,Bio::DB::Taxonomy::entrez,Bio::DB::Biblio::eutils,Bio::Graph::IO::psi_xml',
  'XML::Writer'               => '0.4/Parsing and writing of XML documents/parsing and writing xml/Bio::SeqIO::agave,Bio::SeqIO::game::gameWriter,Bio::SeqIO::chadoxml,Bio::SeqIO::tinyseq,Bio::Variation::IO::xml,Bio::SearchIO::Writer::BSMLResultWriter',
);

# We don't want CPAN, Bundle::BioPerl et al. to know about these for various
# reasons (eg. because they would cause a circular dependency), so will exclude
# them from PREREQ_PM later on
my %exclude_from_prereqs = (
  'Bio::ASN1::EntrezGene' => 'it would cause a circular dependency'
);


#
# Main code
#

# generate string for makefile that allows show_tests and make test_<testname>
# to work
opendir(DIR,"t") || die "No test directory. Cannot continue.";
my @files = readdir(DIR);
shift @files;
shift @files;
my ($tline, $tset);
foreach my $file (@files) {
  $file =~ /(\w+)\.t/ || next;
  $file = $1;
  my $path = catfile('t', $file);
  my $line = "test_$file :: pure_all\n\tPERL_DL_NONLAZY=1 \$(FULLPERL) -I\$(INST_ARCHLIB) -I\$(INST_LIB) -I\$(PERL_ARCHLIB) -I\$(PERL_LIB) -e \'use Test::Harness qw(&runtests \$\$verbose); \$\$verbose=\$(TEST_VERBOSE); runtests \@ARGV;\' $path.t\n";
  $tline .= $line;
  $tset  .= " $file \\\n";
}
$tline = "show_tests :\n\t\@echo \'type make test_<subtest> to run\'\n\t\@echo '$tset'\n$tline\n";
print STDERR "Generated sub tests. go make show_tests to see available subtests\n";

# we may end up creating extra files we would want to clean later; note them here
my @should_be_cleaned;

# ask questions about optional stuff to install/test
my @scripts_to_install = eval {get_scripts_to_install()};
prompt_for_biodbgff();
setup_bioseqfeaturestore();

# notify about optional modules that haven't been installed, and generate hash
# for prereqs needed later in WriteMakefile
my %prereq;
my $error = 0;
print "\n*** Optional dependencies section ***\n";
while (my ($name, $value) = each %packages) {
  my ($ver, $desc, $expl, $module) = split( /\// , $value);
  no strict 'refs';
  # some modules do not define VERSION
  if (! eval "require $name;" ||  (${$name.'::VERSION'} || 0) < $ver) {
    my @modules = split(',', $module);
    my $modules = join("\n * ", @modules);
    print "External Module '$name' >= $ver ($desc) is not installed on this computer.\n The following Bioperl modules need it for $expl:\n * $modules\n";
    $error = 1;
    
    if (exists $exclude_from_prereqs{$name}) {
      print "Additionally, it will not be considered a pre-requisite\n(and therefore will not be automatically installed by CPAN) because:\n $exclude_from_prereqs{$name}.\n\n";
    }
    else {
      print "\n";
    }
  }
  
  unless (exists $exclude_from_prereqs{$name}) {
    $prereq{$name} = $ver;
  }
}
if ($error == 1) {
  print <<QQ_ERROR_QQ;
Information:

  Some Perl modules are missing, listed above and below. The only
  functionality that will be affected is described above;
  the rest of Bioperl will work fine without them.

  The installation of these external packages is very simple. You
  can read more about Bioperl external dependencies in the INSTALL
  file or at:

  http://bioperl.org/wiki/Installing_BioPerl

  Enjoy the rest of Bioperl, which you can use after going 'make install'

QQ_ERROR_QQ
}
else {
	print "All dependencies are installed\n\n";
}

# write the makefile
WriteMakefile(
  NAME         => 'Bio',
  DISTNAME     => 'bioperl',
  VERSION_FROM => "Bio/Root/Version.pm",
  dist         => { COMPRESS     => 'gzip -9f',
                    SUFFIX       => '.gz',
                    DIST_DEFAULT => 'all tardist'
                  },
  'AUTHOR'     => 'Bioperl Team (bioperl-l@bioperl.org)',
  'ABSTRACT'   => 'Bioinformatics Toolkit',
  # clean() is overridden and also gets rid of @should_be_cleaned
  realclean    => { FILES => join(' ', @should_be_cleaned) },
  EXE_FILES    => \@scripts_to_install,
  PL_FILES     => { 'maintenance/symlink_scripts.PLS' => 'maintenance/symlink_scripts.pl' },
  PREREQ_PM    => \%prereq
);

exit;


#
# subroutines
#

sub get_scripts_to_install {
  print "\n*** Script Install Section ****\n";
  print "\nBioperl comes with a number of useful scripts which you may wish to install.\n";
  my $prompt_mode = prompt("Install [a]ll Bioperl scripts, [n]one, or choose groups [i]nteractively?",'a');
  return if $prompt_mode =~ /^[nN]/;
  my @scripts;
  my $script_dir = catfile('./', 'scripts');
  my $dest_dir    = catfile('./', 'scripts_temp');
  rmtree($dest_dir) if -e $dest_dir;
  mkdir($dest_dir) or die "Can't mkdir($dest_dir): $!";
  opendir(F,$script_dir) or die "Can't opendir $script_dir: $!";
  while (my $file_or_dir = readdir(F)) {
    next if $file_or_dir =~ /^\./;
    next if $file_or_dir eq 'CVS';
    my $path = catfile($script_dir, $file_or_dir);
    next unless -d $path;
    next unless prompt_to_install($script_dir, $file_or_dir, $prompt_mode);
    my @files = install_contents($path, $dest_dir);
    push @scripts,@files;
  }
  closedir F;
  return @scripts;
}

sub prompt_to_install {
  my $prefix = shift;
  my $f      = shift;
  my $prompt_mode = shift;
  return 1 if $prompt_mode =~ /^a/i;  # all
  print "\n* Script Directory $f *\n";
  my $path = catfile($prefix, $f, 'TAG');
  if (-e $path && (my $g = IO::File->new($path))) {
    print while <$g>;
  }
  my $result = prompt("Install scripts in $f? [y]es [n]o [q]uit",'n');
  die if $result =~ /^q/i;
  return $result =~ /^y/i;
}

sub install_contents {
  my $dir  = shift;
  my $dest = shift;
  my $bangline = $Config{startperl};

  my @files;
  opendir (D,$dir) or die "Can't open $dir: $!\n";
  while (my $script = readdir(D)) {
    next unless $script =~ /\.PLS$/;
    my $path = catfile($dir, $script);
    my $in  = IO::File->new($path) or die "Can't open $path: $!";
    $script =~ s/\.PLS$/\.pl/;                   # change from .PLS to .pl
    $script =~ s/^/bp_/ unless $script =~ /^bp/; # add the "bp" prefix
    print "\tActivating $script....\n";
    my $path2 = catfile($dest, $script);
    unlink $path2 if -e $path2;
    my $out = IO::File->new(">$path2")  or die "Can't open $path2: $!";
    my $doneit;
    while (<$in>) {
      next if $doneit;
      if (s/^\#\!\S+/$bangline/) {
        $doneit++;
      }
    } continue {
      print $out $_;
    }
    close $in;
    close $out;
    chmod 0555, $path2 or die "Can't change mode of $script to ",MODE(),": $!";
    push @files,$path2;
  }
  closedir D;
  return @files;
}

sub prompt_for_biodbgff {
  if (eval { require DBI } && grep {/mysql/} DBI->available_drivers) {
    print "\n*** Bio::DB::GFF Tests section *** \n";
    my $proceed = prompt("Do you want to run the BioDBGFF tests (requires a mysql, postgres or oracle database) y/n",'n');
    my $path = catfile('t', 'do_biodbgff.tests');
    if ($proceed =~ /^y/i) {
      my $cfg = {dbd_driver => 'mysql'};
      GFFQueryDb($cfg,{'prompt'=>1,'verbose'=>1});
      if (open T,">$path") {
        while (my ($key, $value) = each %$cfg) {
          print T $key, "\t", $value,"\n";
        }
        close T;
      }
      push(@should_be_cleaned, $path);
    }
    else {
      if (-e $path ) {
        unlink $path;
      }
    }
  }
}

# creates $cfg->{$driver}->{test_(db|host|dsn|user|pass)}
sub GFFQueryDb {
  my $cfg     = shift;
  my $options = shift || {};
  my $db     = $cfg->{'description'} || $cfg->{'dbd_driver'};
  my $driver = $cfg->{'dbd_driver'};

  my $prompt = $options->{'prompt'};

  $driver = prompt("Which database driver should I use for testing BioDBGFF (mysql,Pg,Oracle)?",
            $driver) if $prompt;

  my $d = lc $driver;
  $cfg->{'dbd_driver'} = $driver;

  my $test_db =  exists($options->{"$d-test-db"}) ?
    $options->{"$d-test-db"} : ($cfg->{'test_db'} || 'test');

  $test_db = prompt
    ("Which database should I use for testing the $driver driver?",
     $test_db) if $prompt;

  my $test_host = exists($options->{"$d-test-host"}) ?
    $options->{"$d-test-host"} : ($cfg->{'test_host'} || 'localhost');
  $test_host = prompt
    ("On which host is database $test_db running (hostname, ip address or host:port)", $test_host) if $prompt;

  my($test_user, $test_pass);

  $test_user = exists($options->{"$d-test-user"}) ?
    $options->{"$d-test-user"} : ($cfg->{'test_user'} || "undef");
  $test_user = prompt
    ("User name for connecting to database $test_db?", $test_user)
      if $prompt;
  $test_user = undef if $test_user eq 'undef';

  $test_pass = exists($options->{"$d-test-pass"}) ?
    $options->{"$d-test-pass"} : ($cfg->{'test_pass'} || "undef");
  $test_pass = prompt
    ("Password for connecting to database $test_db?", $test_pass)
      if $prompt;
  $test_pass = undef if $test_pass eq 'undef';

  $cfg->{'test_db'}   = $test_db;
  $cfg->{'test_host'} = $test_host;
  $cfg->{'test_user'} = $test_user;
  $cfg->{'test_pass'} = $test_pass;
  if ($test_host eq 'undef'  ||  $test_host eq 'localhost') {
    $test_host = '';
  }
  my $test_dsn;
  if ($d eq 'pg') {
    $test_dsn = "dbi:$driver:dbname=$test_db";
  } else {
    $test_dsn = "dbi:$driver:database=$test_db";
  }
  $cfg->{test_dsn} = $test_dsn;
  if ($test_host) {
    $cfg->{'test_dsn'} .= ";host=$test_host";
  }

  if ($options->{'verbose'}) {
    local $^W=0;  # some uninitialized variable warning coming through
    print("Driver $driver is using the following settings for tests:\n",
      "    Database $test_db\n",
      "    Host     $test_host\n",
      "    DSN      $test_dsn\n",
      "    User     $test_user\n",
      "    Password $test_pass\n");
  }
}

sub setup_bioseqfeaturestore {
  print "\n*** Bio::DB::SeqFeature Tests section *** \n";
  my $path0 = catfile('t', 'BioDBSeqFeature.t');
  if (eval "require DB_File;") {
    print "I see you have Berkeleydb installed. I will create the BDB tests for Bio::DB::SeqFeature::Store...\n";
    my $path = catfile('t', 'BioDBSeqFeature_BDB.t');
    open my $F, ">$path";
    print $F <<END;
system 'perl $path0 -adaptor berkeleydb -create 1 -temp 1';
END
    close $F;
    push(@should_be_cleaned, $path);
  }
  else {
    print "I saw that you didn't have Berkelydb installed, so I didn't create BDB tests\n";
  }

  if (eval { require DBI } && grep {/mysql/} DBI->available_drivers) {
    if (DBI->connect('dbi:mysql:test',undef,undef,{RaiseError=>0,PrintError=>0})) {
      warn "I see you have mysql running and have a connectable 'test' database. I will create the mysql tests for Bio::DB::SeqFeature::Store...\n";
      my $path = catfile('t', 'BioDBSeqFeature_mysql.t');
      open my $F,">$path";
      print $F <<END;
system 'perl $path0 -adaptor DBI::mysql -create 1 -temp 1 -dsn test';
END
      close $F;
      push(@should_be_cleaned, $path);
    }
    else {
      print "I saw that you didn't have a connectable 'test' database, so I didn't create corresponding tests\n";
    }
  }
  else {
    print "I saw that you had no mysql DBI driver, so I didn't create corresponding tests\n";
  }
}


#
# override MakeMaker methods
#

sub MY::postamble {
  my $path = catfile('maintenance', 'symlink_scripts.pl');
  $tline . "\n\n"
         . "link_install ::\n\t\$(PERL) '".cwd()."/".$path."' '\$(INSTALLSITEBIN)'\n\n";
}

package MY;

sub MY::clean {
  return shift->SUPER::clean(@_) . <<END;
\t-rm -rf scripts_temp @should_be_cleaned
END
}

sub MY::manifypods {
  my $self = shift;
  if (1) {
    return "\nmanifypods : pure_all\n\t$self->{NOECHO}\$(NOOP)\n"
  }
  else {
    return $self->SUPER::manifypods(@_);
  }
}

sub MY::install {
  my $inherited = shift->SUPER::install(@_);
  $inherited =~ s/doc_install/doc_install link_install /;
  $inherited;
}
