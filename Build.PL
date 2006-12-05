#!/usr/bin/perl -w

# This is a Module::Build script for Bioperl installation.
# See http://search.cpan.org/~kwilliams/Module-Build/lib/Module/Build.pm

# Uses a custom subclass of Module::Build called ModuleBuildBioperl that
# doesn't get installed

# In the future developers may need to alter the requires and recommends and
# possibly Network sections of ModuleBuildBioperl->new() below, but otherwise
# nothing else here is likely to need changing.

use strict;
use ModuleBuildBioperl;

our @drivers;

# Set up the ModuleBuildBioperl object
my $build = ModuleBuildBioperl->new(
    module_name         => 'Bio',
    dist_name           => 'bioperl',
    dist_version_from   => 'Bio/Root/Version.pm',
    dist_author         => 'Bioperl Team <bioperl-l@bioperl.org>',
    dist_abstract       => 'Bioinformatics Toolkit',
    license             => 'artistic',
    requires            => {
                            'perl'                      => '5.6.1',
                            'IO::String'                => 0,
                            'DB_File'                   => 0
                           },
    build_requires      => {
                            'Test::More'                => 0,
                            'Module::Build'             => 0.2805,
                            'Test::Harness'             => 2.62,
                            'CPAN'                      => 1.81
                           },
    recommends          => [ # does what you would expect of recommends, except more informative output and generates optional_features in META.yml
                             # also note this is an array of hash refs, as the order is important
                            {'Ace'                       => '0/access of ACeDB database/Bio::DB::Ace,Bio::DB::GFF::Adaptor::ace'},
                            {'Bio::ASN1::EntrezGene'     => '0/parsing entrezgene/Bio::SeqIO::entrezgene'},
                            # we actually need 1.01 of Class::AutoClass, but unfortunately it is versioned as 1.0
                            {'Class::AutoClass'          => '1/creating objects/Bio::Graph::SimpleGraph,Bio::Graph::SimpleGraph::Traversal,Bio::Graph::ProteinGraph'},
                            {'Clone'                     => '0/cloning objects/Bio::Graph::ProteinGraph,Bio::Tools::Primer3'},
                            {'Convert::Binary::C'        => '0/strider functionality/Bio::SeqIO::strider'},
                            {'Data::Stag::XMLWriter'     => '0/writing choas xml files/Bio::SeqIO::chaosxml'},
                            {'GD'                        => '1.3/rendering Sequences and Features/Bio::Graphics::Glyph*'},
                            {'GD::SVG'                   => '0/creating SVG images/Bio::Graphics::Panel'},
                            # we specifically want Graph::Directed, but that has no VERSION
                            {'Graph'                     => '0/ontology engine implementation for the GO parser/Bio::Ontology::SimpleGOEngine::GraphAdaptor'},
                            {'HTML::Entities'            => '0/remote analysis POST submissions/Bio::SearchIO::blastxml'},
                            {'HTML::Parser'              => '3/screen scraping www.gdb.org/Bio::DB::GDB'},
                            {'HTTP::Request::Common'     => '0/GenBank+GenPept sequence retrieval, remote http Blast jobs/Bio::DB::*,Bio::Tools::Run::RemoteBlast,Bio::Tools::Analysis::Protein*,Bio::Tools::Analysis::DNA*'},
                            {'LWP::UserAgent'            => '0/remote access/Bio::DB::*,Bio::Tools::Run::RemoteBlast,Bio::Tools::WebBlat,Bio::WebAgent,Bio::Graphics::Glyph::image'},
                            {'Set::Scalar'               => '0/proper operation/Bio::Tree::Compatible'},
                            {'SOAP::Lite'                => '0/XEMBL Services and Bibliographic queries/Bio::DB::XEMBLService,Bio::DB::Biblio::soap'},
                            {'Spreadsheet::ParseExcel'   => '0/parsing Excel files/Bio::SeqIO::excel'},
                            {'Storable'                  => '0/storing sequence objects in local file cache/Bio::DB::FileCache,Bio::SeqFeature::Collection,Bio::PopGen::HtSNP,Bio::PopGen::TagHaplotype,Bio::DB::GFF::Adaptor::berkeleydb'},
                            {'SVG'                       => '2.26/SVG output/Bio::Graphics::Pictogram'},
                            {'SVG::Graph'                => '0.01/creating SVG images/Bio::TreeIO::svggraph'},
                            {'Text::Shellwords'          => '0/test scripts/Bio::Graphics'},
                            {'URI::Escape'               => '0/dealing with web resources/Bio::Tools::WebBlat,Bio::FeatureIO::gff,Bio::FeatureIO::interpro,Bio::DB::Biblio::eutils,Bio::DB::EUtilities::Cookie,Bio::DB::Query::GenBank,Bio::DB::NCBIHelper,Bio::SeqFeature::Annotated'},
                            # this is actually just an undeclared req of XML::DOM::XPath, used by Bio::FeatureIO::interpro
                            {'XML::XPath'                => '0/parsing xml/XML::DOM::XPath,Bio::FeatureIO::interpro'},
                            {'XML::DOM::XPath'           => '0/parsing interpro features/Bio::FeatureIO::interpro'},
                            {'XML::Parser'               => '0/parsing xml/Bio::Biblio::IO::medlinexml'},
                            {'XML::Parser::PerlSAX'      => '0/parsing xml/Bio::SeqIO::tinyseq,Bio::SeqIO::game::gameSubs,Bio::OntologyIO::InterProParser,Bio::ClusterIO::dbsnp'},
                            {'XML::SAX'                  => '0/parsing xml/Bio::SearchIO::blastxml,Bio::SeqIO::tigrxml,Bio::SeqIO::bsml_sax'},
                            {'XML::SAX::Writer'          => '0/writing xml/Bio::SeqIO::tigrxml'},
                            # this is actually just an undeclared req of XML::Simple, used by EUtilities
                            {'XML::SAX::ExpatXS'         => '0/parsing xml/XML::Simple,Bio::DB::EUtilities'},
                            {'XML::Simple'               => '0/parsing xml/Bio::DB::EUtilities'},
                            {'XML::Twig'                 => '0/parsing xml/Bio::Variation::IO::xml,Bio::DB::Taxonomy::entrez,Bio::DB::Biblio::eutils,Bio::Graph::IO::psi_xml'},
                            {'XML::Writer'               => '0.4/parsing and writing xml/Bio::SeqIO::agave,Bio::SeqIO::game::gameWriter,Bio::SeqIO::chadoxml,Bio::SeqIO::tinyseq,Bio::Variation::IO::xml,Bio::SearchIO::Writer::BSMLResultWriter'}
                           ],
    get_options         => {
                            network => { } # not actually used by anything yet, but in the future say perl Build.PL --network
                           },
    auto_features       => {
                            BioDBSeqFeature_BDB   => {
                                                        description      => "BDB tests for Bio::DB::SeqFeature::Store",
                                                        feature_requires => { 'DB_File' => 0 } # feature_requires is like requires, execpt that it doesn't trigger installation
                                                     },
                            BioDBGFF              => {
                                                        description      => "BioDBGFF database tests (will need to answer questions before really enabling)",
                                                        feature_requires => { 'DBI' => 0 },
                                                        excludes_os      => ['mswin'],
                                                        test             => \&test_biodbgff # ModuleBuildBioperl unique requirement that after everything else succeeds, supplied code ref must also return undef
                                                     },
                            BioDBSeqFeature_mysql => {
                                                        description      => "MySQL tests for Bio::DB::SeqFeature::Store",
                                                        feature_requires => { 'DBI' => 0, 'DBD::mysql' => 0 },
                                                        test             => \&test_db
                                                     },
# The following code works, but since no tests in the test suite actually make use of this functionality, don't use it yet
#                            Network               => {
#                                                        description => "Enable tests that need an internet connection",
#                                                        requires    => { 'LWP::UserAgent' => 0 },
#                                                        options     => ['network'], # ModuleBuildBioperl unique requirement that --network was supplied
#                                                        test        => \&ModuleBuildBioperl::test_internet
#                                                     }
                           },
    dynamic_config      => 1
    
    #pm_files           => {} # modules in Bio are treated as if they were in lib and auto-installed
    #script_files       => [] # scripts in scripts directory are installed on-demand
);

# Handle auto features
if ($build->feature('BioDBSeqFeature_BDB')) {
    make_bdb_test();
}
if ($build->feature('BioDBSeqFeature_mysql')) {
    make_dbi_test();
}
{
    $build->notes(network => $build->feature('Network'));
    # then in test script:
    #   use Module::Build;
    #   my $build = Module::Build->current;
    #   my $do_network_tests = $build->notes('network');
}

# Ask questions
$build->choose_scripts;
prompt_for_biodbgff() if $build->feature('BioDBGFF');

# Request that some scripts run post-installation
$build->add_post_install_script('maintenance/symlink_script.pl'); # takes a unix file path regardless of local OS

# Add extra things to MANIFEST.SKIP
$build->add_to_manifest_skip('bioperl.lisp', 'Bio/Tools/WebBlat.pm');

# Create the build script and exit
$build->create_build_script;

exit;


sub make_bdb_test {
    my $path0 = File::Spec->catfile('t', 'BioDBSeqFeature.t');
    my $path = File::Spec->catfile('t', 'BioDBSeqFeature_BDB.t');
    open my $F, ">$path";
    print $F <<END;
system 'perl $path0 -adaptor berkeleydb -create 1 -temp 1';
END
    close $F;
    $build->add_to_cleanup($path);
    $build->add_to_manifest_skip($path);
}

sub test_db {
    eval {require DBI;};  # if not installed, this sub won't actually be called
    unless (eval {DBI->connect('dbi:mysql:test',undef,undef,{RaiseError=>0,PrintError=>0})}) {
        return "Could not connect to test database";
    }
    return;
}

sub make_dbi_test {
    my $path0 = File::Spec->catfile('t', 'BioDBSeqFeature.t');
    my $path = File::Spec->catfile('t', 'BioDBSeqFeature_mysql.t');
    open my $F,">$path";
    print $F <<END;
system 'perl $path0 -adaptor DBI::mysql -create 1 -temp 1 -dsn test';
END
    close $F;
    $build->add_to_cleanup($path);
    $build->add_to_manifest_skip($path);
}

sub test_biodbgff {
    eval {require DBI;};  # if not installed, this sub won't actually be called
    @drivers = DBI->available_drivers;
    unless (grep {/mysql|Pg|Oracle/i} @drivers) {
        return "MySQL, Pg nor Oracle DBI drivers are installed";
    }
    return;
}

sub prompt_for_biodbgff {
    my $proceed = $build->y_n("Do you want to run the BioDBGFF live database tests? y/n", 'n');
    
    if ($proceed) {
        my @driver_choices;
        foreach my $poss ('mysql', 'Pg', 'Oracle') {
            if (grep {/$poss/i} @drivers) {
                my $choice = $poss;
                $choice =~ s/^(.)/[$1]/;
                push(@driver_choices, $choice);
            }
        }
        
        my $cfg = {};
        
        my $driver;
        if (@driver_choices > 1) {
            my ($default) = $driver_choices[0] =~ /\[(.)/;
            $driver = $build->prompt("Which database driver should be used? ".join(" ", @driver_choices), $default);
        }
        else {
            ($driver) = $driver_choices[0] =~ /\[(.)/;
        }
        if ($driver =~ /^[mM]/) {
            $driver = 'mysql';
        }
        elsif ($driver =~ /^[pP]/) {
            $driver = 'Pg';
        }
        elsif ($driver =~ /^[oO]/) {
            $driver = 'Oracle';
        }
        $cfg->{dbd_driver} = $driver;
        
        $cfg->{test_db} = $build->prompt("Which database should I use for testing the $driver driver?", 'test');
        $cfg->{test_host} = $build->prompt("On which host is database $cfg->{test_db} running (hostname, ip address or host:port)", 'localhost');
        my $test_user = $build->prompt("User name for connecting to database $cfg->{test_db}?", 'undef');
        $cfg->{test_user} = $test_user eq 'undef' ? '' : $test_user;
        my $test_pass = $build->prompt("Password for connecting to database $cfg->{test_db}?", 'undef');
        $cfg->{test_pass} = $test_pass eq 'undef' ? '' : $test_pass;
        
        my $use_host = 1;
        if ($cfg->{test_host} eq 'undef' || $cfg->{test_host} eq 'localhost') {
            $use_host = 0;
        }
        
        my $test_dsn;
        if ($driver eq 'Pg') {
            $test_dsn = "dbi:$driver:dbname=$cfg->{test_db}";
        }
        else {
            $test_dsn = "dbi:$driver:database=$cfg->{test_db}";
        }
        if ($use_host) {
            $test_dsn .= ";host=$cfg->{test_host}";
        }
        $cfg->{test_dsn} = $test_dsn;
        
        my $path = File::Spec->catfile('t', 'do_biodbgff.tests');
        if (open T,">$path") {
            while (my ($key, $value) = each %$cfg) {
                print T $key, "\t", $value,"\n";
            }
            close T;
        }
        
        $build->add_to_cleanup($path);
        $build->add_to_manifest_skip($path);
        
        $build->log_info("  - will run the BioDBGFF tests with database driver '$driver' and these settings:\n",
                        "    Database $cfg->{test_db}\n",
                        "    Host     $cfg->{test_host}\n",
                        "    DSN      $test_dsn\n",
                        "    User     $test_user\n",
                        "    Password $test_pass\n");
    }
    else {
        $build->log_info("  - will not run the BioDBGFF live database tests\n");
    }
    
    $build->log_info("\n");
}
