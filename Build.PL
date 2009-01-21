#!/usr/bin/perl -w

# This is a Module::Build script for Bioperl installation.
# See http://search.cpan.org/~kwilliams/Module-Build/lib/Module/Build.pm

# Uses a custom subclass of Module::Build called Bio::Root::Build

# In the future developers may need to alter the requires and recommends
# sections of Bio::Root::Build->new() below, but otherwise nothing else here is
# likely to need changing.

use strict;
use lib '.';
use Bio::Root::Build;

our @drivers;

my $mysql_ok = 0;

# Set up the Bio::Root::Build object
my $build = Bio::Root::Build->new(
    module_name         => 'Bio',
    dist_name           => 'BioPerl',
    dist_version_from   => 'Bio/Root/Version.pm',
    dist_author         => 'BioPerl Team <bioperl-l@bioperl.org>',
    dist_abstract       => 'Bioinformatics Toolkit',
    license             => 'perl',
    requires            => {
                            'perl'                      => '5.6.1',
                            'IO::String'                => 0,
                            'DB_File'                   => 0,
                            'Data::Stag'                => 0.11, # Bio::SeqIO::swiss, we can change to 'recommend' if needed
                            'Scalar::Util'              => 0 # not in Perl 5.6.1, arrived in core in 5.7.3
                           },
    build_requires      => {
                            'Test::More'                => 0,
                            'Module::Build'             => 0.2805,
                            'Test::Harness'             => 2.62,
                            'CPAN'                      => 1.81
                           },
    recommends          => { # does what you would expect of recommends, except more informative output and generates optional_features in META.yml
                            'Ace'                       => '0/access of ACeDB database/Bio::DB::Ace,Bio::DB::GFF::Adaptor::ace',
                            'Algorithm::Munkres'        => '0/Phylogenetic Networks/Bio::PhyloNetwork',
                            'Array::Compare'            => '0/Phylogenetic Networks/Bio::PhyloNetwork',
                            # this won't actually install due to circular dep, but we have no way of doing a post-install
                            # the [circular dependency!] specifies it is only installed on explicit request for this specific module,
                            # not when simply choosing to install 'all' modules
                            'Bio::ASN1::EntrezGene'     => '0/parsing entrezgene/Bio::SeqIO::entrezgene [circular dependency!]',
                            # we actually need 1.01 of Class::AutoClass, but unfortunately it is versioned as 1.0
                            'Clone'                     => '0/cloning objects/Bio::Tools::Primer3',
                            'Convert::Binary::C'        => '0/strider functionality/Bio::SeqIO::strider',
                            # we specifically want Graph::Directed, but that has no VERSION
                            'Graph'                     => '0/ontology engine implementation for the GO parser/Bio::PhyloNetwork',
                            'GraphViz'                  => '0/Phylogenetic Network Visulization/Bio::PhyloNetwork::GraphViz',
                            'HTML::Entities'            => '0/remote analysis POST submissions/Bio::SearchIO::blastxml',
                            'HTML::HeadParser'          => '3/parsing <HEAD> section of HTML docs/Bio::Tools::Analysis::DNA::ESEfinder',
                            'HTTP::Request::Common'     => '0/GenBank+GenPept sequence retrieval, remote http Blast jobs/Bio::DB::*,Bio::Tools::Run::RemoteBlast,Bio::Tools::Analysis::Protein*,Bio::Tools::Analysis::DNA*',
                            'LWP::UserAgent'            => '0/remote access/Bio::DB::*,Bio::Tools::Run::RemoteBlast,Bio::WebAgent',
                            'Math::Random'              => '0/Random Phylogenetic Networks/Bio::PhyloNetwork::RandomFactory',
                            'PostScript::TextBlock'     => '0/EPS output/Bio::Tree::Draw::Cladogram',
                            'Set::Scalar'               => '0/proper operation/Bio::Tree::Compatible',
                            'SOAP::Lite'                => '0/Bibliographic queries/Bio::DB::Biblio::soap',
                            'Spreadsheet::ParseExcel'   => '0/parsing Excel files/Bio::SeqIO::excel',
                            'Storable'                  => '2.05/storing sequence objects in local file cache/Bio::DB::FileCache,Bio::SeqFeature::Collection,Bio::PopGen::HtSNP,Bio::PopGen::TagHaplotype,Bio::DB::GFF::Adaptor::berkeleydb',
                            'SVG::Graph'                => '0.01/creating SVG images/Bio::TreeIO::svggraph',
                            'Text::ParseWords'          => '0/test scripts/Bio::DB::SeqFeature::Store::FeatureFileLoader',
                            'URI::Escape'               => '0/dealing with web resources/Bio::FeatureIO::gff,Bio::FeatureIO::interpro,Bio::DB::Biblio::eutils,Bio::DB::EUtilParameters,Bio::DB::Query::GenBank,Bio::DB::NCBIHelper,Bio::SeqFeature::Annotated',
                            'XML::DOM::XPath'           => '0.13/parsing interpro features/Bio::FeatureIO::interpro',
                            'XML::Parser'               => '0/parsing xml/Bio::Biblio::IO::medlinexml',
                            'XML::Parser::PerlSAX'      => '0/parsing xml/Bio::SeqIO::tinyseq,Bio::SeqIO::game::gameSubs,Bio::OntologyIO::InterProParser,Bio::ClusterIO::dbsnp',
                            'XML::SAX'                  => '0.15/parsing xml/Bio::SearchIO::blastxml,Bio::SeqIO::tigrxml,Bio::SeqIO::bsml_sax',
                            'XML::SAX::Writer'          => '0/writing xml/Bio::SeqIO::tigrxml',
                            'XML::Twig'                 => '0/parsing xml/Bio::Variation::IO::xml,Bio::DB::Taxonomy::entrez,Bio::DB::Biblio::eutils',
                            'XML::Writer'               => '0.4/parsing and writing xml/Bio::SeqIO::agave,Bio::SeqIO::game::gameWriter,Bio::SeqIO::chadoxml,Bio::SeqIO::tinyseq,Bio::Variation::IO::xml,Bio::SearchIO::Writer::BSMLResultWriter',
    },
    get_options         => {
                            network => { } # say 'perl Build.PL --network' to manually request network tests
                           },
    auto_features       => {
                            BioDBSeqFeature_BDB   => {
                                                        description      => "BDB tests for Bio::DB::SeqFeature::Store",
                                                        feature_requires => { 'DB_File' => 0 } # feature_requires is like requires, except that it doesn't trigger installation
                                                     },
                            BioDBGFF              => {
                                                        description      => "BioDBGFF database tests (will need to answer questions before really enabling)",
                                                        feature_requires => { 'DBI' => 0 },
                                                        excludes_os      => ['mswin'],
                                                        test             => \&test_biodbgff # Bio::Root::Build unique requirement that after everything else succeeds, supplied code ref must also return undef
                                                     },
                            BioDBSeqFeature_mysql => {
                                                        description      => "MySQL tests for Bio::DB::SeqFeature::Store",
                                                        feature_requires => { 'DBI' => 0, 'DBD::mysql' => 0 },
                                                        test             => \&test_db_sf
                                                     },
                            Network               => {
                                                        description => "Enable tests that need an internet connection",
                                                        requires    => { 'LWP::UserAgent' => 0 },
                                                        test        => \&Bio::Root::Build::test_internet
                                                     }
                           },
    dynamic_config      => 1,
    create_makefile_pl  => 'passthrough',
    recursive_test_files => 1,
    # Extra files needed for BioPerl modules
    xml_files => {'./Bio/DB/HIV/lanl-schema.xml' => 'lib/Bio/DB/HIV/lanl-schema.xml'},
    
    #pm_files           => {} # modules in Bio are treated as if they were in lib and auto-installed
    #script_files       => [] # scripts in scripts directory are installed on-demand
);

my $accept = $build->args->{accept};

my $proceed = prompt_for_biodb($accept) if $build->feature('BioDBGFF') || $build->feature('BioDBSeqFeature_mysql');

# Handle auto features
if ($proceed && $build->feature('BioDBSeqFeature_BDB')) {
    # will return without doing anything if user chose not to run tests during
    # prompt_for_biodb() above
    make_bdb_test();
}
if ($proceed && ($build->feature('BioDBSeqFeature_mysql'))) {
    make_dbi_test();
}

# Ask questions
$build->choose_scripts($accept);
#prompt_for_biodbgff($accept) if $build->feature('BioDBGFF');
{
    if ($build->args('network')) {
        if ($build->feature('Network')) {
            $build->notes(network => 1);
            $build->log_info("  - will run internet-requiring tests\n");
        }
        else {
            $build->notes(network => 0);
            $build->log_info("  - will not run network tests since I seem to be missing essential network functionality\n");
        }
    }
    else {
        $build->prompt_for_network($accept) if $build->feature('Network');
    }
    # then in test script:
    #   use Bio::Root::Build;
    #   my $build = Module::Build->current;
    #   my $do_network_tests = $build->notes('network');
}

# Request that some scripts run post-installation
$build->add_post_install_script('maintenance/symlink_script.pl'); # takes a unix file path regardless of local OS

# Add extra things to MANIFEST.SKIP
$build->add_to_manifest_skip('bioperl.lisp');

# Add additional files here
$build->add_build_element('xml');

# Create the build script and exit
$build->create_build_script;

exit;

sub make_bdb_test {
    my $path0 = File::Spec->catfile('t', 'LocalDB', 'SeqFeature.t');
    my $path = File::Spec->catfile('t', 'LocalDB','SeqFeature_BDB.t');
    unlink($path) if (-e $path); 
    open(my $F, ">", $path) || die "Can't create test file\n";
    print $F <<END;
system '$^X $path0 -adaptor berkeleydb -create 1 -temp 1';
END
    close $F;
    $build->add_to_cleanup($path);
    $build->add_to_manifest_skip($path);
}

sub test_db_sf {
    eval {require DBI;};  # if not installed, this sub won't actually be called
    @drivers = DBI->available_drivers;
    unless (grep {/mysql/i} @drivers) {
        $mysql_ok = 0;
        return "Only MySQL driver is supported for Bio::DB::SeqFeature RDMS tests";
    }
    $mysql_ok = 1;
    return;
}

sub make_dbi_test {
    my $dsn = $build->notes('test_dsn') || return;
    my $path0 = File::Spec->catfile('t', 'LocalDB', 'SeqFeature.t');
    my $driver = $build->notes('dbd_driver');
    my $path = File::Spec->catfile('t', 'LocalDB', 'SeqFeature_mysql.t' );
    my $test_db = $build->notes('test_db');
    my $user = $build->notes('test_user');
    my $pass = $build->notes('test_pass');
    open my $F,">$path";
    my $str = "$path0 -adaptor DBI::$driver -create 1 -temp 1 -dsn $dsn";
    $str .= " -user $user" if $user;
    $str .= " -password $pass" if $pass;
    print $F <<END;
system '$^X $str';
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

sub prompt_for_biodb {
    my $accept = shift;
    my $proceed = $accept ? 0 : $build->y_n("Do you want to run the Bio::DB::GFF or ".
                              "Bio::DB::SeqFeature::Store live database tests? ".
                              "y/n", 'n');
    
    if ($proceed) {
        my @driver_choices;
        foreach my $poss ('mysql', 'Pg', 'Oracle') {
            if (grep {/$poss/i} @drivers) {
                my $choice = $poss;
                $choice =~ s/^(.)/[$1]/;
                push(@driver_choices, $choice);
            }
        }
        
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
        
        my $test_db = $build->prompt("Which database should I use for testing the $driver driver?\n".
                                     "This database should already be present but doesn't have to ".
                                     "be preloaded for any schema", 'test');
        my $test_host = $build->prompt("On which host is database '$test_db' running (hostname, ip address or host:port)", 'localhost');
        my $test_user = $build->prompt("User name for connecting to database '$test_db'?", 'undef');
        my $test_pass = $build->prompt("Password for connecting to database '$test_db'?", 'undef');
        
        my $use_host = 1;
        if ($test_host eq 'undef' || $test_host eq 'localhost') {
            $use_host = 0;
        }
        
        my $test_dsn;
        if ($driver eq 'Pg') {
            $test_dsn = "dbi:$driver:dbname=$test_db";
            $mysql_ok = 0;
        }
        else {
            $test_dsn = "dbi:$driver:database=$test_db";
            $mysql_ok = 0;
        }
        if ($use_host) {
            $test_dsn .= ";host=$test_host";
        }
        
        $build->notes(dbd_driver => $driver);
        $build->notes(test_db => $test_db);
        $build->notes(test_host => $test_host);
        $build->notes(test_user => $test_user eq 'undef' ? undef : $test_user);
        $build->notes(test_pass => $test_pass eq 'undef' ? undef : $test_pass);
        $build->notes(test_dsn => $test_dsn);
        
        $build->log_info("  - will run tests with database driver '$driver' and these settings:\n",
                        "    Database $test_db\n",
                        "    Host     $test_host\n",
                        "    DSN      $test_dsn\n",
                        "    User     $test_user\n",
                        "    Password $test_pass\n");
        $build->log_info("  - will not run the BioDBSeqFeature live ".
                         "database tests (requires MySQL driver)\n") unless ($driver eq 'mysql');
    }
    else {
        $build->log_info("  - will not run the BioDBGFF or BioDBSeqFeature live database tests\n");
    }
    
    $build->log_info("\n");
    return $proceed;
}
