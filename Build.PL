#!/usr/bin/perl

# This is a Module::Build script for Bioperl installation.
# See http://search.cpan.org/~kwilliams/Module-Build/lib/Module/Build.pm

# Uses a custom subclass of Module::Build called Bio::Root::Build

# In the future developers may need to alter the requires and recommends
# sections of Bio::Root::Build->new() below, but otherwise nothing else here is
# likely to need changing.

use strict;
use lib '.';
use Bio::Root::Build;

# XML::SAX::RTF doesn't work with BioPerl, at all, nada, zilch.
#
# Since we're running into this now on CPAN Testers, catch it up front and
# deal with it.
#
# See: https://rt.cpan.org/Ticket/Display.html?id=5943
#      https://redmine.open-bio.org/issues/2975

{ eval { require XML::SAX; 1; };

unless ($@) {
    if (grep {$_->{Name} =~ 'XML::SAX::RTF'} @{XML::SAX->parsers()}) {
        warn <<WARN;

############################# WARNING #############################

XML::SAX::RTF is not XML::SAX-compliant but is registered as an
XML::SAX parser.  If used as the primary parser, modules requiring
XML::SAX will NOT work. Please install another XML::SAX-compliant
module and modify your local ParserDetails.ini file per XML::SAX
docs to remove references to XML::SAX::RTF.

############################# WARNING #############################

WARN
    sleep 2;
    }
}

}

my %recommends = (
    # AcePerl support is deprecated, per LDS - cjfields - 5-5-2011

    #'Ace'                       => [0,
    #    'Access of ACeDB database/Bio::DB::Ace,Bio::DB::GFF::Adaptor::ace'],

    'Algorithm::Munkres'        => [0,
        'Phylogenetic Networks/Bio::PhyloNetwork'],

    'Array::Compare'            => [0,
        'Phylogenetic Networks/Bio::PhyloNetwork'],

    # this won't actually install due to circular dep, but we have no way of
    # doing a post-install the [circular dependency!] specifies it is only
    # installed on explicit request for this specific module, not when simply
    # choosing to install 'all' modules

    #'Bio::ASN1::EntrezGene'     => [0,
    #    'Parsing entrezgene/Bio::SeqIO::entrezgene [circular dependency!]'],

    'Bio::Phylo'                => [0,
        'NeXML related modules/Bio::AlignIO::nexml,Bio::Nexml::Factory,'.
        'Bio::SeqIO::nexml,Bio::TreeIO::nexml'],

    'Clone'                     => [0,
        'Cloning objects/Bio::Root::Root,Bio::Tools::Primer3'],

    'Convert::Binary::C'        => [0,
        'Strider functionality/Bio::SeqIO::strider'],

    'DB_File'                   => [0,
        'Database functionality/Bio::Assemby,Bio::DB'],

    'Error'                     => [0,
        'OO-based exception handling (very optional)/Bio::Root::Exception'],

    'GD'                        => [0,
        'Alignment graphic output/Bio::Align::Graphics'],

    'Graph'                     => [0.50,
        'Phylogenetic Networks, ontology engine implementation, contig analysis'.
        '/Bio::PhyloNetwork,Bio::Ontology::SimpleGOEngine::GraphAdaptor,'.
        'Bio::Assembly::Tools::ContigSpectrum'],

    'GraphViz'                  => [0,
        'Phylogenetic Network Visualization/Bio::PhyloNetwork::GraphViz'],

    'HTML::Entities'            => [0,
        'Remote analysis POST submissions/Bio::SearchIO::blastxml'],

    'HTML::HeadParser'          => [3,
        'Parsing <HEAD> section of HTML docs/Bio::Tools::Analysis::DNA::ESEfinder'],

    'HTML::TableExtract'        => [0,
        'Parsing HTML tables/Bio::DB::SeqVersion::gi'],

    'HTTP::Request::Common'     => [0,
        'GenBank+GenPept sequence retrieval, remote http Blast jobs'.
        '/Bio::DB::*,Bio::Tools::Run::RemoteBlast,'.
        'Bio::Tools::Analysis::Protein*,Bio::Tools::Analysis::DNA*'],

    'IO::Scalar'               => [0,
        'Deal with non-seekable filehandles/Bio::Tools::GuessSeqFormat'],

    'List::MoreUtils'           => [0,
        'Back- or reverse-translation of sequences/'.
        'Bio::Tools::SeqPattern,Bio::Tools::SeqPattern::BackTranslate'],

    'LWP::UserAgent'            => [0,
        'Remote access/Bio::DB::*,Bio::Tools::Run::RemoteBlast,Bio::WebAgent'],

    'PostScript::TextBlock'     => [0,
        'EPS output/Bio::Tree::Draw::Cladogram'],

    'Set::Scalar'               => [0,
        'Proper operation/Bio::Tree::Compatible'],

    'SOAP::Lite'                => [0,
        'Bibliographic queries/Bio::DB::Biblio::soap'],

    'Sort::Naturally'           => [0,
        'Sort lexically, but sort numeral parts numerically/'.
        'Bio::Assembly::IO::ace,Bio::Assembly::IO::tigr'],

    'Spreadsheet::ParseExcel'   => [0,
        'Parsing Excel files/Bio::SeqIO::excel'],

    'Storable'                  => [2.05,
        'Storing sequence objects in local file cache/'.
        'Bio::DB::FileCache,Bio::SeqFeature::Collection,Bio::PopGen::HtSNP,'.
        'Bio::PopGen::TagHaplotype,Bio::DB::GFF::Adaptor::berkeleydb,Bio::Root::Root'],

    'SVG'                       => [2.26,
        'Creating SVG images/Bio::Draw::Pictogram'],

    'SVG::Graph'                => [0.01,
        'Creating SVG images/Bio::TreeIO::svggraph'],

    'Text::ParseWords'          => [0,
        'Test scripts/Bio::DB::SeqFeature::Store::FeatureFileLoader'],

    'XML::DOM'                  => [0,
        'Parsing XML/Bio::SeqIO::bsml,Bio::SeqIO::interpro'],

    'XML::DOM::XPath'           => [0,
        'Parsing XML/Bio::SeqIO::interpro'],

    'XML::LibXML'               => [0,
        'Parsing XML/Bio::SeqIO::seqxml,Bio::TreeIO::phyloxml'],

    'XML::Parser'               => [0,
        'Parsing XML/Bio::Biblio::IO::medlinexml'],

    'XML::Parser::PerlSAX'      => [0,
        'Parsing XML/Bio::SeqIO::tinyseq,Bio::SeqIO::game::gameSubs,',
        'Bio::OntologyIO::InterProParser,Bio::ClusterIO::dbsnp'],

    'XML::SAX'                  => [0.15,
        'Parsing XML/Bio::SearchIO::blastxml,Bio::SeqIO::tigrxml,Bio::SeqIO::bsml_sax'],

    'XML::SAX::Writer'          => [0,
        'Writing XML/Bio::SeqIO::tigrxml'],

    'XML::Simple'               => [0,
        'Reading custom XML/Bio::Tools::EUtilities,Bio::DB::HIV,Bio::DB::Query::HIVQuery'],

    'XML::Twig'                 => [0,
        'Parsing XML/Bio::Variation::IO::xml,Bio::DB::Taxonomy::entrez,'.
        'Bio::DB::Biblio::eutils'],

    'XML::Writer'               => [0.4,
        'Parsing and writing XML/Bio::SeqIO::agave,Bio::SeqIO::game::gameWriter,'.
        'Bio::SeqIO::chadoxml,Bio::SeqIO::tinyseq,Bio::Variation::IO::xml,'.
        'Bio::SearchIO::Writer::BSMLResultWriter'],

    'YAML'                      => [0,
        'GenBank->GFF3/bp_genbank2gff3.pl'],
);

my $mysql_ok = 0;

my @drivers = available_drivers();

# Set up the Bio::Root::Build object
my $build = Bio::Root::Build->new(
    module_name         => 'Bio',
    dist_name           => 'BioPerl',
    dist_version_from   => 'Bio/Root/Version.pm',
    dist_author         => 'BioPerl Team <bioperl-l@bioperl.org>',
    dist_abstract       => 'Bioinformatics Toolkit',
    license             => 'perl',
    no_index            => {'x_dir'       => [qw(examples/root/lib)]},
    requires            => {
                            'perl'                      => '5.6.1',
                            'IO::String'                => 0,    # why is this required?
                            'Data::Stag'                => 0.11, # Bio::SeqIO::swiss, we can change to 'recommend' if needed
                            'Scalar::Util'              => 0,    # not in Perl 5.6.1, arrived in core in 5.7.3
                            'ExtUtils::Manifest'        => '1.52', # allows spaces in file names
                           },

    build_requires      => {
                            'CPAN'                      => 1.81,
                            'Module::Build'             => 0.2805,
                            'Test::Harness'             => 2.62,
                            'Test::Most'                => 0,
                            'URI::Escape'               => 0
                           },

    recommends          => {
        # reverted to a simple Module::Build-compatible hash, but we keep
        # additional data in the %recommends hash above. May be converted to
        # something simpler if there aren't complaints down the line.
        map {$_ => $recommends{$_}[0]} sort keys %recommends
    },

    get_options         => {
        accept  => { },
        network => { } # say 'perl Build.PL --network' to manually request network tests
                           },

    auto_features       => {
        'EntrezGene'    => {
            description => "Presence of Bio::ASN1::EntrezGene",
            requires    => { 'Bio::ASN1::EntrezGene' => 0 } # feature_requires is like requires, except that it doesn't trigger installation
        },

        'DB_File Tests' => {
            description => "BDB tests for Bio::DB::SeqFeature::Store",
            requires    => { 'DB_File' => 0 } # feature_requires is like requires, except that it doesn't trigger installation
         },

        'Bio::DB::GFF Tests' => {
            description => "Bio::DB::GFF database tests (will need to answer questions before really enabling)",
            requires    => { 'DBI' => 0 },
         },

        'MySQL Tests'   => {
            description => "MySQL-related tests for Bio::DB::SeqFeature::Store",
            requires    => { 'DBI' => 0, 'DBD::mysql' => 0 },
         },

        'Pg Tests'      => {
            description => "PostgreSQL-related tests for Bio::DB::SeqFeature::Store",
            requires    => { 'DBI' => 0, 'DBD::Pg' => 0},
         },

        'SQLite Tests'  => {
            description => "SQLite-related tests for Bio::DB::SeqFeature::Store",
            requires    => { 'DBI' => 0, 'DBD::SQLite' => 0},
         },

        'Network Tests' => {
            description => "Enable tests that need an internet connection",
            requires    => { 'LWP::UserAgent' => 0 },
        }
                           },
    dynamic_config      => 1,
    #create_makefile_pl  => 'passthrough',
    recursive_test_files => 1,

    # Extra files needed for BioPerl modules
    xml_files => {'./Bio/DB/HIV/lanl-schema.xml' => 'lib/Bio/DB/HIV/lanl-schema.xml'},

    #pm_files           => {} # modules in Bio are treated as if they were in lib and auto-installed
    #script_files       => [] # scripts in scripts directory are installed on-demand
);

my $accept = $build->args('accept');

# how much do I hate this?  Let me count the ways.....
#if (!$build->feature('EntrezGene')) {
#    warn <<WARN;
############################## WARNING #############################
#
#Bio::ASN1::EntrezGene not found. This is an *optional* module;
#however, because it has a circular dependency with BioPerl we do not
#include it on our list of recommended modules.
#
#If you require EntrezGene functionality, you can install
#Bio::ASN1::EntrezGene after BioPerl has finished installing.
#
############################## WARNING #############################
#
#WARN
#    sleep 3;
#}

my $proceed = prompt_for_biodb($accept)
    if $build->feature('Bio::DB::GFF') || $build->feature('MySQL Tests') ||
      $build->feature('Pg Tests')      || $build->feature('SQLite Tests');

# Handle auto features
if ($proceed && $build->feature('DB_File Tests')) {
    # will return without doing anything if user chose not to run tests during
    make_bdb_test();
}
if ($proceed && ($build->feature('MySQL Tests') ||
                 $build->feature('Pg Tests') ||
                 $build->feature('SQLite Tests'))) {
    make_dbi_test();
}

# Ask questions
$build->choose_scripts($accept);

if ($build->args('network')) {
    if ($build->feature('Network Tests')) {
        $build->notes(network => 1);
        $build->log_info("  - will run internet-requiring tests\n");
    }
    else {
        $build->notes(network => 0);
        $build->log_info("  - Missing LWP::UserAgent, can't run network tests\n");
    }
}
else {
    $build->prompt_for_network($accept) if $build->feature('Network Tests');
}

# Add additional files here
$build->add_build_element('xml');

# Create the build script and exit
$build->create_build_script;

exit;

########################## Helper subs ##########################

sub make_bdb_test {
    my $path0 = File::Spec->catfile('t', 'LocalDB', 'SeqFeature.t');
    my $path = File::Spec->catfile('t', 'LocalDB','SeqFeature_BDB.t');
    unlink($path) if (-e $path);
    open my $F, '>', $path or die "Could not write test file '$path': $!\n";
    print $F <<END;
system '$^X $path0 -adaptor berkeleydb -create 1 -temp 1';
END
    close $F;
    $build->add_to_cleanup($path);
    #$build->add_to_manifest_skip($path);
}

sub available_drivers {
    eval {require DBI; 1;};  # if not installed, this sub won't actually be called
    return if $@;
    @drivers = DBI->available_drivers;
    unless (grep {/mysql|Pg|SQLite/i} @drivers) {
        $mysql_ok = 0;
        return "Only MySQL, Postgres and SQLite DBI drivers supported for Bio::DB::SeqFeature RDMS tests";
    }
    $mysql_ok = 1;
    return @drivers;
}

sub make_dbi_test {
    my $dsn = $build->notes('test_dsn') || return;
    my $path0 = File::Spec->catfile('t', 'LocalDB', 'SeqFeature.t');
    my $driver = $build->notes('dbd_driver');
    my $path = File::Spec->catfile('t', 'LocalDB',
                ($driver eq 'mysql')   ? 'SeqFeature_mysql.t'  :
                ($driver eq 'SQLite')  ? 'SeqFeature_SQLite.t'  :
                'SeqFeature_Pg.t');
    my $test_db = $build->notes('test_db');
    my $user = $build->notes('test_user');
    my $pass = $build->notes('test_pass');
    open my $F, '>', $path or die "Could not write file '$path' for DBI test: $!\n";
    my $str = "$path0 -adaptor DBI::$driver -create 1 -temp 1 -dsn \"$dsn\"";
    $str .= " -user $user" if $user;
    $str .= " -password $pass" if $pass;
    print $F <<END;
system '$^X $str';
END
    close $F;
    $build->add_to_cleanup($path);
    $build->add_to_cleanup($test_db) if $driver eq 'SQLite';
    #$build->add_to_manifest_skip($path);
}

sub test_biodbgff {
    eval {require DBI;};  # if not installed, this sub won't actually be called
    return if $@;
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
        foreach my $poss ('SQLite', 'mysql', 'Pg', 'Oracle') {
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
        elsif ($driver =~ /^[sS]/) {
            $driver = 'SQLite';
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
        if ($driver eq 'Pg' || $driver eq 'SQLite') {
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
                         "database tests (requires MySQL or Pg driver)\n") unless ($driver eq 'mysql' or $driver eq 'Pg');
    }
    else {
        $build->log_info("  - will not run the BioDBGFF or BioDBSeqFeature live database tests\n");
    }

    $build->log_info("\n");
    return $proceed;
}
