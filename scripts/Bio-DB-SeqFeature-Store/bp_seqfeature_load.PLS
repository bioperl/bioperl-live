#!/usr/bin/perl

use strict;
use warnings;

## Used to output the 'usage' message
use Pod::Usage;

## Used to parse command line options
use Getopt::Long;

## Used to create temporary files, if necessary
use File::Spec;

## BioPerl!
use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;



## The available options. Note, these defaults are 'hard coded' into
## the USAGE POD, so if you change one of the defaults (you shouldn't),
## you should update the USAGE.

my $DSN			= 'dbi:mysql:test';
my $SFCLASS		= 'Bio::DB::SeqFeature';
my $ADAPTOR		= 'DBI::mysql';
my $NAMESPACE;
my $VERBOSE		= 1;
my $FAST		= 0;
my $TMP			= File::Spec->tmpdir();
my $IGNORE_SEQREGION	= 0;
my $CREATE		= 0;
my $USER 		= '';
my $PASS 		= '';
my $COMPRESS		= 0;
my $INDEX_SUB		= 1;
my $NOALIAS_TARGET	= 0;
my $SUMMARY_STATS	= 0;
my $NOSUMMARY_STATS  = 0;

## Two flags based on http://stackoverflow.com/questions/1232116
## how-to-create-pod-and-use-pod2usage-in-perl
my $opt_help;
my $opt_man;

GetOptions( 'd|dsn=s'			=> \$DSN,
	    's|seqfeature=s'		=> \$SFCLASS,
	    'n|namespace=s'		=> \$NAMESPACE,
	    'a|adaptor=s'		=> \$ADAPTOR,
	    'v|verbose!'		=> \$VERBOSE,
	    'f|fast'			=> \$FAST,
	    'T|temporary-directory=s'	=> \$TMP,
	    'i|ignore-seqregion'	=> \$IGNORE_SEQREGION,
	    'c|create'			=> \$CREATE,
	    'u|user=s'			=> \$USER,
	    'p|password=s'		=> \$PASS,
	    'z|zip'			=> \$COMPRESS,
	    'S|subfeatures!'		=> \$INDEX_SUB,

	    ## Any good single letter choices here?
	    'noalias-target'		=> \$NOALIAS_TARGET,
	    'summary'			=> \$SUMMARY_STATS,
        'N|nosummary'    => \$NOSUMMARY_STATS,

	    ## I miss '--help' when it isn't there!
	    'h|help!'			=> \$opt_help,
	    'm|man!'			=> \$opt_man,
	  )
  or pod2usage( -message =>
		"\nTry 'bp_seqfeature_load.pl --help' for more information\n",
		-verbose => 0,
		-exitval => 2,
	      );

## Should we output usage information?
pod2usage( -verbose => 1 ) if $opt_help;
pod2usage( -verbose => 2 ) if $opt_man;

## Did we get any files to process?
@ARGV
  or pod2usage( -message =>
		"\nYou need to pass some GFF or fasta files to load\n",
		-verbose => 0,
		-exitval => 2,
	      );



## POD

=head1 NAME

bp_seqfeature_load.pl - Load GFF into a SeqFeature database

=head1 DESCRIPTION

Pass any number of GFF or fasta format files (or GFF with embedded
fasta) to load the features and sequences into a SeqFeature
database. The database (and adaptor) to use is specified on the
command line. Use the --create flag to create a new SeqFeature
database.

=head1 SYNOPSIS

 bp_seqfeature_load.pl [options] gff_or_fasta_file1 [gff_or_fasta_file2 [...]]

Try 'bp_seqfeature_load.pl --help' or '--man' for more information.

=head1 OPTIONS

=over 4

=item -d, --dsn

DBI data source (default dbi:mysql:test)

=item -n, --namespace

The table prefix to use (default undef) Allows several independent
sequence feature databases to be stored in a single database

=item -s, --seqfeature

The type of SeqFeature to create... RTSC (default Bio::DB::SeqFeature)

=item -a, --adaptor

The storage adaptor (class) to use (default DBI::mysql)

=item -v, --verbose

Turn on verbose progress reporting (default true) Use --noverbose to
switch this off.

=item -f, --fast

Activate fast loading. (default 0) Only available for some adaptors.

=item -T, --temporary-directory

Specify temporary directory for fast loading (default
File::Spec->tmpdir())

=item -i, --ignore-seqregion

If true, then ignore ##sequence-region directives in the GFF3 file
(default, create a feature for each region)

=item -c, --create

Create the database and reinitialize it (default false) Note, this
will erase previous database contents, if any.

=item -u, --user

User to connect to database as

=item -p, --password

Password to use to connect to database

=item -z, --zip

Compress database tables to save space (default false)

=item -S, --subfeatures

Turn on indexing of subfeatures (default true) Use --nosubfeatures to
switch this off.

=item --summary

Generate summary statistics for coverage graphs (default false) This
can be run on a previously loaded database or during the load. It will
default to true if --create is used.

=item -N, --nosummary

Do not generate summary statistics to save some space and load time (default if
--create is not specified, use this option to explicitly turn off summary
statistics when --create is specified)

=item --noalias-target

Don't create an Alias attribute whose value is the target_id in a
Target attribute (if the feature contains a Target attribute, the
default is to create an Alias attribute whose value is the target_id
in the Target attribute)

=back

Please see http://www.sequenceontology.org/gff3.shtml for information
about the GFF3 format. BioPerl extends the format slightly by adding a
##index-subfeatures directive. Set this to a true value if you wish
the database to be able to retrieve a feature's individual parts (such
as the exons of a transcript) independently of the top level feature:

  ##index-subfeatures 1

It is also possible to control the indexing of subfeatures on a
case-by-case basis by adding "index=1" or "index=0" to the feature's
attribute list. This should only be used for subfeatures.

Subfeature indexing is true by default. Set to false (0) to save lots
of database space and speed performance. You may use --nosubfeatures
to force this.

=cut





if ($FAST) {
  -d $TMP && -w $TMP
    or die "Fast loading is requested, but I cannot write into the directory $TMP";
  $DSN .= ";mysql_local_infile=1" if $ADAPTOR =~ /mysql/i && $DSN !~ /mysql_local_infile/;
}

my @options;
@options = ($USER,$PASS) if $USER || $PASS;

my $store = Bio::DB::SeqFeature::Store->new
(
    -dsn        => $DSN,
    -namespace  => $NAMESPACE,
    -adaptor    => $ADAPTOR,
    -tmpdir     => $TMP,
    -user       => $USER,
    -pass       => $PASS,
    -write      => 1,
    -create     => $CREATE,
    -compress   => $COMPRESS,
)
or die "Couldn't create connection to the database";

$store->init_database('erase') if $CREATE;
$SUMMARY_STATS++               if $CREATE; # this is a good thing

my $loader = Bio::DB::SeqFeature::Store::GFF3Loader->new
(
    -store              => $store,
    -sf_class           => $SFCLASS,
    -verbose            => $VERBOSE,
    -tmpdir             => $TMP,
    -fast               => $FAST,
    -ignore_seqregion   => $IGNORE_SEQREGION,
    -index_subfeatures  => $INDEX_SUB,
    -noalias_target     => $NOALIAS_TARGET,
    -summary_stats      => $NOSUMMARY_STATS ? 0 : $SUMMARY_STATS,
)
or die "Couldn't create GFF3 loader";

# on signals, give objects a chance to call their DESTROY methods
$SIG{TERM} = $SIG{INT} = sub {  undef $loader; undef $store; die "Aborted..."; };

$loader->load(@ARGV);

exit 0;
