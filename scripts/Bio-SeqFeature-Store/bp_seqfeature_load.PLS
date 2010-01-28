#!/usr/bin/perl
# $Id$

use strict;

use Getopt::Long;
use File::Spec;
use Bio::DB::SeqFeature::Store::GFF3Loader;
use Bio::DB::SeqFeature::Store;

my $DSN         = 'dbi:mysql:test';
my $SFCLASS     = 'Bio::DB::SeqFeature';
my $ADAPTOR     = 'DBI::mysql';
my $NAMESPACE;
my $VERBOSE  = 1;
my $FAST     = 0;
my $TMP      = File::Spec->tmpdir();
my $IGNORE_SEQREGION   = 0;
my $CREATE   = 0;
my $USER     = '';
my $PASS     = '';
my $COMPRESS = 0;
my $INDEX_SUB = 1;

GetOptions(
       'dsn=s'          => \$DSN,
       'seqfeature=s'   => \$SFCLASS,
       'namespace=s'    => \$NAMESPACE,
       'adaptor=s'      => \$ADAPTOR,
       'verbose!'       => \$VERBOSE,
       'fast'           => \$FAST,
       'ignore-seqregion'           => \$IGNORE_SEQREGION,
       'T|temporary-directory=s'    => \$TMP,
       'create'         => \$CREATE,
       'user=s'         => \$USER,
       'password=s'     => \$PASS,
       'zip'            => \$COMPRESS,
       'S|subfeatures!' => \$INDEX_SUB,
       ) || die <<END;
Usage: $0 [options] gff_file1 gff_file2...
  Options:
          -d --dsn        The database name ($DSN)
          -n --namespace  The table prefix, excluding ($NAMESPACE)
          -s --seqfeature The type of SeqFeature to create ($SFCLASS)
          -a --adaptor    The storage adaptor to use ($ADAPTOR)
          -v --verbose    Turn on verbose progress reporting
             --noverbose  Turn off verbose progress reporting
          -f --fast       Activate fast loading (only some adaptors)
          -T --temporary-directory  Specify temporary directory for fast loading ($TMP)
          -c --create     Create the database and reinitialize it (will erase contents)
          -u --user       User to connect to database as
          -p --password   Password to use to connect to database
          -S --subfeatures   Turn on indexing of subfeatures (default)
             --nosubfeatures Turn off indexing of subfeatures
          -i --ignore-seqregion
                          If true, then ignore ##sequence-region directives in the
                          GFF3 file (default, create a feature for each region)
          -z --zip        If true, database tables will be compressed to save space

Please see http://www.sequenceontology.org/gff3.shtml for information
about the GFF3 format. BioPerl extends the format slightly by adding 
a ##index-subfeatures directive. Set this to a true value if you wish 
the database to be able to retrieve a feature's individual parts (such as the
exons of a transcript) independently of the top level feature:

  ##index-subfeatures 1

It is also possible to control the indexing of subfeatures on a case-by-case
basis by adding "index=1" or "index=0" to the feature's attribute list. This
should only be used for subfeatures.

Subfeature indexing is true by default. Set to false (0) to save lots
of database space and speed performance. You may use --nosubfeatures
to force this.  
END

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

my $loader = Bio::DB::SeqFeature::Store::GFF3Loader->new
(
    -store              => $store,
    -sf_class           => $SFCLASS,
    -verbose            => $VERBOSE,
    -tmpdir             => $TMP,
    -fast               => $FAST,
    -ignore_seqregion   => $IGNORE_SEQREGION,
    -index_subfeatures  => $INDEX_SUB,
)
or die "Couldn't create GFF3 loader";

# on signals, give objects a chance to call their DESTROY methods
$SIG{TERM} = $SIG{INT} = sub {  undef $loader; undef $store; die "Aborted..."; };

$loader->load(@ARGV);

exit 0;

