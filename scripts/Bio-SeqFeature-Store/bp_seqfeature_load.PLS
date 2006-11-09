#!/usr/bin/perl
# $Id: bp_seqfeature_load.PLS,v 1.8 2006-11-09 06:57:12 jason Exp $

use strict;

use Getopt::Long;
use File::Spec;
use Bio::DB::SeqFeature::Store::GFF3Loader;
use Bio::DB::SeqFeature::Store;

my $DSN         = 'dbi:mysql:test';
my $SFCLASS     = 'Bio::DB::SeqFeature';
my $ADAPTOR     = 'DBI::mysql';
my $VERBOSE  = 1;
my $FAST     = 0;
my $TMP      = File::Spec->tmpdir();
my $CREATE   = 0;
my $USER     = '';
my $PASS     = '';

GetOptions(
	   'dsn=s'       => \$DSN,
	   'seqfeature=s'  => \$SFCLASS,
	   'adaptor=s'   => \$ADAPTOR,
	   'verbose!'    => \$VERBOSE,
	   'fast'       => \$FAST,
	   'T|temporary-directory' => \$TMP,
	   'create'      => \$CREATE,
	   'user=s'      => \$USER,
	   'password=s'  => \$PASS
	   ) || die <<END;
Usage: $0 [options] gff_file1 gff_file2...
  Options:
          -d --dsn        The database name ($DSN)
          -s --seqfeature The type of SeqFeature to create ($SFCLASS)
          -a --adaptor    The storage adaptor to use ($ADAPTOR)
          -v --verbose    Turn on verbose progress reporting
             --noverbose  Turn off verbose progress reporting
          -f --fast       Activate fast loading (only some adaptors)
          -T --temporary-directory  Specify temporary directory for fast loading ($TMP)
          -c --create     Create the database and reinitialize it (will erase contents)
          -u --user       User to connect to database as
          -p --password   Password to use to connect to database
END

if ($FAST) {
  -d $TMP && -w $TMP
    or die "Fast loading is requested, but I cannot write into the directory $TMP";
}

my @options;
@options = ($USER,$PASS) if $USER || $PASS;

my $store = Bio::DB::SeqFeature::Store->new(
					    -dsn     => $DSN,
					    -adaptor => $ADAPTOR,
					    -tmpdir  => $TMP,
					    -user    => $USER,
					    -pass    => $PASS,
					    -write    => 1,
					    -create   => $CREATE)
  or die "Couldn't create connection to the database";

$store->init_database('erase') if $CREATE;

my $loader = Bio::DB::SeqFeature::Store::GFF3Loader->new(-store    => $store,
							 -sf_class => $SFCLASS,
							 -verbose  => $VERBOSE,
							 -tmpdir   => $TMP,
							 -fast     => $FAST)
  or die "Couldn't create GFF3 loader";

# on signals, give objects a chance to call their DESTROY methods
$SIG{TERM} = $SIG{INT} = sub {  undef $loader; undef $store; die "Aborted..."; };

$loader->load(@ARGV);

exit 0;

