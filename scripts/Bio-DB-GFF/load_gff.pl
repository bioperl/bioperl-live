#!/usr/bin/perl

use strict;
use lib '../blib/lib';
use Bio::DB::GFF;
use Getopt::Long;

my ($DSN,$ADAPTOR,$CREATE,$USER,$PASSWORD,$FASTA,$UPGRADE);

GetOptions ('dsn:s'       => \$DSN,
	    'adaptor:s'   => \$ADAPTOR,
	    'user:s'      => \$USER,
	    'password:s'  => \$PASSWORD,
            'fasta:s'     => \$FASTA,
            'upgrade'     => \$UPGRADE,
	    create        => \$CREATE) or die <<USAGE;
Usage: $0 [options] <gff file 1> <gff file 2> ...
Load a Bio::DB::GFF database from GFF files.

 Options:
   --dsn     <dsn>       Data source (default dbi:mysql:test)
   --adaptor <adaptor>   Schema adaptor (default dbi::mysqlopt)
   --user    <user>      Username for mysql authentication
   --pass    <password>  Password for mysql authentication
   --fasta   <path>      Fasta file or directory containing fasta files for the DNA
   --create              Force creation and initialization of database
   --upgrade             Upgrade existing database to current schema

If no GFF files are specified, or the file name is given as "-", then
the data is taken from standard input.  Compressed files (.gz, .Z, .bz2)
are uncompressed automatically.

To load FASTA files without GFF data, run like this:

 load_gff.pl --database my_database --fasta my_fasta_file_path </dev/null

If you get a lot of errors on loading, it may be that the database
schema is out of sync.  Run with the --upgrade option to update the
schema in a non-destructive manner.

USAGE
;

# some local defaults
$DSN     ||= 'dbi:mysql:test';
$ADAPTOR ||= 'dbi::mysqlopt';

my @auth;
push @auth,(-user=>$USER)     if defined $USER;
push @auth,(-pass=>$PASSWORD) if defined $PASSWORD;

my $db = Bio::DB::GFF->new(-adaptor=>$ADAPTOR,-dsn => $DSN,@auth)
  or die "Can't open database: ",Bio::DB::GFF->error,"\n";

if ($CREATE) {
  $db->initialize(1);
} elsif ($UPGRADE) {
  warn qq(expect to see several "table already exists" messages\n);
  $db->initialize(0);
  my $dbi = $db->dbh;  # get the raw database handle
  my ($count) = $dbi->selectrow_array('SELECT COUNT(*) FROM fnote');
  if (defined($count) && $count > 0) {
    warn qq(fnote table detected.  Translating into fattribute table.  This may take a while.\n);
    $dbi->do("INSERT INTO fattribute VALUES (1,'Note')") or die "failed: ",$dbi->errstr;
    $dbi->do("INSERT INTO fattribute_to_feature (fid,fattribute_id,fattribute_value) SELECT fnote.fid,1,fnote FROM fnote") or die "failed: ",$dbi->errstr;
    warn qq(Schema successfully upgraded.  You might want to drop the fnote table when you're sure everything's working.\n);
  }
}

for my $file (@ARGV) {
  warn "$file: loading...\n";
  my $loaded = $db->load_gff($file);
  warn "$file: $loaded records loaded\n";
}

if ($FASTA) {
  warn "Loading fasta ",(-d $FASTA?"directory":"file"), " $FASTA\n";
  my $loaded = $db->load_fasta($FASTA);
  warn "$FASTA: $loaded records loaded\n";
}

