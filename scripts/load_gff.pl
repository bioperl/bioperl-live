#!/usr/bin/perl

use strict;
use lib '../blib/lib';
use Bio::DB::GFF;
use Getopt::Long;

my ($DSN,$ADAPTOR,$CREATE,$USER,$PASSWORD);

GetOptions ('dsn:s'       => \$DSN,
	    'adaptor:s'   => \$ADAPTOR,
	    'user:s'      => \$USER,
	    'password:s'  => \$PASSWORD,
	    create    => \$CREATE) or die <<USAGE;
Usage: $0 [options] <gff file 1> <gff file 2> ...
Load a Bio::DB::GFF database from GFF files.

 Options:
   --dsn     <dsn>       Data source (default dbi:mysql:test)
   --adaptor <adaptor>   Schema adaptor (default dbi::mysqlopt)
   --user    <user>      Username for mysql authentication
   --pass    <password>  Password for mysql authentication
   --create              Force creation and initialization of database

If no GFF files are specified, or the file name is given as "-", then
the data is taken from standard input.  Compressed files (.gz, .Z, .bz2)
are uncompressed automatically.
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

$db->initialize(1) if $CREATE;
@ARGV = '-' unless @ARGV;

for my $file (@ARGV) {
  warn "$file: loading...\n";
  my $loaded = $db->load($file);
  warn "$file: $loaded records loaded\n";
}

