#!/usr/bin/perl

use strict;
use lib '../blib/lib';
use Bio::DB::GFF;
use Getopt::Long;

=head1 NAME

load_gff.pl - Load a Bio::DB::GFF database from GFF files.

=head1 SYNOPSIS

  % load_gff.pl -d testdb -f dna.fa features1.gff features2.gff ...

=head1 DESCRIPTION

This script loads a Bio::DB::GFF database with the features contained
in a list of GFF files.  You must use the exact variant of GFF
described in L<Bio::DB::GFF>.  Various command-line options allow you
to control which database to load, whether to load DNA (in fasta
format) as well, and whether to allow an existing database to be
overwritten.

This script uses the Bio::DB::GFF interface, and so works with all
database adaptors currently supported by that module (MySQL, Oracle,
PostgreSQL soon).  However, it is slow.  For faster loading, see the
MySQL-specific L<bulk_load_gff.pl> and <fast_load_gff.pl> scripts.

NOTE: If no arguments are provided, or the filename is given as "-"
then the input is taken from standard input. Compressed files (.gz,
.Z, .bz2) are automatically uncompressed.  Fasta files must end in .fa
optionally followed by a compression suffix in order to be recognized.

To load FASTA files without GFF data, run like this:

 load_gff.pl --database my_database --fasta my_fasta_file_path </dev/null

If you get a lot of errors on loading, it may be that the database
schema is out of sync.  Run with the --upgrade option to update the
schema in a non-destructive manner.

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options.
e.g. -d instead of --database.

   --dsn     <dsn>       Data source (default dbi:mysql:test)
   --adaptor <adaptor>   Schema adaptor (default dbi::mysqlopt)
   --user    <user>      Username for mysql authentication
   --pass    <password>  Password for mysql authentication
   --fasta   <path>      Fasta file or directory containing fasta files for the DNA
   --create              Force creation and initialization of database
   --upgrade             Upgrade existing database to current schema

=head1 SEE ALSO

L<Bio::DB::GFF>, L<bulk_load_gff.pl>, L<load_gff.pl>

=head1 AUTHOR

Lincoln Stein <lstein@cshl.org>.

Copyright (c) 2002 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut


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

