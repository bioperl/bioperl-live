#!/usr/bin/perl

use strict;
use warnings;
use lib '../blib/lib';
use Bio::DB::GFF;
use Getopt::Long;

=head1 NAME

bp_load_gff.pl - Load a Bio::DB::GFF database from GFF files.

=head1 SYNOPSIS

  % bp_load_gff.pl -d testdb -u user -p pw
     --dsn 'dbi:mysql:database=dmel_r5_1;host=myhost;port=myport'
        dna1.fa dna2.fa features1.gff features2.gff ...

=head1 DESCRIPTION

This script loads a Bio::DB::GFF database with the features contained
in a list of GFF files and/or FASTA sequence files.  You must use the
exact variant of GFF described in L<Bio::DB::GFF>.  Various
command-line options allow you to control which database to load and
whether to allow an existing database to be overwritten.

This script uses the Bio::DB::GFF interface, and so works with all
database adaptors currently supported by that module (MySQL, Oracle,
PostgreSQL soon).  However, it is slow.  For faster loading, see the
MySQL-specific L<bp_bulk_load_gff.pl> and L<bp_fast_load_gff.pl> scripts.

=head2 NOTES

If the filename is given as "-" then the input is taken from standard
input. Compressed files (.gz, .Z, .bz2) are automatically
uncompressed.

FASTA format files are distinguished from GFF files by their filename
extensions.  Files ending in .fa, .fasta, .fast, .seq, .dna and their
uppercase variants are treated as FASTA files.  Everything else is
treated as a GFF file.  If you wish to load -fasta files from STDIN,
then use the -f command-line swith with an argument of '-', as in 

    gunzip my_data.fa.gz | bp_fast_load_gff.pl -d test -f -

On the first load of a database, you will see a number of "unknown
table" errors.  This is normal.

About maxfeature: the default value is 100,000,000 bases.  If you have
features that are close to or greater that 100Mb in length, then the
value of maxfeature should be increased to 1,000,000,000, or another
power of 10.

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options.
e.g. -d instead of --database.

   --dsn     <dsn>       Data source (default dbi:mysql:test)
   --adaptor <adaptor>   Schema adaptor (default dbi::mysqlopt)
   --user    <user>      Username for mysql authentication
   --pass    <password>  Password for mysql authentication
   --fasta   <path>      Fasta file or directory containing fasta files for the DNA
   --create              Force creation and initialization of database
   --maxfeature          Set the value of the maximum feature size (default 100 Mb; must be a power of 10)
   --group               A list of one or more tag names (comma or space separated)
                          to be used for grouping in the 9th column.
   --upgrade             Upgrade existing database to current schema
   --gff3_munge          Activate GFF3 name munging (see Bio::DB::GFF)
   --quiet               No progress reports
   --summary             Generate summary statistics for drawing coverage histograms.
                           This can be run on a previously loaded database or during
                           the load.

=head1 SEE ALSO

L<Bio::DB::GFF>, L<bulk_load_gff.pl>, L<bp_load_gff.pl>

=head1 AUTHOR

Lincoln Stein, lstein@cshl.org

Copyright (c) 2002 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

my ($DSN,$ADAPTOR,$CREATE,$USER,$PASSWORD,$FASTA,$UPGRADE,$MAX_BIN,$GROUP_TAG,$MUNGE,$QUIET,$SUMMARY_STATS);

GetOptions ('dsn:s'                  => \$DSN,
	    'adaptor:s'              => \$ADAPTOR,
	    'u|user:s'               => \$USER,
	    'p|password:s'           => \$PASSWORD,
            'fasta:s'                => \$FASTA,
            'upgrade'                => \$UPGRADE,
            'maxbin|maxfeature:s'    => \$MAX_BIN,
	    'group:s'                => \$GROUP_TAG,
	    'gff3_munge'             => \$MUNGE,
	    'quiet'                  => \$QUIET,
	    'summary'                => \$SUMMARY_STATS,
	    'create'                 => \$CREATE) or (system('pod2text',$0), exit -1);

# some local defaults
$DSN     ||= 'dbi:mysql:test';
$ADAPTOR ||= 'dbi::mysqlopt';
$MAX_BIN ||= 1_000_000_000; # to accomodate human-sized chromosomes

my @args;
push @args,(-user=>$USER)     if defined $USER;
push @args,(-pass=>$PASSWORD) if defined $PASSWORD;
push @args,(-preferred_groups=>[split(/[,\s+]+/,$GROUP_TAG)]) if defined $GROUP_TAG;
push @args,(-create=>1)       if $CREATE;
push @args,(-write=>1);

my $db = Bio::DB::GFF->new(-adaptor=>$ADAPTOR,-dsn => $DSN,@args)
  or die "Can't open database: ",Bio::DB::GFF->error,"\n";

$db->gff3_name_munging(1) if $MUNGE;

if ($CREATE) {
    $SUMMARY_STATS++;
    $MAX_BIN ? $db->initialize(-erase=>1,-MAX_BIN=>$MAX_BIN) :
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

my (@gff,@fasta);
foreach (@ARGV) {
  if (/\.(fa|fasta|dna|seq|fast)$/i) {
    push @fasta,$_;
  } else {
    push @gff,$_;
  }
}

for my $file (@gff) {
  warn "$file: loading...\n";
  my $loaded = $db->load_gff($file,!$QUIET);
  warn "$file: $loaded records loaded\n";
}

unshift @fasta,$FASTA if defined $FASTA;

for my $file (@fasta) {
  warn "Loading fasta ",(-d $file?"directory":"file"), " $file\n";
  my $loaded = $db->load_fasta($file,!$QUIET);
  warn "$file: $loaded records loaded\n";
}

if ($SUMMARY_STATS) {
    warn "Building summary statistics for coverage histograms...\n";
    $db->build_summary_statistics;
}
