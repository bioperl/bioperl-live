#!/usr/bin/env perl
# AUTHOR: malcolm.cook@stowers-institute.org

use strict;
use warnings;

use Carp;
use Getopt::Long;
use File::Spec;
use Bio::DB::SeqFeature::Store;

#use Carp::Always;

my $DSN;
my $ADAPTOR;
my $VERBOSE  = 1;
my $USER     = '';
my $PASS     = '';
my @gff3opt;

GetOptions(
	   'dsn=s'       => \$DSN,
	   'adaptor=s'   => \$ADAPTOR,
	   'user=s'      => \$USER,
	   'password=s'  => \$PASS,
	   'gff3opt=i{,}'    => \@gff3opt,
	  ) || die <<END;
Usage: $0 [options] -- [WHICH FEATURES]
  output GFF3 for selected database features
  Options:
          -d --dsn        The database name ($DSN)
          -a --adaptor    The storage adaptor to use ($ADAPTOR)
          -u --user       User to connect to database as
          -p --password   Password to use to connect to database
          -g --gff3opt    flag options to gff3_string (i.e.: pass -gffopt 1 to recurse)

         WHICH FEATURES: any remaining options after '--' will select
         a subset of features. The arguments are identical to those
         accepted by Bio::DB::SeqFeature::Store->features().

END

$ADAPTOR     ||= 'DBI::mysql';
$DSN         ||= $ADAPTOR eq 'DBI::mysql' ? "mysql_read_default_file=$ENV{HOME}/.my.cnf" : '';

my $store = Bio::DB::SeqFeature::Store->new(
					    -dsn     => $DSN,
					    -adaptor => $ADAPTOR,
					    -user    => $USER,
					    -pass    => $PASS,
					   )
  or die "Couldn't create connection to the database";

# on signals, give objects a chance to call their DESTROY methods
$SIG{TERM} = $SIG{INT} = sub {   undef $store; die "Aborted..."; };

my $seq_stream = $store->get_seq_stream(@ARGV)  or die "failed to get_seq_stream(@ARGV)"; 
while (my $seq = $seq_stream->next_seq) {
  ### 20100725 // genehack
  # Try to call a gff3_string() method, but fall back to gff_string() if $seq
  # doesn't support that. Note that gff_string() is required per
  # Bio::SeqFeatureI, while gff3_string() is not. Currently, only
  # Bio::SeqFeature::Lite implements gff3_string().
  if ( $seq->can( 'gff3_string' )) {
    print $seq->gff3_string(@gff3opt) . "\n";
  }
  elsif ( $seq->can( 'gff_string' )) {
    # since we intend on getting a GFF3 string, make sure to pass the version
    $seq->gff_format->gff_version(3);
    print $seq->gff_string() . "\n";
  }
  else {
    confess "sequence object $seq does not support gff3_string() or gff_string() methods!"
  }
}

exit 0;
