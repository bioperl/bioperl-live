#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use Getopt::Long;
use Bio::DB::GFF;

=head1 NAME

bp_meta_gff.pl - Get/set Bio::DB::GFF meta-data

=head1 SYNOPSIS

  # set the following meta data values
  % bp_meta_gff.pl -d testdb tag1=value1 tag2=value2

  # get the indicated meta data value
  % bp_meta_gff.pl -d testdb tag1 tag2

=head1 DESCRIPTION

This script gets or sets metadata in a Bio::DB::GFF database.  Not all
adaptors support this operation!  To set a series of tags, pass a set
of tag=value pairs to the script.  To get the contents of a series of
tags, pass the bare tag names.

The output from the get operation will be an easily parseable set of
tag=value pairs, one per line.

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options.
e.g. -d instead of --database.

   --database <dsn>      Mysql database name (default dbi:mysql:test)
   --adaptor <adaptor>   Mysql adaptor (default dbi::mysqlopt)
   --user    <user>      Username for mysql authentication
   --pass    <password>  Password for mysql authentication

=head1 SEE ALSO

L<Bio::DB::GFF>

=head1 AUTHOR

Lincoln Stein, lstein@cshl.org

Copyright (c) 2002 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

my ($DSN,$ADAPTOR,$USER,$PASSWORD);

GetOptions ('database:s'    => \$DSN,
	    'adaptor:s'     => \$ADAPTOR,
	    'user:s'      => \$USER,
	    'password:s'  => \$PASSWORD,
	   ) or (system('pod2text', $0), exit -1);

$DSN     ||= 'dbi:mysql:test';
$ADAPTOR ||= 'dbi::mysqlopt';

my @args;
push @args,(-user=>$USER)     if defined $USER;
push @args,(-pass=>$PASSWORD) if defined $PASSWORD;

my $db = Bio::DB::GFF->new(-adaptor=>$ADAPTOR,-dsn => $DSN,@args)
  or die "Can't open database: ",Bio::DB::GFF->error,"\n";

for my $pair (@ARGV) {
  my ($tag,$value) = split /=/,$pair;
  if ($value) {  # set operation
    $db->meta($tag,$value);
    unless ($db->meta($tag) eq $value) {
      print STDERR "value for '$tag' not set; perhaps this adaptor does not support meta data?\n";
    }
  } else {
    print "$tag=",$db->meta($tag),"\n";
  }
}

__END__
