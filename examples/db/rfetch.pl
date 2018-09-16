#!/usr/bin/perl

#
# Fetch sequence data via OBDA registry system
#
# usage: rfetch -i <file_with_accession_list> -a -v -d embl -s start -e end
#

use Bio::DB::Registry;
use Bio::SeqIO;
use Getopt::Long;
use strict;

my $database = 'embl_biosql';
my $start    = undef;
my $end      = undef;
my $format   = 'fasta';
my $file     = undef;
my $acc      = undef;
my $verbose  = undef;

&GetOptions(
	    'd|database:s' => \$database,
	    's|start:i' => \$start,
	    'e|end:i'   => \$end,
	    'f|format:s' => \$format,
	    'i|input:s' => \$file,
	    'a|acc'     => \$acc,
	    'v|verbose' => \$verbose,
	   );


my $registry = Bio::DB::Registry->new();

my $db = $registry->get_database($database);

my $seqout = Bio::SeqIO->new( '-format' => $format, '-fh' => \*STDOUT);

my @ids;

if( defined $file ) {
  open my $F, '<', $file or die "Could not read file '$file': $!\n";
  while( <$F> ) {
    my ($id) = split;
    push(@ids,$id);
  }
  close $F;
} else {
  @ids = @ARGV;
}

foreach my $id ( @ids ) {
  my $seq;
  if( $verbose ){
    print STDERR "fetching $id\n";
  }

  if( $acc ) {
    $seq = $db->get_Seq_by_acc($id);
  } else {
    $seq = $db->get_Seq_by_id($id);
  }

  if( defined $start && defined $end ) {
    $seq = $seq->trunc($start,$end);
  }

  $seqout->write_seq($seq);
}
