#!/usr/local/bin/perl



use Bio::DB::Registry;
use Getopt::Long;
use strict;

my $database = 'embl';
my $start    = undef;
my $end      = undef;
my $format   = 'fasta';


&GetOptions(
	   'd|database:s' => \$database,
	   's|start:i' => \$start,
	   'e|end:i'   => \$end,
	   'f|format:s' => \$format
	   );


my $registry = Bio::DB::Registry->new();

my $db = $registry->get_database($database);

my $seqout = Bio::SeqIO->new( '-format' => $format, '-fh' => \*STDOUT);

foreach my $id ( @ARGV ) {
  my $seq = $db->get_by_Seq_id($id);

  if( defined $start && defined $end ) {
    $seq = $seq->trunc($start,$end);
  }

  $seqout->write_seq($seq);
}
