package Bio::DB::Flat::BDB::fasta;

use strict;
use Bio::DB::Flat::BDB;
use vars '@ISA';

@ISA = qw(Bio::DB::Flat::BDB);

sub parse_one_record {
  my $self  = shift;
  my $fh    = shift;

  undef $self->{fasta_stored_id} if exists $self->{fasta_stored_fh}
    && $fh ne $self->{fasta_stored_fh} ;
  $self->{fasta_stored_fh} = $fh;

  while (<$fh>) {		# don't try this at home
    if (/^>(\S+)/) {
      my $id = $self->{fasta_stored_id};
      $self->{fasta_stored_id} = $1;
      next unless defined $id;
      return ($id,-length($_));
    }
  }
  # we get here at the end of the file
  return $self->{fasta_stored_id};
}

sub default_file_format { "fasta" }

1;
