package Bio::DB::Flat::BDB::embl;

use strict;
use Bio::DB::Flat::BDB;
use vars '@ISA';

@ISA = qw(Bio::DB::Flat::BDB);

sub parse_one_record {
  my $self  = shift;
  my $fh    = shift;
  my $parser =
    $self->{embl_cached_parsers}{fileno($fh)} ||= Bio::SeqIO->new(-fh=>$fh,-format=>$self->default_file_format);
  my $seq = $parser->next_seq;
  my $ids = $self->seq_to_ids($seq);
  return $ids;
}

sub seq_to_ids {
  my $self = shift;
  my $seq  = shift;

  my $display_id = $seq->display_id;
  my $accession  = $seq->accession_number;
  my %ids;
  $ids{ID}       = $display_id;
  $ids{ACC}      = $accession   if defined $accession;
  return \%ids;
}

sub default_primary_namespace {
  return "ID";
}

sub default_secondary_namespaces {
  return qw(ACC);
}

sub default_file_format { "embl" }


1;
