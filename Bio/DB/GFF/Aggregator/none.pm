package Bio::DB::GFF::Aggregator::none;

use strict;
use carp 'croak';
use base 'Bio::DB::GFF::AggregatorI';
our $VERSION = '0.10';

sub disaggregate {
  my $self  = shift;
  my $types = shift;
  $types;          # no change
}

sub aggregate {
  my $self = shift;
  my $features = shift;
  $features;   # no change
}

1;
