package Bio::DB::GFF::Aggregator::none;

use strict;
use Bio::DB::GFF::Aggregator;
use vars qw($VERSION @ISA);

@ISA = qw(Bio::DB::GFF::Aggregator);
$VERSION = '0.10';

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
