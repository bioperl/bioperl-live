package Bio::DB::GFF::Adaptor::dbi::iterator;
use strict;

use constant STH         => 0;
use constant CALLBACK    => 1;
use constant CACHE       => 2;

*next_seq = \&next_feature;

sub new {
  my $class = shift;
  my ($sth,$callback) = @_;
  return bless [$sth,$callback,[]],$class;
}

sub next_feature {
  my $self = shift;
  return shift @{$self->[CACHE]} if @{$self->[CACHE]};
  my $sth = $self->[STH] or return;
  my $callback = $self->[CALLBACK];

  my $features;
  while (1) {
    if (my @row = $sth->fetchrow_array) {
      $features = $callback->(@row);
      last if $features;
    } else {
      $sth->finish;
      undef $self->[STH];
      $features = $callback->();
      last;
    }
  }
  $self->[CACHE] = $features or return;
  shift @{$self->[CACHE]};
}

1;
