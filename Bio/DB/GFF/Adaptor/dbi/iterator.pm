package Bio::DB::GFF::Adaptor::dbi::iterator;
use strict;

use constant STH         => 0;
use constant CALLBACK    => 1;

*next_seq = \&next_feature;

sub new {
  my $class = shift;
  my ($sth,$callback) = @_;
  return bless [$sth,$callback],$class;
}

sub next_feature {
  my $self = shift;
  my $sth = $self->[STH] or return;
  my $callback = $self->[CALLBACK];

  my $feature;
  while (1) {
    if (my @row = $sth->fetchrow_array) {
      $feature = $callback->(@row);
      last if $feature;
    } else {
      $sth->finish;
      undef $self->[STH];
      $feature = $callback->();
      last;
    }
  }
  return $feature;
}

1;
