=head1 NAME

Bio::DB::GFF::Adaptor::dbi::iterator - iterator for Bio::DB::GFF::Adaptor::dbi

=head1 SYNOPSIS

For internal use only

=head1 DESCRIPTION

This is an internal module that is used by the Bio::DB::GFF DBI
adaptor to return an iterator across a sequence feature query.  The
object has a single method, next_feature(), that returns the next
feature from the query.  The method next_seq() is an alias for
next_feature().

=head1 BUGS

None known yet.

=head1 SEE ALSO

L<Bio::DB::GFF>,

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

package Bio::DB::GFF::Adaptor::dbi::iterator;
use strict;
use Bio::Root::Version;

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
