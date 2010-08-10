=head1 NAME

Bio::DB::GFF::Adaptor::berkeleydb::iterator - iterator for Bio::DB::GFF::Adaptor::berkeleydb

=head1 SYNOPSIS

For internal use only

=head1 DESCRIPTION

This is an internal module that is used by the Bio::DB::GFF in-memory
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

package Bio::DB::GFF::Adaptor::berkeleydb::iterator;
use strict;
use DB_File qw(R_FIRST R_NEXT);

# this module needs to be cleaned up and documented
use Bio::Root::Version;
*next_seq = \&next_feature;

sub new {
  my $class = shift;
  my ($data,$callback,$tmpfile) = @_;
  return bless {data     => $data,
		callback => $callback,
		tmpfile  => $tmpfile,
                cache    => []},$class;
}

sub next_feature {
  my $self = shift;
  return shift @{$self->{cache}} if @{$self->{cache}};

  my $data     = $self->{data} or return;
  my $callback = $self->{callback};

  my $features;
  my $db = tied(%$data);
  my ($key,$value);

  for (my $status = $db->seq($key,$value,$self->{iter}++ ? R_NEXT : R_FIRST);
       $status == 0;
       $status = $db->seq($key,$value,R_NEXT)) {
    my @feature       = split ($;,$value);
    $features   = $callback->(@feature);
    last if $features;
  }

  unless ($features) {
    $features = $callback->();
    undef $self->{data};
    undef $self->{cache};
    unlink $self->{tmpfile};
  }

  $self->{cache} = $features or return;
  shift @{$self->{cache}};
}

1;
