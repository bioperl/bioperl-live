package Bio::DB::GFF::Aggregator;

use strict;
use Bio::DB::GFF::Util::Rearrange;  # for rearrange()
use Bio::DB::GFF::Feature;
use vars qw($VERSION @ISA);

$VERSION = '0.10';
@ISA = qw(Bio::Root::RootI);

sub new {
  my $class = shift;
  my ($args) = rearrange([],@_);
  $args = {} unless $args;
  return bless $args,$class;
}

# this is called at the beginning to turn the pseudo-type 
# into its component feature types
sub disaggregate {
  my $self  = shift;
  my $types = shift;
  my $factory = shift;

  my $sub_features = Bio::DB::GFF->parse_types($self->get_part_names);
  my $main_feature = Bio::DB::GFF->parse_types($self->get_main_name);

  my (@synthetic_types,@result);
  foreach (@$types) {
    my ($method,$source) = @$_;
    if (lc($method) eq $self->method) { # e.g. "transcript"
      push @synthetic_types,map { [$_->[0],$_->[1] || $source] } @$sub_features,@$main_feature;
    } else {
      push @result,$_;
    }
  }

  # remember what we're searching for
  $self->components(\@synthetic_types);
  push @result,@synthetic_types;
  \@result;
}

sub match_sub {
  my $self    = shift;
  my $factory = shift;

  my $types_to_aggregate = $self->components();  # saved from disaggregate call
  return unless @$types_to_aggregate;
  return $factory->make_match_sub($types_to_aggregate);
}

sub components {
  my $self = shift;
  my $d = $self->{components};
  $self->{components} = shift if @_;
  return unless ref $d;
  return wantarray ? @$d : $d;
}

sub aggregate {
  my $self = shift;
  my $features = shift;
  my $factory  = shift;

  $self->throw("aggregate() method must be overridden by a subclass");
}

sub method {
  shift->throw("method() method must be overriden by a subclass");
}

sub get_part_names {
  my $self = shift;
  if (exists $self->{PARTS}) {
    return ref $self->{PARTS} ? @{$self->{PARTS}} : $self->{PARTS};
  } else {
    return $self->part_names;
  }
}

sub get_main_name {
  my $self = shift;
  return $self->{MAIN} if exists $self->{MAIN};
  return $self->main_name;
}

sub part_names {
  return;
}
sub main_name {
  return;
}

1;
