package Bio::DB::GFF::Aggregator::alignment;

use strict;

use Bio::DB::GFF::Aggregator;
use vars qw($VERSION @ISA);

@ISA = qw(Bio::DB::GFF::Aggregator);
$VERSION = '0.10';

# we look for features of type Sequence and add them to a pseudotype transcript
sub aggregate {
  my $self = shift;
  my $features = shift;
  my $factory  = shift;

  my $matchsub = $self->match_sub($factory) or return;
  my $method = $self->get_method;

  my (%alignments,%targets);

  warn "running aligner aggregator" if $factory->debug;
  for my $feature (@$features) {

    next unless $matchsub->($feature);

    my $group  = $feature->{group};
    my $source = $feature->source;
    unless (exists $alignments{$group,$source}) {
      my $type = Bio::DB::GFF::Typename->new($method,$source);

      my $f = $feature->clone;
      # this is a violation of OO encapsulation, but need to do it this way
      # to achieve desired performance
      @{$f}{qw(type score phase)} = ($type,undef,undef);

      $alignments{$group,$source} = $f or next;
    }

    my $main = $alignments{$group,$source};
    $main->add_subfeature($feature);
  }

  warn "running aligner adjuster" if $factory->debug;
  my @result;
  for my $alignment (values %alignments) {
    $alignment->adjust_bounds;
    push @result,$alignment;
  }
  warn "aligner done" if $factory->debug;
  \@result;
}

sub method { 'alignment' }
sub get_part_names {
  my $self = shift;
  return @{$self->{parts}} if exists $self->{parts};
  return qw(similarity);
}

1;
