package Bio::DB::GFF::Aggregator::clone;

use strict;

use Bio::DB::GFF::Aggregator;
use vars qw($VERSION @ISA);

$VERSION = '0.10';
@ISA = qw(Bio::DB::GFF::Aggregator);

# we look for features of type Sequence and add them to a pseudotype transcript
sub aggregate {
  my $self = shift;
  my $features = shift;
  my $factory  = shift;

  my $matchsub    = $self->match_sub($factory) or return $features;

  my (%clones,@result);
  for my $feature (@$features) {
    next unless $matchsub->($feature);

    if ($feature->method eq 'Sequence' && $feature->source eq 'Genomic_canonical') {
      $clones{$feature->group}{canonical} = $feature;
    } elsif ($feature->method eq 'Clone_left_end') {
      $clones{$feature->group}{left} = $feature;
    } elsif ($feature->method eq 'Clone_right_end') {
      $clones{$feature->group}{right} = $feature;
    }
  } continue {
    push @result,$feature;  # in case someone else wants to look at the components
  }

  for my $clone (keys %clones) {
    my $canonical = $clones{$clone}{canonical} or next;

    # the genomic_canonical doesn't tell us where the clone starts and stops
    # so don't assume it
    my $duplicate = $canonical->clone;   # make a duplicate of the feature
    my ($start,$stop) = $duplicate->strand > 0 ? ('start','stop') : ('stop','start');
    @{$duplicate}{$start,$stop} =(undef,undef);

    $duplicate->{$start} = $clones{$clone}{left}{$start}  if exists $clones{$clone}{left};
    $duplicate->{$stop}  = $clones{$clone}{right}{$stop}  if exists $clones{$clone}{right};
    $duplicate->method($self->method);
    push @result,$duplicate;
  }

  \@result;
}

sub method { 'clone' }
sub get_part_names {
  my $self = shift;
  return @{$self->{parts}} if exists $self->{parts};
  return qw(Clone_left_end Clone_right_end Sequence:genomic_canonical);
}

1;
