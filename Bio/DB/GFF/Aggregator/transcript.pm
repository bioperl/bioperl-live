package Bio::DB::GFF::Aggregator::transcript;

use strict;
use Carp 'croak';
use base 'Bio::DB::GFF::Aggregator';
our $VERSION = '0.10';

# we look for features of type Sequence and add them to a pseudotype transcript
sub aggregate {
  my $self = shift;
  my $features = shift;
  my $factory  = shift;

  my $main_method = $self->main_name;
  my $matchsub    = $self->match_sub($factory) or return $features;

  my (@result,%transcripts);
  for my $feature (@$features) {
    next unless $matchsub->($feature);
    if ($feature->type =~ /$main_method/) {
      $transcripts{$feature->group}{base} ||= $feature;
    } else {
      push @{$transcripts{$feature->group}{subparts}},$feature;
    }
  } continue {
    push @result,$feature;   # in case someone else wants to look
  }

  # aggregate transcripts
  my $pseudo_method = $self->method;

  foreach (keys %transcripts) {
    next unless exists $transcripts{$_}{base};
    my $base = $transcripts{$_}{base};
    $base->method($pseudo_method);
    foreach (@{$transcripts{$_}{subparts}}) {
      $base->add_subfeature($_);
    }
    $base->adjust_bounds;
    push @result,$base;
  }

  \@result;
}

sub method { 'transcript' }

sub part_names {
  return qw(intron exon CDS);
}

sub main_name {
  return 'Sequence:curated';
}

1;
