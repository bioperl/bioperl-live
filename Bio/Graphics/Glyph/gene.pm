package Bio::Graphics::Glyph::gene;

# $Id$

use strict;
use base 'Bio::Graphics::Glyph::processed_transcript';

sub extra_arrow_length {
  my $self = shift;
  return 0 unless $self->{level} == 1;
  local $self->{level} = 0;  # fake out superclass
  return $self->SUPER::extra_arrow_length;
}

sub pad_right {
  my $self = shift;
  my $strand = $self->feature->strand;
  $strand *= -1 if $self->{flip};
  my $pad    = $self->SUPER::pad_right;
  return $pad unless defined($strand) && $strand > 0;
  my $al = $self->arrow_length;
  return $al > $pad ? $al : $pad;
}

sub bump {
  my $self = shift;
  return 1 if $self->{level} == 0; # top level bumps, other levels don't unless specified in config
  return $self->SUPER::bump;
}

sub label {
  my $self = shift;
  if ($self->label_transcripts && $self->feature->primary_tag eq 'mRNA') { # the mRNA
    return $self->_label;
  } else {
    return $self->SUPER::label;
  }
}

sub label_position {
  my $self = shift;
  return 'top' if $self->{level} == 0;
  return 'left';
}

sub label_transcripts {
  my $self = shift;
  return $self->option('label transcripts');
}

sub draw_connectors {
  my $self = shift;
  return if $self->feature->primary_tag eq 'gene';
  $self->SUPER::draw_connectors(@_);
}

sub maxdepth {
  my $self = shift;
  my $md   = $self->Bio::Graphics::Glyph::maxdepth;
  return $md if defined $md;
  return 2;
}


sub _subseq {
  my $class   = shift;
  my $feature = shift;
  return $feature->get_SeqFeatures('mRNA') if $feature->primary_tag eq 'gene';
  return $feature->get_SeqFeatures();# 'CDS',"5'_UTR","3'_UTR");
}

1;

__END__
