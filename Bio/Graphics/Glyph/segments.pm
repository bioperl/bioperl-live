package Bio::Graphics::Glyph::segments;
#$Id$

use strict;
use Bio::Location::Simple;
use Bio::Graphics::Glyph::generic;
use Bio::Graphics::Glyph::segmented_keyglyph;
use vars '@ISA','$VERSION';

use constant RAGGED_START_FUZZ => 25;  # will show ragged ends of alignments
                                       # up to this many bp.

@ISA = qw( Bio::Graphics::Glyph::segmented_keyglyph
	   Bio::Graphics::Glyph::generic
	 );
$VERSION = '1.06';
my %complement = (g=>'c',a=>'t',t=>'a',c=>'g',n=>'n',
		  G=>'C',A=>'T',T=>'A',C=>'G',N=>'N');

sub pad_left {
  my $self = shift;
  return $self->SUPER::pad_left unless $self->option('draw_target') && $self->option('ragged_start') && $self->dna_fits;
  return $self->SUPER::pad_left unless $self->level > 0;
  my $target = $self->feature->hit;
  return $self->SUPER::pad_left unless $target->start<$target->end && $target->start < RAGGED_START_FUZZ;
  return ($target->start-1) * $self->scale;
}

sub pad_right {
  my $self = shift;
  return $self->SUPER::pad_right unless $self->level > 0;
  return $self->SUPER::pad_right unless $self->option('draw_target') && $self->option('ragged_start') && $self->dna_fits;
  my $target = $self->feature->hit;
  return $self->SUPER::pad_right unless $target->end < $target->start && $target->start < RAGGED_START_FUZZ;
  return ($target->end-1) * $self->scale;
}

# group sets connector to 'solid'
sub connector {
  my $self = shift;
  return $self->SUPER::connector(@_) if $self->all_callbacks;
  return ($self->SUPER::connector(@_) || 'solid');
}

# never allow our components to bump
sub bump {
  my $self = shift;
  return $self->SUPER::bump(@_) if $self->all_callbacks;
  return 0;
}

sub fontcolor {
  my $self = shift;
  return $self->SUPER::fontcolor unless $self->option('draw_target') || $self->option('draw_dna');
  return $self->SUPER::fontcolor unless $self->dna_fits;
  return $self->bgcolor;
}

sub draw_component {
  my $self = shift;
  my ($draw_dna,$draw_target) = ($self->option('draw_dna'),$self->option('draw_target'));
  return $self->SUPER::draw_component(@_)
    unless $draw_dna || $draw_target;
  return $self->SUPER::draw_component(@_) unless $self->dna_fits;

  my $dna = $draw_target ? eval {$self->feature->hit->seq}
                         : eval {$self->feature->seq};
  return $self->SUPER::draw_component(@_) unless length $dna > 0;  # safety

  my $show_mismatch = $draw_target && $self->option('show_mismatch');
  my $genomic = eval {$self->feature->seq} if $show_mismatch;

  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->bounds(@_);

  # adjust for nonaligned left end (for ESTs...)  The size given here is roughly sufficient
  # to show a polyA end or a c. elegans trans-spliced leader.
  my $offset = 0;
  eval {  # protect against data structures that don't implement the target() method.
    if ($draw_target && $self->option('ragged_start')){
      my $target = $self->feature->hit;
      if ($target->start < $target->end && $target->start < RAGGED_START_FUZZ  && $self->{partno} == 0) {
	$offset = $target->start - 1;
	if ($offset > 0) {
	  $dna       = $target->subseq(1-$offset,0)->seq . $dna;
	  $genomic   = $self->feature->subseq(1-$offset,0)->seq         . $genomic;
	  $x1        -= $offset * $self->scale;
	}
      }
      elsif ($target->end < $target->start && $target->end < RAGGED_START_FUZZ && $self->{partno} == $self->{total_parts}) {
	$offset = $target->end - 1;
	if ($offset > 0) {
	  $dna       .= $target->factory->get_dna($target,$offset,1);
	  $genomic   = $self->feature->subseq(-$offset,0)->seq . $genomic;
	  $x2        += $offset * $self->scale;
	  $offset = 0;
	}
      }
    }
  };

  $self->draw_dna($gd,$offset,$dna,$genomic,$x1,$y1,$x2,$y2);
}

sub draw_dna {
  my $self = shift;

  my ($gd,$start_offset,$dna,$genomic,$x1,$y1,$x2,$y2) = @_;
  my $pixels_per_base = $self->scale;
  my $complement      = $self->feature->strand < 0;

  my @bases   = split '',$dna;
  my @genomic = split '',$genomic;
  my $color = $self->fgcolor;
  my $font  = $self->font;
  my $lineheight = $font->height;
  my $fontwidth  = $font->width;
  $y1 -= $lineheight/2 - 3;
  my $pink = $self->factory->translate_color('lightpink');

  my $start  = $self->map_no_trunc($self->feature->start-$start_offset);
  my $offset = int (($x1-$start-1)/$pixels_per_base);

  for (my $i=$offset;$i<@bases;$i++) {
    my $x = int($start + $i * $pixels_per_base+0.5);
    next if $x+1 < $x1;
    last if $x+1 > $x2;
    if ($genomic[$i] && lc($bases[$i]) ne lc($complement ? $complement{$genomic[@genomic - $i - 1]} : $genomic[$i])) {
      $self->filled_box($gd,$x,$y1+3,$x+$fontwidth,$y1+$lineheight-3,$pink,$pink);
    }
    $gd->char($font,$x+1,$y1,$complement ? $complement{$bases[$i]} || $bases[$i] : $bases[$i],$color);
  }


}

# Override _subseq() method to make it appear that a top-level feature that
# has no subfeatures appears as a feature that has a single subfeature.
# Otherwise at high mags gaps will be drawn as components rather than
# as connectors.  Because of differing representations of split features
# in Bio::DB::GFF::Feature and Bio::SeqFeature::Generic, there is
# some breakage of encapsulation here.
sub _subseq {
  my $self    = shift;
  my $feature = shift;
  my @subseq  = $self->SUPER::_subseq($feature);
  return @subseq if @subseq;
  if ($self->level == 0 && !@subseq && !eval{$feature->compound}) {
    my($start,$end) = ($feature->start,$feature->end);
    ($start,$end) = ($end,$start) if $start > $end; # to keep Bio::Location::Simple from bitching
    #    return Bio::Location::Simple->new(-start=>$start,-end=>$end);
    return $self->feature;
  } else {
    return;
  }
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::segments - The "segments" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph is used for drawing features that consist of discontinuous
segments.  Unlike "graded_segments" or "alignment", the segments are a
uniform color and not dependent on the score of the segment.

=head2 OPTIONS

The following options are standard among all Glyphs.  See
L<Bio::Graphics::Glyph> for a full explanation.

  Option      Description                      Default
  ------      -----------                      -------

  -fgcolor      Foreground color	       black

  -outlinecolor	Synonym for -fgcolor

  -bgcolor      Background color               turquoise

  -fillcolor    Synonym for -bgcolor

  -linewidth    Line width                     1

  -height       Height of glyph		       10

  -font         Glyph font		       gdSmallFont

  -connector    Connector type                 0 (false)

  -connector_color
                Connector color                black

  -label        Whether to draw a label	       0 (false)

  -description  Whether to draw a description  0 (false)

  -strand_arrow Whether to indicate            0 (false)
                 strandedness

  -draw_dna     If true, draw the dna residues 0 (false)
                 when magnification level
                 allows.

  -draw_target  If true, draw the dna residues 0 (false)
                 of the TARGET sequence when
                 magnification level allows.
                 SEE NOTE.

  -ragged_start When combined with -draw_target, 0 (false)
                 draw a few bases beyond the end
                 of the alignment.  SEE NOTE.

  -show_mismatch When combined with -draw_target, 0 (false)
                 highlights mismatched bases in
                 pink.  SEE NOTE.

The -draw_target and -ragged_start options only work with seqfeatures
that implement the hit() method (Bio::SeqFeature::SimilarityPair).
The -ragged_start option is mostly useful for looking for polyAs and
cloning sites at the beginning of ESTs and cDNAs.  Currently there is
no way of activating ragged ends.  The length of the ragged starts is
hard-coded at 25 bp, and the color of mismatches is hard-coded as
light pink.

=head1 BUGS

Please report them.

=head1 SEE ALSO


L<Bio::Graphics::Panel>,
L<Bio::Graphics::Glyph>,
L<Bio::Graphics::Glyph::arrow>,
L<Bio::Graphics::Glyph::cds>,
L<Bio::Graphics::Glyph::crossbox>,
L<Bio::Graphics::Glyph::diamond>,
L<Bio::Graphics::Glyph::dna>,
L<Bio::Graphics::Glyph::dot>,
L<Bio::Graphics::Glyph::ellipse>,
L<Bio::Graphics::Glyph::extending_arrow>,
L<Bio::Graphics::Glyph::generic>,
L<Bio::Graphics::Glyph::graded_segments>,
L<Bio::Graphics::Glyph::heterogeneous_segments>,
L<Bio::Graphics::Glyph::line>,
L<Bio::Graphics::Glyph::pinsertion>,
L<Bio::Graphics::Glyph::primers>,
L<Bio::Graphics::Glyph::rndrect>,
L<Bio::Graphics::Glyph::segments>,
L<Bio::Graphics::Glyph::ruler_arrow>,
L<Bio::Graphics::Glyph::toomany>,
L<Bio::Graphics::Glyph::transcript>,
L<Bio::Graphics::Glyph::transcript2>,
L<Bio::Graphics::Glyph::translation>,
L<Bio::Graphics::Glyph::triangle>,
L<Bio::DB::GFF>,
L<Bio::SeqI>,
L<Bio::SeqFeatureI>,
L<Bio::Das>,
L<GD>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
