package Bio::Graphics::Glyph::generic;

use strict;
use Bio::Graphics::Glyph;
use vars '@ISA';
@ISA = 'Bio::Graphics::Glyph';

my %complement = (g=>'c',a=>'t',t=>'a',c=>'g',
		  G=>'C',A=>'T',T=>'A',C=>'G');

# new options are 'label'       -- short label to print over glyph
#                 'description'  -- long label to print under glyph
# label and description can be flags or coderefs.
# If a flag, label will be taken from seqname, if it exists or primary_tag().
#            description will be taken from source_tag().

sub pad_top {
  my $self = shift;
  my $top  = $self->option('pad_top');
  return $top if defined $top;
  my $pad = $self->SUPER::pad_top;
  $pad   += $self->labelheight if $self->label;
  $pad;
}
sub pad_bottom {
  my $self = shift;
  my $bottom  = $self->option('pad_bottom');
  return $bottom if defined $bottom;
  my $pad = $self->SUPER::pad_bottom;
  $pad   += $self->labelheight if $self->description;
  $pad;
}
sub pad_right {
  my $self = shift;
  my $pad = $self->SUPER::pad_right;
  my $label_width       = length($self->label||'') * $self->font->width;
  my $description_width = length($self->description||'') * $self->font->width;
  my $max = $label_width > $description_width ? $label_width : $description_width;
  my $right = $max - $self->width;
  return $pad > $right ? $pad : $right;
}

sub labelheight {
  my $self = shift;
  return $self->{labelheight} ||= $self->font->height;
}
sub label {
  my $self = shift;
  return if $self->{overbumped};  # set by the bumper when we have hit bump limit
  return unless $self->{level} == 0;
  return exists $self->{label} ? $self->{label}
                               : ($self->{label} = $self->_label);
}
sub description {
  my $self = shift;
  return if $self->{overbumped}; # set by the bumper when we have hit bump limit
  return unless $self->{level} == 0;
  return exists $self->{description} ? $self->{description}
                                     : ($self->{description} = $self->_description);
}
sub _label {
  my $self = shift;

  # allow caller to specify the label
  my $label = $self->option('label');
  return unless defined $label;
  return $label unless $label eq '1';
  return "1"    if $label eq '1 '; # 1 with a space


  # figure it out ourselves
  my $f = $self->feature;

  return $f->display_name if $f->can('display_name');
  return $f->info         if $f->can('info');   # deprecated API
  return $f->seq_id       if $f->can('seq_id');
  return eval{$f->primary_tag};
}
sub _description {
  my $self = shift;

  # allow caller to specify the long label
  my $label = $self->option('description');
  return unless defined $label;
  return $label unless $label eq '1';
  return "1"   if $label eq '1 ';

  return $self->{_description} if exists $self->{_description};
  return $self->{_description} = $self->get_description($self->feature);
}

sub get_description {
  my $self = shift;
  my $feature = shift;

  # common places where we can get descriptions
  return join '; ',$feature->notes if $feature->can('notes');
  return $feature->desc            if $feature->can('desc');

  my $tag = $feature->source_tag;
  return undef if $tag eq '';
  $tag;
}

sub draw {
  my $self = shift;
  $self->SUPER::draw(@_);
  $self->draw_label(@_)       if $self->option('label');
  $self->draw_description(@_) if $self->option('description');
}

sub draw_label {
  my $self = shift;
  my ($gd,$left,$top,$partno,$total_parts) = @_;
  my $label = $self->label or return;
  my $x = $self->left + $left;
  $x = $self->panel->left + 1 if $x <= $self->panel->left;
  my $font = $self->option('labelfont') || $self->font;
  $gd->string($font,
	      $x,
	      $self->top + $top,
	      $label,
	      $self->fontcolor);
}
sub draw_description {
  my $self = shift;
  my ($gd,$left,$top,$partno,$total_parts) = @_;
  my $label = $self->description or return;
  my $x = $self->left + $left;
  $x = $self->panel->left + 1 if $x <= $self->panel->left;
  $gd->string($self->font,
	      $x,
	      $self->bottom - $self->pad_bottom + $top,
	      $label,
	      $self->font2color);
}

sub dna_fits {
  my $self = shift;

  my $pixels_per_base = $self->scale;
  my $font            = $self->font;
  my $font_width      = $font->width;

  return $pixels_per_base >= $font_width;
}

sub arrowhead {
  my $self = shift;
  my $gd   = shift;
  my ($x,$y,$height,$orientation) = @_;

  my $fg = $self->set_pen;
  my $style = $self->option('arrowstyle') || 'regular';

  if ($style eq 'filled') {
    my $poly = new GD::Polygon;
    if ($orientation >= 0) {
      $poly->addPt($x-$height,$y-$height);
      $poly->addPt($x,$y);
      $poly->addPt($x-$height,$y+$height,$y);
    } else {
      $poly->addPt($x+$height,$y-$height);
      $poly->addPt($x,$y);
      $poly->addPt($x+$height,$y+$height,$y);
    }
    $gd->filledPolygon($poly,$fg);
  } else {
    if ($orientation >= 0) {
      $gd->line($x-$height,$y-$height,$x,$y,$fg);
      $gd->line($x,$y,$x-$height,$y+$height,$fg);
    } else {
      $gd->line($x+$height,$y-$height,$x,$y,$fg);
      $gd->line($x,$y,$x+$height,$y+$height,$fg);
    }
  }
}

sub arrow {
  my $self = shift;
  my $gd   = shift;
  my ($x1,$x2,$y) = @_;

  my $fg     = $self->set_pen;
  my $height = $self->height/3;

  $gd->line($x1,$y,$x2,$y,$fg);
  $self->arrowhead($gd,$x2,$y,$height,+1) if $x1 < $x2;
  $self->arrowhead($gd,$x2,$y,$height,-1) if $x2 < $x1;
}

sub reversec {
  $_[1]=~tr/gatcGATC/ctagCTAG/;
  return scalar reverse $_[1];
}

1;

=head1 NAME

Bio::Graphics::Glyph::generic - The "generic" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This is identical to the "box" glyph.  It is the default glyph used
when not otherwise specified.

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

  -pad_top      Top padding                    0

  -pad_bottom   Bottom padding                 0

  -label        Whether to draw a label	       0 (false)

  -description  Whether to draw a description  0 (false)

  -strand_arrow Whether to indicate            0 (false)
                 strandedness

-pad_top and -pad_bottom allow you to insert some blank space between
the glyph's boundary and its contents.  This is useful if you are
changing the glyph's height dynamically based on its feature's score.

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

Allen Day E<lt>day@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
