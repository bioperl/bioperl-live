package Bio::Graphics::Glyph::dna;

use strict;
use Bio::Graphics::Glyph::generic;
use vars '@ISA';
@ISA = qw(Bio::Graphics::Glyph::generic);

# turn off description
sub description { 0 }

# turn off label
# sub label { 1 }

sub height {
  my $self = shift;
  my $font = $self->font;
  return $self->dna_fits ? 2*$font->height
       : $self->do_gc    ? $self->SUPER::height
       : 0;
}

sub pixels_per_base {
  my $self = shift;

  my $width           = $self->width;
  my $length          = $self->feature->length;

  return $width/($self->feature->length-1);
}

sub dna_fits {
  my $self = shift;

  my $pixels_per_base = $self->pixels_per_base;
  my $font            = $self->font;
  my $font_width      = $font->width;

  return $pixels_per_base >= $font_width;
}

sub do_gc {
  my $self = shift;
  my $do_gc = $self->option('do_gc');
  return  if defined($do_gc) && !$do_gc;
  return  1;
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->bounds(@_);

  my $dna        = eval { $self->feature->seq };
  $dna or return;

  if ($self->dna_fits) {
    $self->draw_dna($gd,$dna,$x1,$y1,$x2,$y2);
  } elsif ($self->do_gc) {
    $self->draw_gc_content($gd,$dna,$x1,$y1,$x2,$y2);
  }
}

sub draw_dna {
  my $self = shift;

  my ($gd,$dna,$x1,$y1,$x2,$y2) = @_;
  my $pixels_per_base = $self->pixels_per_base;

  my @bases = split '',$dna;
  my $color = $self->fgcolor;
  my $font  = $self->font;
  my $lineheight = $font->height;

  my %complement = (g=>'c',a=>'t',t=>'a',c=>'g',
		    G=>'C',A=>'T',T=>'A',C=>'G');
  for (my $i=0;$i<@bases;$i++) {
    my $x = $x1 + $i * $pixels_per_base;
    $gd->char($font,$x,$y1,$bases[$i],$color);
    $gd->char($font,$x,$y1+$lineheight,$complement{$bases[$i]}||$bases[$i],$color);
  }

}

sub draw_gc_content {
  my $self     = shift;
  my $gd       = shift;
  my $dna = shift;
  my ($x1,$y1,$x2,$y2) = @_;

  my $bin_size = length($dna) / ($self->option('gc_bins') || 100);
  $bin_size = 100 if $bin_size < 100;

  my @bins;
  for (my $i = 0; $i < length($dna) - $bin_size; $i+= $bin_size) {
    my $subseq  = substr($dna,$i,$bin_size);
    my $gc      = $subseq =~ tr/gcGC/gcGC/;
    my $content = $gc/$bin_size;
    push @bins,$content;
  }
  push @bins,0.5 unless @bins;  # avoid div by zero
  my $bin_width  = ($x2-$x1)/@bins;
  my $bin_height = $y2-$y1;
  my $fgcolor    = $self->fgcolor;
  my $bgcolor    = $self->factory->translate_color($self->panel->gridcolor);
  my $axiscolor  = $self->color('axis_color') || $fgcolor;

  $gd->line($x1,  $y1,        $x1,  $y2,        $axiscolor);
  $gd->line($x2-1,$y1,        $x2-1,$y2,        $axiscolor);
  $gd->line($x1,  $y1,        $x1+3,$y1,        $axiscolor);
  $gd->line($x1,  $y2,        $x1+3,$y2,        $axiscolor);
  $gd->line($x1,  ($y2+$y1)/2,$x1+3,($y2+$y1)/2,$axiscolor);
  $gd->line($x2-3,$y1,        $x2-1, $y1,       $axiscolor);
  $gd->line($x2-3,$y2,        $x2-1, $y2,       $axiscolor);
  $gd->line($x2-3,($y2+$y1)/2,$x2-1,($y2+$y1)/2,$axiscolor);
  $gd->line($x1+5,$y2,        $x2-5,$y2,        $bgcolor);
  $gd->line($x1+5,($y2+$y1)/2,$x2-5,($y2+$y1)/2,$bgcolor);
  $gd->line($x1+5,$y1,        $x2-5,$y1,        $bgcolor);
  $gd->string($self->font,$x1+5,$y1,'% gc',$axiscolor);

  for (my $i = 0; $i < @bins; $i++) {
    my $bin_start  = $x1+$i*$bin_width;
    my $bin_stop   = $bin_start + $bin_width;
    my $y          = $y2 - ($bin_height*$bins[$i]);
    $gd->line($bin_start,$y,$bin_stop,$y,$fgcolor);
    $gd->line($bin_stop,$y,$bin_stop,$y2 - ($bin_height*$bins[$i+1]),$fgcolor)
      if $i < @bins-1;
  }
}

sub make_key_feature {
  my $self = shift;
  my $offset = $self->panel->offset;
  my $feature =
    Bio::Graphics::Feature->new(-start=> $offset,
				-stop => $offset + 10,
				-name => 'DNA/GC content',
				-strand => '+1',
			       );
  $feature;
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::dna - The "dna" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws DNA sequences.  At high magnifications, this glyph
will draw the actual base pairs of the sequence (both strands).  At
low magnifications, the glyph will plot the GC content.

For this glyph to work, the feature must return a DNA sequence string
in response to the dna() method.

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

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description               Default
  ------      -----------               -------

  -do_gc      Whether to draw the GC      true
              graph at low mags

  -gc_bins    Fixed number of intervals   100
              to sample across the
              panel.

  -axis_color Color of the vertical axes  fgcolor
              in the GC content graph

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

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
