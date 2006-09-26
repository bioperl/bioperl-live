package Bio::Graphics::Glyph::protein;

use strict;
use base qw(Bio::Graphics::Glyph::generic);

# turn off description
sub description { 0 }

# turn off label
# sub label { 1 }

sub height {
  my $self = shift;
  my $font = $self->font;
  return $self->dna_fits ? 2 * $font->height
       : $self->do_kd    ? $self->SUPER::height
       : 0;
}

sub do_kd {
  my $self = shift;
  my $do_kd = $self->option('do_kd');
  return  if defined($do_kd) && !$do_kd;
  return  1;
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->bounds(@_);

  my $protein = eval { $self->feature->seq };
  $protein = $protein->seq
      if ref($protein) and $protein->can('seq'); # to catch Bio::PrimarySeqI objects
  $protein or return;

  # workaround for my misreading of interface -- LS
  $protein = $protein->seq if ref($protein) && $protein->can('seq');

  if ($self->dna_fits) {
    $self->draw_protein($gd,$protein,$x1,$y1,$x2,$y2);
  } elsif ($self->do_kd) {
    $self->draw_kd_plot($gd,$protein,$x1,$y1,$x2,$y2);
  }
}

sub draw_protein {
  my $self = shift;

  my ($gd,$protein,$x1,$y1,$x2,$y2) = @_;
  my $pixels_per_base = $self->scale;
  
  my $feature = $self->feature;

  my @bases = split '', $protein;
  my $color = $self->fgcolor;
  my $font  = $self->font;
  my $lineheight = $font->height;
  $y1 -= $lineheight/2 - 3;

  my $start   = $self->panel->left + $self->map_pt($feature->start);
  my $end     = $self->panel->left + $self->map_pt($feature->end);

  my $offset  = int(($x1-$start-1)/$pixels_per_base);

  for (my $i=$offset;$i<@bases;$i++) {
    my $x = $start + $i * $pixels_per_base;
    next if $x+1 < $x1;
    last if $x > $x2;
    $gd->char($font,$x+2,$y1,$bases[$i],$color);
  }

}

sub draw_kd_plot {
  my $self     = shift;
  my $gd       = shift;
  my $protein = shift;
  my ($x1,$y1,$x2,$y2) = @_;

  my $kd_window = $self->option('kd_window') || 9;

  # Calculate the KD plot ...

  my %scores = ( I => 4.5,
		 V => 4.2,
		 L => 3.8,
		 F => 2.8,
		 C => 2.5,
		 M => 1.9,
		 A => 1.8,
		 G => -0.4,
		 T => -0.7,
		 W => -0.9,
		 S => -0.8,
		 Y => -1.3,
		 P => -1.6,
		 H => -3.2,
		 E => -3.5,
		 Q => -3.5,
		 D => -3.5,
		 N => -3.5,
		 K => -3.9,
		 R => -4.5,
	       );
		 
  my @datapoints;
  my @seq = split('', uc($protein));

  $kd_window = $kd_window < scalar(@seq) ? $kd_window : scalar(@seq);

  my $maxkd = 4.5;
  my $minkd = -4.5;

  my $kd = 0;
  for (my $i = 0 ; $i < @seq && $i < $kd_window ; $i++) {
    $kd += $scores{$seq[$i]} || 0;
  }

  my $content = $kd / $kd_window;
  push @datapoints, $content;

  for (my $i = $kd_window; $i < @seq; $i++) {
    $kd -= $scores{$seq[$i-$kd_window]} || 0;
    $kd += $scores{$seq[$i]} || 0;
    $content = $kd / $kd_window;
    push @datapoints, $content;
  }

  my $scale = $maxkd - $minkd;
  foreach (my $i = 0; $i < @datapoints; $i++) {
    $datapoints[$i] = ($datapoints[$i] - $minkd) / $scale;
  }

  # Calculate values that will be used in the layout
  
  my $bin_height = $y2-$y1;
  my $fgcolor    = $self->fgcolor;
  my $bgcolor    = $self->factory->translate_color($self->panel->gridcolor);
  my $axiscolor  = $self->color('axis_color') || $fgcolor;

  # Draw the axes
  
  $gd->line($x1,  $y1,        $x1,  $y2,        $axiscolor);
  $gd->line($x2-2,$y1,        $x2-2,$y2,        $axiscolor);
  $gd->line($x1,  $y1,        $x1+3,$y1,        $axiscolor);
  $gd->line($x1,  $y2,        $x1+3,$y2,        $axiscolor);
  $gd->line($x1,  ($y2+$y1)/2,$x1+3,($y2+$y1)/2,$axiscolor);
  $gd->line($x2-4,$y1,        $x2-1, $y1,       $axiscolor);
  $gd->line($x2-4,$y2,        $x2-1, $y2,       $axiscolor);
  $gd->line($x2-4,($y2+$y1)/2,$x2-1,($y2+$y1)/2,$axiscolor);
  $gd->line($x1+5,$y2,        $x2-5,$y2,        $bgcolor);
  $gd->line($x1+5,($y2+$y1)/2,$x2-5,($y2+$y1)/2,$bgcolor);
  $gd->line($x1+5,$y1,        $x2-5,$y1,        $bgcolor);
  $gd->string($self->font,$x1+5,$y1,'Kyte-Doolittle hydropathy plot',$axiscolor)
      if $bin_height > $self->font->height*2;

  $gd->string($self->font,$x2-20,$y1,$maxkd,$axiscolor) 
    if $bin_height > $self->font->height*2.5;
  $gd->string($self->font,$x2-20,$y2-$self->font->height,$minkd,$axiscolor) 
    if $bin_height > $self->font->height*2.5;

  my $graphwidth = $x2 - $x1;
  my $scale = $graphwidth / (@datapoints + $kd_window - 1);
  for (my $i = 1; $i < @datapoints; $i++) {
    my $x = $i + $kd_window / 2;
    my $xlo = $x1 + ($x - 1) * $scale;
    my $xhi = $x1 + $x * $scale;
    my $y = $y2 - ($bin_height*$datapoints[$i]);
    $gd->line($xlo, $y2 - ($bin_height*$datapoints[$i-1]), 
	      $xhi, $y, 
	      $fgcolor);
  }
}

sub make_key_feature {
  my $self = shift;
  my @gatc = qw(A C D E F G H I K L M N P Q R S T V W Y);
  my $offset = $self->panel->offset;
  my $scale = 1/$self->scale;  # base pairs/pixel

  my $start = $offset+1;
  my $stop  = $offset+100*$scale;
  my $feature =
    Bio::Graphics::Feature->new(-start=> $start,
				-stop => $stop,
				-seq  => join('',map{$gatc[rand 4]} (1..500)),
				-name => $self->option('key'),
				-strand => '+1',
			       );
  $feature;
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::protein - The "protein" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws protein sequences.  At high magnifications, this
glyph will draw the actual amino acids of the sequence.  At low
magnifications, the glyph will plot the Kyte-Doolite hydropathy.  By
default, the KD plot will use a window size of 9 residues, but this
can be changed by specifying the kd_window option.

For this glyph to work, the feature must return a protein sequence
string in response to the seq() method.

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

  -hilite       Highlight color                undef (no color)

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description               Default
  ------      -----------               -------

  -do_kd      Whether to draw the Kyte-  true
              Doolittle hydropathy plot
              at low mags

  -kd_window  Size of the sliding window  9
  	      to use in the KD hydropathy 
	      calculation.

  -axis_color Color of the vertical axes  fgcolor
              in the KD hydropathy plot


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

Aaron J. Mackey, based on the "dna" glyphy by Lincoln Stein
E<lt>lstein@cshl.orgE<gt> and Peter Ashton E<lt>pda@sanger.ac.ukE<gt>.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
