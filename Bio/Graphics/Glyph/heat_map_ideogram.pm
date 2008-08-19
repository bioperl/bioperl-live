package Bio::Graphics::Glyph::heat_map_ideogram;

# $Id: heat_map_ideogram.pm,v 1.5 2006/10/18 01:46:14 sheldon_mckay Exp $
# Glyph to draw chromosome heat_map ideograms

use strict qw/vars refs/;
use GD;

use base qw(Bio::Graphics::Glyph::ideogram Bio::Graphics::Glyph::heat_map);

sub draw {
  my $self = shift;

  my @parts = $self->parts;

  @parts = $self if !@parts && $self->level == 0;
  return $self->SUPER::draw(@_) unless @parts;

  $self->{single}++ if @parts == 1;

  $self->calculate_gradient(\@parts);

  # adjust for label and description
  my ($gd,$x,$y) = @_;
  $x    += $self->left + $self->pad_left;
  $y += $self->top  + $self->pad_top;

  # Draw centromeres and telomeres last
  my @last;
  for my $part (@parts) {
    push @last, $part and next if 
	$part->feature->method eq 'centromere' ||
	$part->feature->start <= 1 ||
	$part->feature->stop >= $self->panel->end - 1000;

    $self->draw_component($part,$gd,$x,$y);
  }

  for my $part (@last) {
    my $tile = $self->create_tile('right') 
	if $part->feature->method eq 'centromere';
    $self->draw_component($part,$gd,$x,$y);
  }

  $self->draw_label(@_)       if $self->option('label');
  $self->draw_description(@_) if $self->option('description');
}

sub draw_component {
  my $self  = shift;
  my $glyph = shift;
  my $gd    = shift;
  my ( $x1, $y1, $x2, $y2 ) = $glyph->bounds(@_);
  # force odd width so telomere arcs are centered 
  $y2 ++ if ($y2 - $y1) % 2;
  
  my $arcradius = $self->option('arcradius') || 7;
  my $feature   = $glyph->feature;
  my $score     = $feature->score;
  my $is_cent   = 1 if $feature->method eq 'centromere';
  my $fgcolor   = $self->fgcolor;
  my $bgcolor;

  # skip normal cytobands
  return if $feature->attributes('stain') && !$is_cent;
     
  # Set the bgcolor
  unless ($is_cent || defined $score || defined $self->score_range ) {
    my @rgb = @{$self->low_rgb};
    $bgcolor = $self->color_index(@rgb);
  }
  else {
    my @rgb = $self->calculate_color($score);
    $bgcolor = $self->color_index(@rgb);
  }

  # bgcolorindex must return true
  $bgcolor ||= $self->adjust_bgcolor;

  # Is this a centromere?
  if ( $is_cent ) {
    $fgcolor = $self->color_index(0,0,0);

    if ( $self->panel->image_class =~ /SVG/ ) {
      $bgcolor  = $gd->colorResolve( 102, 102, 153 );
      $self->draw_centromere( $gd, $x1, $y1, $x2, $y2, $bgcolor, $fgcolor );
    }
    else {
      $self->draw_centromere( $gd, $x1, $y1, $x2, $y2, gdTiled, $fgcolor );
    }
  }
  # a telomere?
  elsif ( $feature->start <= 1 ) {
    # left (top)
    my $status = 1 unless $self->panel->flip;
    # if this is a full-length chromosome?
    $status = -1 if $feature->stop >= $self->panel->end - 1000;

    $self->draw_telomere( $gd, $x1, $y1, $x2, $y2, $bgcolor, $fgcolor,
			  $arcradius, $status );
  }
  elsif ( $feature->stop >= $self->panel->end - 1000 ) {
    # right (bottom)
    my $status = 1 if $self->panel->flip;
    $self->draw_telomere( $gd, $x1, $y1, $x2, $y2, $bgcolor, $fgcolor,
			  $arcradius, $status );
  }
  # or a regular band?
  else {
    $self->draw_cytoband( $gd, $x1, $y1, $x2, $y2, $bgcolor, $fgcolor );
  }
}


# Nudge the color over just a bit if the color index
# is 0 (panel bgcolor).  This overcomes default bgcolor
# and fgcolor when the index does not return true
sub adjust_bgcolor {
  my $self = shift;
  my $gd   = $self->panel->gd;
  my @rgb = $self->panel->rgb($self->panel->bgcolor);

  for (@rgb) {
    $_++ if $_  < 255;
    $_-- if $_ == 255;
  }
  
  return $gd->colorResolve(@rgb);
}

sub fgcolor {
  my $self = shift;
  my $clr  = $self->option('fgcolor') || 
             $self->option('outline') ||
             'black';
  return $self->panel->translate_color($clr);

}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::heat_map_ideogram - The "heat_map_ideogram" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws a chromosome ideogram using scored features instead
of cytobands.  It is a hybrid of the heat_map and ideograms glyphs
and accepts options for both.  A typical usage would be to pair this
glyph with an aggregator that groups scored features such as blast hits
or gene_density bins, etc with a centromere.  The result is a chromosome
ideogram that has bands whose colors vary porportionate to the feature
score.

=head2 OPTIONS

See L<Bio::Graphics::Glyph> for a full explanation of standard options.

See L<Bio::Graphics::Glyph::heat_map> for an explanation of heat_map options.

See L<Bio::Graphics::Glyph::ideogram> for an explanation of ideogram options.

=head1 BUGS

Please report them.

=head1 SEE ALSO

L<Bio::Graphics::Panel>,
L<Bio::Graphics::Glyph>,
L<Bio::DB::GFF>,
L<Bio::SeqI>,
L<Bio::SeqFeatureI>,
L<Bio::Das>,
L<GD>

=head1 AUTHOR

Sheldon McKay E<lt>mckays@cshl.eduE<gt>

Copyright (c) 2006 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut







