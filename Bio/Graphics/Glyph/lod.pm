package Bio::Graphics::Glyph::lod;

use strict;
use Bio::Graphics::Glyph::generic;
use vars '@ISA', '$VERSION';
@ISA = qw(Bio::Graphics::Glyph::generic);
$VERSION = 1.00;

# turn off description
# sub description { 0 }

# turn off label
# sub label { 0 }

# we are unbumpable!
sub bump {
  return 0;
}

sub pad_left {
  my $self = shift;
  my $left = $self->SUPER::pad_left;
  return $left unless ( $self->option('point') && ( $self->level != -1 ) );
  my $extra = $self->option('height')/3;
  return $extra > $left ? $extra : $left;
}

sub pad_right {
  my $self = shift;
  my $right = $self->SUPER::pad_right;
  return $right unless ( $self->option('point') && ( $self->level != -1 ) );
  my $extra = $self->option('height')/3;
  return $extra > $right ? $extra : $right;
}

sub height {
  my $self = shift;
  return ( ( $self->level == -1 ) ?
           $self->option('graph_height') :
           $self->SUPER::height );
}

#sub width {
#  my $self = shift;
#  return ( ( $self->level == -1 ) ?
#           $self->SUPER::width :
#           $self->option('height') || $self->SUPER::width );
#}

# track sets connector to empty
sub connector {
  my $self = shift;
  if( $self->level == -1 ) {
    return $self->SUPER::connector(@_) if $self->all_callbacks;
    return 'none';
  } else {
    return $self->SUPER::connector( @_ );
  }
}

## We override this to do the simplest thing, as if this track glyph
## doesn't have any parts.  That's because we do our own kind of
## layout in the draw(..) method.
sub layout {
  my $self = shift;
  return $self->{layout_height} if exists $self->{layout_height};

  if( $self->level == -1 ) {
    return
      $self->{layout_height} =
        $self->height + $self->pad_top + $self->pad_bottom;
  }
  return $self->SUPER::layout( @_ );
}

sub draw {
  my $self = shift;
  my ($gd,$dx,$dy,$partno,$total_parts) = @_;
  my ($left,$top,$right,$bottom) = $self->calculate_boundaries($dx,$dy);

  if( $self->level() == -1 ) {
    my @parts = $self->parts;
    my $divisor = ( $self->option( 'lod_score_divisor' ) || 1 );
    my $max = $self->option( 'max_lod_score' ) || 1;
    my $min = $self->option( 'min_lod_score' ) || -1;
    ## TODO: Calculate max & min if they're not given?
    my ( $score, $score_percent, $score_height );
    for( my $part_i = 0; $part_i < @parts; $part_i++ ) {
      $score = ( $parts[ $part_i ]->score() / $divisor );
      ## TODO: REMOVE
      #warn "part score is $score";
      $score_percent = ( ( $score - $min ) / ( $max - $min ) );
      $score_height  = ( ( $bottom - $top ) * $score_percent );
      $parts[ $part_i ]->{ '_y_position' } = $bottom - $score_height;
    } # End for each part, set its _y_position
    $self->draw_lod_graph($gd,$left,$top,$right,$bottom);
    for( my $part_i = 0; $part_i < scalar( @parts ); $part_i++ ) {
      ## Don't draw the offscreen ones.
      next if( $parts[ $part_i ]->{'_offscreen'} );

      $parts[ $part_i ]->draw_component( $gd, $dx, $dy, 0, 1 );
    } # End for each part, draw it.
  } else {
    $self->throw( "lod::draw(..) should only be called on the track level (-1).  It was called on a glyph with level ".$self->level()."." );
  }
} # draw(..)

## Only draw the key label
sub draw_label {
  my $self = shift;
  
  if( $self->label !~ /LOD/ ) {
    return;
  }
  $self->SUPER::draw_label( @_ );
  return;
} # draw_label

sub score {
  my $self = shift;
  if( $self->level == -1 ) {
    return 0;
  }
  $self->feature->score;
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ( $x1, $y1, $x2, $y2 ) = $self->calculate_boundaries( @_ );

  ## TODO: REMOVE
  #warn "lod::draw_component(..): level is ".$self->level()." feature is ".$self->feature().", and its subfeatures are ( ".join( ', ', $self->feature()->features() )." )";

  my $y = $self->{ '_y_position' };
  $self->draw_microsatellite($gd,$x1,$y,$x2,$y);
}

sub draw_lod_graph {
  my $self     = shift;
  my $gd       = shift;

  my ($left,$top,$right,$bottom) = @_;

  my $fgcolor    = $self->fgcolor;
  my $bgcolor    = $self->factory->translate_color($self->panel->gridcolor);
  my $axis_fgcolor  = $self->color('axis_fgcolor') || $fgcolor;
  my $axis_bgcolor  = $self->color('axis_bgcolor') ||
    $self->factory->translate_color($self->panel->gridcolor);
  my $axis_text_color = $self->color('axis_text_color') || $axis_fgcolor;
  my $zeroaxis_fgcolor  = $self->color('zero_axis_fgcolor') || $axis_fgcolor;
  my $zeroaxis_bgcolor  = $self->color('zero_axis_bgcolor') || $axis_bgcolor;

  my $max = $self->option( 'max_lod_score' ) || 1;
  my $min = $self->option( 'min_lod_score' ) || -1;
  my $zero_percent = ( ( 0 - $min ) / ( $max - $min ) );
  my $zero_height  = ( ( $bottom - $top ) * $zero_percent );

  # Draw the axis
  # Left
  $gd->line($left,  $top,        $left,  $bottom,        $axis_fgcolor);
  # Right
  $gd->line($right-1,$top,        $right-1,$bottom,        $axis_fgcolor);
  # Top-Left axis
  #$gd->line($left,  $top,        $left+3,$top,        $axis_fgcolor);
  # Top-Right axis
  #$gd->line($right-4,$top,        $right-1, $top,       $axis_fgcolor);
  # Top bgcolor axis
  #$gd->line($left+5,$top,        $right-5,$top,        $axis_bgcolor);
  # Bottom-Right axis
  #$gd->line($right-4,$bottom,        $right-1, $bottom,       $axis_fgcolor);
  # Bottom-Left axis
  #$gd->line($left,  $bottom,        $left+3,$bottom,        $axis_fgcolor);
  # Bottom bgcolor axis
  #$gd->line($left+5,$bottom,        $right-5,$bottom,        $axis_bgcolor);
  ## TODO: Make this, too, configurable ( draw every whole number versus every two, or every half, or what-have-you. )
  my $half_font_height = 1 + ( .5 * $self->font()->height() );
  my ( $wn_percent, $wn_height );
  for( my $whole_number =
       ( ( int( $min ) < $min ) ? int( $min ) + 1 : int( $min ) ); # ceil($min)
       $whole_number < $max;
       $whole_number++ ) {
    if( $whole_number == 0 ) {
      # Zero-Left axis
      $gd->line($left,$bottom-$zero_height,$left+3,$bottom-$zero_height,$zeroaxis_fgcolor);
      # Zero-Right axis
      $gd->line($right-4,$bottom-$zero_height,$right-1,$bottom-$zero_height,$zeroaxis_fgcolor);
      # Zero bgcolor axis
      $gd->line($left+10,$bottom-$zero_height,$right-5,$bottom-$zero_height,$zeroaxis_bgcolor);
      # Zero axis label
      if( $zero_height > 0 ) {
        $gd->string($self->font,$left+5,$bottom-($zero_height+$half_font_height),0,$axis_text_color);
      }
    } else {
      $wn_percent = ( ( $whole_number - $min ) / ( $max - $min ) );
      $wn_height  = ( ( $bottom-$top ) * $wn_percent );
      if( $axis_text_color == $axis_fgcolor ) {
        # Whole Number bgcolor axis
        $gd->line($left+5,$bottom-$wn_height,$right-5,$bottom-$wn_height,$axis_bgcolor);
      } else {
        # Whole Number-Left axis
        $gd->line($left,  $bottom-$wn_height,$left+3,$bottom-$wn_height,$axis_fgcolor);
        # Whole Number-Right axis
        $gd->line($right-4,$bottom-$wn_height,$right-1,$bottom-$wn_height,$axis_fgcolor);
        # Whole Number bgcolor axis
        $gd->line($left+5,$bottom-$wn_height,$right-5,$bottom-$wn_height,$axis_bgcolor);
      }
      # Whole Number axis label
      $gd->string($self->font,$left+5,$bottom-($wn_height+$half_font_height),$whole_number,$axis_text_color);
    }
  }

  # Draw the track label
  $gd->string($self->font,$left+15,$top,$self->option('lod_label'),$axis_fgcolor);

  my @parts = $self->parts;
  return unless( @parts );

  # Prepare to draw the overlay line, if the overlay_line_color is defined.
  my ( $overlay_line_color, $gd_width, $gd_height );
  if( $self->option('overlay_line_color') ) {
    $overlay_line_color  = $self->color('overlay_line_color');
    ( $gd_width, $gd_height ) = $gd->getBounds();
  }

  my $last_x = $left;
  ## Leave $last_y undefined if you don't want the left line to go to
  ## 0 when there's nothing offscreen to the left.
  my $last_y;# = ( $bottom - $zero_height );
  my ( $x, $y );
  for my $part (sort { $a->left <=> $b->left } @parts) {
    ## TODO: REMOVE
    #warn "Part is at ".$part->feature->toRelRangeString( 1 );
    $x = $part->left + ( $part->width / 2 );
    $y = $part->{'_y_position'};
    ## Why do we need this next line?
    #next unless( $x >= $last_x );
    #last if( $x > $right );
    if( defined $last_y ) {
      $gd->line( $last_x, $last_y, $x, $y, $fgcolor );
    }
    # Draw the overlay line, if the overlay_line_color is defined.
    if( $self->option('overlay_line_color') ) {
      $gd->line( $x, 0, $x, $gd_height, $overlay_line_color );
    }
    $last_x = $x;
    $last_y = $y;
  }
  ## TODO: REMOVE
  #warn "Done with draw_lod_graph(..)";
} # draw_lod_graph(..)

sub draw_microsatellite {
  my $self     = shift;
  my $gd       = shift;
  my ($x1,$y1,$x2,$y2) = @_;

  ## TODO: Make outline color settable in the options, but ensure that
  ## it's not the same as the fgcolor!!!
  my $outline_color = $self->factory->translate_color('black');
  my $orient = $self->option('orient') || 'S';

  # find the center and vertices
  my $xmid = ($x1+$x2)/2;
  my $ymid = ($y1+$y2)/2;

  my ($vx1,$vy1,$vx2,$vy2,$vx3,$vy3);

  #make an equilateral
  my $radius = $self->option('height');
  # Round using sprintf.
  my $half_width  = sprintf( '%d', ( $radius / 2 ) );
  my $half_height = sprintf( '%d', ( ( $half_width * sqrt( 3 ) ) / 2 ) );
  $x1 = $xmid - $half_width;
  $x2 = $xmid + $half_width;
  $y1 = $ymid - $half_height;
  $y2 = $ymid + $half_height;

  ## TODO: REMOVE?  Testing fix to boxes().
  $self->{'box_left'}  = $x1;
  $self->{'box_right'} = $x2;
  $self->{'box_top'}   = $y1;
  $self->{'box_bottom'}= $y2;

  ## TODO: REMOVE
  #warn "DRAWING TRIANGLE xmid = $xmid, width = $radius, ymid = $ymid, height = $radius.  Color is $outline_color.";

  if   ($orient eq 'S'){$vx1=$x1;$vy1=$y1;$vx2=$x2;$vy2=$y1;$vx3=$xmid;$vy3=$y2;}
  elsif($orient eq 'N'){$vx1=$x1;$vy1=$y2;$vx2=$x2;$vy2=$y2;$vx3=$xmid;$vy3=$y1;}
  elsif($orient eq 'W'){$vx1=$x2;$vy1=$y1;$vx2=$x2;$vy2=$y2;$vx3=$xmid;$vy3=$ymid;}
  elsif($orient eq 'E'){$vx1=$x1;$vy1=$y1;$vx2=$x1;$vy2=$y2;$vx3=$xmid;$vy3=$ymid;}

  # now draw the triangle
  $gd->line($vx1,$vy1,$vx2,$vy2,$outline_color);
  $gd->line($vx2,$vy2,$vx3,$vy3,$outline_color);
  $gd->line($vx3,$vy3,$vx1,$vy1,$outline_color);

  ## TODO: Make the fuzzy & regular fillcolors options from the .conf file.

  my $fill_color;
  if( $self->feature->has_tag( 'fuzzy' ) &&
      $self->feature->get_tag_values( 'fuzzy' ) ) {
    $fill_color = $self->factory->translate_color('white');
  } else {
    $fill_color = $self->factory->translate_color('red'); #= $self->bgcolor;
  }
  if( $fill_color ) {
    #$gd->fill($xmid,$ymid,$fill_color);
    $gd->fillToBorder($xmid,$ymid,$outline_color,$fill_color) if $orient eq 'S' || $orient eq 'N';
    $gd->fillToBorder($x1+1,$ymid,$outline_color,$fill_color) if $orient eq 'E';
    $gd->fillToBorder($x2-1,$ymid,$outline_color,$fill_color) if $orient eq 'W';
  }

} # draw_microsatellite

# return boxes surrounding each part
sub boxes {
  my $self = shift;

  unless( $self->level == -1 ) {
    $self->throw( "lod::boxes(..) should only be called on the track level (-1).  It was called on a glyph with level ".$self->level()."." );
  }

  my @result;
  my $parts = $self->parts or return;
  my $part;
  for( my $part_i = 0; $part_i < scalar( @$parts ); $part_i++ ) {
    $part = $parts->[ $part_i ];
    push( @result,
          [ $part->feature(),
            $part->{'box_left'},
            $part->{'box_top'},
            $part->{'box_right'},
            $part->{'box_bottom'} ] );
  }
  return wantarray ? @result : \@result;
} # boxes(..)

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::lod - The Microsattelite LOD score track/glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

## TODO: Update the docs.  These are DNA docs, not lod docs.

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
