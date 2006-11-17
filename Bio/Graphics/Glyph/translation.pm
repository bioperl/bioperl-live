package Bio::Graphics::Glyph::translation;

use strict;
use Bio::Graphics::Util qw(frame_and_offset);
use base qw(Bio::Graphics::Glyph::generic);

my %default_colors = qw(
			frame0f  cornflowerblue
			frame1f  blue
			frame2f  darkblue
			frame0r  magenta
			frame1r  red
			frame2r  darkred
		       );

# turn off description
sub description { 0 }

# turn off label
# sub label { 1 }

sub default_color {
  my ($self,$key) = @_;
  return $self->factory->translate_color($default_colors{$key});
}

sub height {
  my $self = shift;
  my $font = $self->font;
  my $lines = $self->translation_type eq '3frame' ? 3
            : $self->translation_type eq '6frame' ? 6
            : 1;
  return $self->protein_fits ? $lines*$font->height
       : $self->SUPER::height;
}

sub pixels_per_base {
  my $self = shift;
  return $self->scale;
}

sub pixels_per_residue {
  my $self = shift;
  return $self->scale * 3;
}

sub gridcolor {
  my $self = shift;
  my $color = $self->option('gridcolor') || 'lightgrey';
  $self->factory->translate_color($color);
}

sub show_sequence {
  my $self = shift;
  my $show_sequence = $self->option('show_sequence');
  return 1 unless defined $show_sequence;  # default to true
  return $show_sequence;
}

sub triletter_code {
  my $self = shift;
  my $triletter_code = $self->option("triletter_code");
  return 0 unless defined $triletter_code; # default to false
  return $triletter_code;
}

sub longprotein_fits {
  my $self = shift;
  return unless $self->show_sequence;

  my $pixels_per_residue = $self->pixels_per_residue;
  my $font               = $self->font;
  my $font_width         = $font->width * 4; # not 3; leave room for whitespace

  return $pixels_per_residue >= $font_width;
}

sub translation_type {
  my $self = shift;
  return $self->option('translation') || '1frame';
}

sub arrow_height {
  my $self = shift;
  $self->option('arrow_height') || 1;
}

sub show_stop_codons {
  my $self = shift;
  my $show = $self->option('stop_codons');
  return $show if defined $show;
  return 1;
}

sub show_start_codons {
  my $self = shift;
  my $show = $self->option('start_codons');
  return $show if defined $show;
  return 0;
}

sub strand {
  my $self = shift;
  return $self->option('strand') || '+1';
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->bounds(@_);

  my $type   = $self->translation_type;
  my $strand = $self->strand;

  my @strands =  $type eq '6frame' ? (1,-1)
	       : $strand > 0       ? (1)
	       : -1;
  my @phase = (0,1,2);
  for my $s (@strands) {
    for (my $i=0; $i < @phase; $i++) {
      $self->draw_frame($self->feature,$s,$i,$phase[$i],$gd,$x1,$y1,$x2,$y2);
    }
  }

}

sub draw_frame {
  my $self = shift;
  my ($feature,$strand,$base_offset,$phase,$gd,$x1,$y1,$x2,$y2) = @_;
  my ($seq,$pos);
  $seq = $feature->seq or return; # no sequence, arggh.

  my $strand0 = $strand;
  $strand *= -1 if $self->{flip};

  $pos = $strand < 0 ? $feature->end : $feature->start;

  my ($frame,$offset) = frame_and_offset($pos,$strand,$phase);
  # warn "frame=$frame, phase=$phase";

  my ($x1_orig,$x2_orig) = ($x1,$x2);  # remember this for arrowheads

  ($strand >= 0 ? $x1 : $x2) += $self->pixels_per_base * $offset;
  my $y0 = $y1;
  my $lh;
  if ($self->translation_type eq '6frame') {
    $lh = $self->height / 6;
    $y1 += $lh * $frame;
    $y1 += $self->height/2 if $strand < 0;
  } else {
    $lh = $self->height / 3;
    $y1 += $lh * $frame;
  }

  $y1  = $y0 + ($self->height - ($y1-$y0)) - $lh if $self->{flip};

  $y2 = $y1;

  my $codon_table = $self->option('codontable') || $self->option('geneticcode') || 1;

  # the dreaded difference between a Bio::SeqFeature and a Bio::Seq

  my $realseq  = $self->get_seq($seq);
  return unless $realseq;
  $realseq    = $realseq->revcom if $strand < 0;

  my $protein = $realseq->translate(undef,undef,$base_offset,$codon_table)->seq;

  my $k       = $strand >= 0     ? 'f' : 'r';

  my $color   = $self->color("frame$frame$k") ||
                $self->color("frame$frame") ||
                $self->default_color("frame$frame$k") || $self->fgcolor;

  my $awo = 0;
  if ($self->protein_fits) {
    $self->draw_protein(\$protein,$strand,$color,$gd,$x1,$y1,$x2,$y2);
    $awo += $self->font->height/2;
  } else {
    $self->draw_orfs(\$protein,$strand,$color,$gd,$x1,$y1,$x2,$y2);
  }

  $strand0 > 0 ? $self->arrowhead($gd,$x2_orig+5,$y1+$awo,3,+1)
               : $self->arrowhead($gd,$x1_orig-5,$y1+$awo,3,-1)

}

sub draw_protein {
  my $self = shift;
  my ($protein,$strand,$color,$gd,$x1,$y1,$x2,$y2) = @_;
  my $pixels_per_base = $self->pixels_per_base;
  my $font   = $self->font;
  my $flip   = $self->{flip};
  my $left   = $self->panel->left;
  my $right  = $self->panel->right;

  my $longprotein = $self->triletter_code && $self->longprotein_fits;

  my %abbrev = ( A => "Ala", B => "Asx", C => "Cys", D => "Asp",
		 E => "Glu", F => "Phe", G => "Gly", H => "His",
		 I => "Ile", J => "???", K => "Lys", L => "Leu",
		 M => "Met", N => "Asn", O => "???", P => "Pro",
		 Q => "Gln", R => "Arg", S => "Ser", T => "Thr",
		 U => "Sec", V => "Val", W => "Trp", X => "Xaa",
		 Y => "Tyr", Z => "Glx", '*' => " * ",
	       );

  my @residues = split '',$$protein;
  my $fontwidth = $font->width;
  for (my $i=0;$i<@residues;$i++) {
    my $x = $strand > 0
      ? $x1 + 3 * $i * $pixels_per_base
      : $x2 - 3 * $i * $pixels_per_base - $pixels_per_base;
    next if $x+1 < $x1;
    last if $x > $x2;
    if ($flip) {
      $x -= $pixels_per_base - $font->width - 1; #align right, not left
      if ($longprotein) {
	$gd->string($font,$right-($x-$left+$pixels_per_base)+1,$y1,$abbrev{$residues[$i]},$color);
      } else {
	$gd->char($font,$right-($x-$left+$pixels_per_base)+2,$y1,$residues[$i],$color);
      }
    } else {
      if ($longprotein) {
	$gd->string($font, $x+1, $y1, $abbrev{$residues[$i]}, $color);
      } else {
	$gd->char($font,$x+2,$y1,$residues[$i],$color);
      }
    }
  }
}

sub draw_orfs {
  my $self     = shift;
  my ($protein,$strand,$color,$gd,$x1,$y1,$x2,$y2) = @_;
  my $pixels_per_base = $self->pixels_per_base * 3;
  $y1++;
  my $right  = $self->panel->right;
  my $left   = $self->panel->left;
  my $flip   = $self->{flip};

  my $gcolor = $self->gridcolor;
  $gd->line($x1,$y1,$x2,$y1,$gcolor);

  if ($self->show_stop_codons) {
    my $stops  = $self->find_codons($protein,'*');

    for my $stop (@$stops) {
      my $pos = $strand > 0 
	? $x1 + $stop * $pixels_per_base
        : $x2 - $stop * $pixels_per_base;
      next if $pos+1 < $x1;
      last if $pos   > $x2;
      if ($flip) {
	$gd->line($right-($pos-$left),$y1-2,$right-($pos-$left),$y1+2,$color);
      } else {
	$gd->line($pos,$y1-2,$pos,$y1+2,$color);
      }
    }
  }

  my $arrowhead_height = $self->arrow_height;

  if ($self->show_start_codons) {
    my $starts  = $self->find_codons($protein,'M');

    for my $start (@$starts) {
      my $pos = $strand > 0 
	? $x1 + $start * $pixels_per_base
        : $x2 - $start * $pixels_per_base;
      next if $pos+1 < $x1;
      last if $pos   > $x2;
      $pos = $self->{flip} ? $right - $pos : $pos;

      # little arrowheads at the start codons
      $strand > 0 ? $self->arrowhead($gd,$pos-$arrowhead_height,$y1,
				     $arrowhead_height,+1)
	          : $self->arrowhead($gd,$pos+$arrowhead_height,$y1,
				     $arrowhead_height,-1)
    }
  }
  $strand *= -1 if $flip;

}

sub find_codons {
  my $self    = shift;
  my $protein = shift;
  my $codon   = shift || '*';
  my $pos = -1;
  my @stops;
  while ( ($pos = index($$protein,$codon,$pos+1)) >= 0) {
    push @stops,$pos;
  }
  \@stops;
}

sub make_key_feature {
  my $self = shift;
  my @gatc = qw(g a t c);
  my $offset = $self->panel->offset;
  my $scale = 1/$self->scale;  # base pairs/pixel
  my $start = $offset;
  my $stop  = $offset + 100 * $scale;
  my $seq   = join('',map{$gatc[rand 4]} (1..500));
  my $feature =
    Bio::Graphics::Feature->new(-start=> $start,
				-end  => $stop,
				-seq  => $seq,
				-name => $self->option('key')
			       );
  $feature;
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::translation - The "6-frame translation" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws the conceptual translation of DNA sequences.  At high
magnifications, it simply draws lines indicating open reading frames.
At low magnifications, it draws a conceptual protein translation.
Options can be used to set 1-frame, 3-frame or 6-frame translations.

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

  Option        Description                 Default
  ------        -----------                 -------

  -translation  Type of translation to      1frame
                perform.  One of "1frame",
                "3frame", or "6frame"

  -strand       Forward (+1) or reverse (-1) +1
                translation.

  -frame0       Color for the first frame    fgcolor

  -frame1       Color for the second frame   fgcolor

  -frame2       Color for the third frame    fgcolor

  -gridcolor    Color for the horizontal     lightgrey
                lines of the reading frames

  -start_codons Draw little arrowheads       0 (false)
                indicating start codons

  -stop_codons  Draw little vertical ticks   1 (true)
                indicating stop codons

  -arrow_height Height of the start codon    1
                arrowheads

  -show_sequence Show the amino acid sequence 1 (true)
                if there's room.

  -triletter_code Show the 3-letter amino acid 0 (false)
                code if there's room

  -codontable   Codon table to use           1 (see Bio::Tools::CodonTable)

=head1 SUGGESTED STANZA FOR GENOME BROWSER

This produces a nice gbrowse display in which the DNA/GC Content glyph
is sandwiched between the forward and reverse three-frame
translations.  The frames are color-coordinated with the example
configuration for the "cds" glyph.

 [TranslationF]
 glyph        = translation
 global feature = 1
 frame0       = cadetblue
 frame1       = blue
 frame2       = darkblue
 height       = 20
 fgcolor      = purple
 strand       = +1
 translation  = 3frame
 key          = 3-frame translation (forward)

 [DNA/GC Content]
 glyph        = dna
 global feature = 1
 height       = 40
 do_gc        = 1
 fgcolor      = red
 axis_color   = blue

 [TranslationR]
 glyph        = translation
 global feature = 1
 frame0       = darkred
 frame1       = red
 frame2       = crimson
 height       = 20
 fgcolor      = blue
 strand       = -1
 translation  = 3frame
 key          = 3-frame translation (reverse)

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
