package Bio::Graphics::Glyph::cds;

use strict;
use Bio::Graphics::Glyph::segments;
use Bio::Graphics::Util qw(frame_and_offset);
use Bio::Tools::CodonTable;
use Bio::Graphics::Glyph::translation;
use vars '@ISA','$VERSION';
@ISA = qw(Bio::Graphics::Glyph::segments Bio::Graphics::Glyph::translation);
$VERSION = '1.01';

sub connector   { 0 };
sub description {
  my $self = shift;
  return if $self->level;
  return $self->SUPER::description;
};

# figure out (in advance) the color of each component
sub draw {
  my $self = shift;
  my ($gd,$left,$top) = @_;

  my @parts = $self->parts;

  return $self->SUPER::draw(@_) unless @parts;

  my $fits = $self->protein_fits;

  if (!$fits) {
    # draw the staff (musically speaking)
    my ($x1,$y1,$x2,$y2) = $self->bounds($left,$top);
    my $height = ($y2-$y1)/3;
    my $grid  = $self->gridcolor;
    for (0..2) {
      my $offset = $y1+$height*$_+1;
      $gd->line($x1,$offset,$x2,$offset,$grid);
    }
  }

  $self->{cds_part2color} ||= {};
  my $fill   = $self->bgcolor;
  my $strand = $self->feature->strand;

  # figure out the colors of each part
  # sort minus strand features backward
  @parts = sort {$b->left <=> $a->left} @parts if $strand < 0;
  my $translate_table = Bio::Tools::CodonTable->new;

  for (my $i=0; $i < @parts; $i++) {
    my $part    = $parts[$i];
    my $feature = $part->feature;
    my $pos     = $strand > 0 ? $feature->start : $feature->end;
    my $phase           = eval {$feature->phase} || 0;
    my ($frame,$offset) = frame_and_offset($pos,
					   $feature->strand,
					   -$phase);
    my $suffix = $strand < 0 ? 'r' : 'f';
    my $key = "frame$frame$suffix";
    $self->{cds_frame2color}{$key} ||= $self->color($key) || $fill;
    $part->{cds_partcolor} = $self->{cds_frame2color}{$key};
    $part->{cds_frame}     = $frame;
    $part->{cds_offset}    = $offset;

    next unless $fits;

    # do in silico splicing in order to find the codon that
    # arises from the splice
    my $protein = $part->feature->translate(undef,undef,$phase)->seq;
    $part->{cds_translation}  = $protein;

  BLOCK: {
      length $protein >= $feature->length/3           and last BLOCK;
      ($feature->length - $phase) % 3 == 0            and last BLOCK;

      my $next_part    = $parts[$i+1]
	or do {
	  $part->{cds_splice_residue} = '?';
	  last BLOCK; };

      my $next_feature = $next_part->feature         or  last BLOCK;
      my $next_phase   = eval {$next_feature->phase} or  last BLOCK;
      my $splice_codon = '';
      my $left_of_splice  = substr($feature->seq,-$next_phase,$next_phase);
      my $right_of_splice = substr($next_feature->seq,0,3-$next_phase);
      $splice_codon = $left_of_splice . $right_of_splice;
      length $splice_codon == 3                      or last BLOCK;
      my $amino_acid = $translate_table->translate($splice_codon);
      $part->{cds_splice_residue} = $amino_acid;
    }
  }

  $self->Bio::Graphics::Glyph::generic::draw($gd,$left,$top);
}


# draw the notes on the staff
sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->bounds(@_);

  my $color = $self->{cds_partcolor};
  my $feature   = $self->feature;
  my $frame     = $self->{cds_frame};

  unless ($self->protein_fits) {
    my $height = ($y2-$y1)/3;
    my $offset = $y1 + $height*$frame;
    $gd->filledRectangle($x1,$offset,$x2,$offset+2,$color);
    return;
  }

  # we get here if there's room to draw the primary sequence
  my $font  = $self->font;
  my $pixels_per_residue = $self->pixels_per_residue;
  my $strand = $feature->strand;

  # have to remap feature start and end into pixel coords in order to:
  # 1) correctly align the amino acids with the nucleotide seq
  # 2) correct for the phase offset
  my $start = $self->map_no_trunc($feature->start + $self->{cds_offset});
  my $stop  = $self->map_no_trunc($feature->end   + $self->{cds_offset});

  my @residues = split '',$self->{cds_translation};

  push @residues,$self->{cds_splice_residue} if $self->{cds_splice_residue};
  for (my $i=0;$i<@residues;$i++) {
    my $x = $strand > 0 ? $start + $i * $pixels_per_residue
                        : $stop  - $i * $pixels_per_residue;
    $gd->char($font,$x,$y1,$residues[$i],$color) if $x >= $x1 && $x <= $x2;
  }
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::cds - The "cds" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws features that are associated with a protein coding
region.  At high magnifications, draws a series of boxes that are
color-coded to indicate the frame in which the translation occurs.  At
low magnifications, draws the amino acid sequence of the resulting
protein.  Amino acids that are created by a splice are optionally
shown in a distinctive color.

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

In addition, the alignment glyph recognizes the following
glyph-specific options:

  Option      Description                  Default
  ------      -----------                  -------

  -frame0f    Color for first (+) frame    background color

  -frame1f    Color for second (+) frame   background color

  -frame2f    Color for third (+) frame    background color

  -frame0r    Color for first (-) frame    background color

  -frame1r    Color for second (-) frame   background color

  -frame2r    Color for third (-) frame    background color

  -gridcolor  Color for the "staff"        lightslategray

=head1 SUGGESTED STANZA FOR GENOME BROWSER

Using the "coding" aggregator, this produces a nice gbrowse display.

 [CDS]
 feature      = coding
 glyph        = cds
 frame0f      = cadetblue
 frame1f      = blue
 frame2f      = darkblue
 frame0r      = darkred
 frame1r      = red
 frame2r      = crimson
 description  = 0
 height       = 13
 label        = CDS frame
 key          = CDS
 citation     = This track shows CDS reading frames.

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
