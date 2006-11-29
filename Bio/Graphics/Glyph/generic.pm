package Bio::Graphics::Glyph::generic;

use strict;
use Bio::Graphics::Util qw(frame_and_offset);
use base qw(Bio::Graphics::Glyph);

my %complement = (g=>'c',a=>'t',t=>'a',c=>'g',
		  G=>'C',A=>'T',T=>'A',C=>'G');

# new options are 'label'       -- short label to print over glyph
#                 'description'  -- long label to print under glyph
# label and description can be flags or coderefs.
# If a flag, label will be taken from seqname, if it exists or primary_tag().
#            description will be taken from source_tag().

sub height {
  my $self = shift;
  my $h    = $self->SUPER::height;
  return $h unless
    $self->option('draw_translation') && $self->protein_fits
      or
	$self->option('draw_dna') && $self->dna_fits;
  my $fh = $self->font->height + 2;
  return $h > $fh ? $h : $fh;
}

sub pad_top {
  my $self = shift;
  my $top  = $self->option('pad_top');
  return $top if defined $top;
  my $pad = $self->SUPER::pad_top;
  $pad   += $self->labelheight if $self->label && $self->label_position eq 'top';
  $pad;
}
sub pad_bottom {
  my $self = shift;
  my $bottom  = $self->option('pad_bottom');
  return $bottom if defined $bottom;
  my $pad = $self->SUPER::pad_bottom;
  $pad   += $self->labelheight if $self->description;
  $pad   += $self->labelheight if $self->part_labels && $self->label_position eq 'top';
  $pad;
}
sub pad_right {
  my $self = shift;
  my $pad = $self->SUPER::pad_right;
  my $label_width       = $self->label_position eq 'top' ? $self->labelwidth : 0;
  my $description_width = $self->descriptionwidth;
  my $max = $label_width > $description_width ? $label_width : $description_width;
  my $right = $max - $self->width;
  return $pad > $right ? $pad : $right;
}
sub pad_left {
  my $self = shift;
  my $pad = $self->SUPER::pad_left;
  return $pad unless $self->label_position eq 'left' && $self->label;
  $pad += $self->labelwidth;
  $pad;
}
sub labelfont {
  my $self = shift;
  return $self->getfont('label_font',$self->font);
}
sub descfont {
  my $self = shift;
  return $self->getfont('desc_font',$self->font);
}
sub labelwidth {
  my $self = shift;
  return $self->{labelwidth} ||= length($self->label||'') * $self->font->width;
}
sub descriptionwidth {
  my $self = shift;
  return $self->{descriptionwidth} ||= length($self->description||'') * $self->font->width;
}
sub labelheight {
  my $self = shift;
  return $self->{labelheight} ||= $self->font->height;
}
sub label_position {
  my $self = shift;
  return $self->{labelposition} ||= $self->option('label_position') || 'top';
}
sub label {
  my $self = shift;
  return if $self->{overbumped};  # set by the bumper when we have hit bump limit
  return unless $self->subpart_callbacks;  # returns true if this is level 0 or if subpart callbacks allowed
  return $self->_label if $self->{level} >= 0;
  return exists $self->{label} ? $self->{label}
                               : ($self->{label} = $self->_label);
}
sub description {
  my $self = shift;
  return if $self->{overbumped}; # set by the bumper when we have hit bump limit
  return unless $self->subpart_callbacks;  # returns true if this is level 0 or if subpart callbacks allowed
  return $self->_description if $self->{level} > 0;
  return exists $self->{description} ? $self->{description}
                                     : ($self->{description} = $self->_description);
}

sub part_labels {
  my $self = shift;
  my @parts = $self->parts;
  return ($self->{level} == 0) && @parts && @parts>1 && $self->option('part_labels');
}

sub part_label_merge {
  shift->option('part_label_merge');
}

sub maxdepth {
  my $self = shift;
  my $maxdepth =  $self->option('maxdepth');
  return $maxdepth if defined $maxdepth;
  return 1;
}

sub _label {
  my $self = shift;

  # allow caller to specify the label
  my $label = $self->option('label');

  return unless defined $label;
  return "1"    if $label eq '1 '; # 1 with a space
  return $label unless $label eq '1';

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
  return "1"   if $label eq '1 ';
  return $label unless $label eq '1';

  return $self->{_description} if exists $self->{_description};
  return $self->{_description} = $self->get_description($self->feature);
}

sub get_description {
  my $self = shift;
  my $feature = shift;

  # common places where we can get descriptions
  return join '; ',$feature->notes if $feature->can('notes');
  return $feature->desc            if $feature->can('desc');

  if ($feature->can('has_tag')) {
    return join '; ',$feature->get_tag_values('note')        if $feature->has_tag('note');
    return join '; ',$feature->get_tag_values('description') if $feature->has_tag('description');
  }

  my $tag = $feature->source_tag;
  return if $tag eq '';
  $tag;
}

sub draw {
  my $self = shift;
  my ($gd,$left,$top,$partno,$total_parts) = @_;

  local($self->{partno},$self->{total_parts});
  @{$self}{qw(partno total_parts)} = ($partno,$total_parts);

  $self->calculate_cds()      if $self->option('draw_translation') && $self->protein_fits;

  $self->SUPER::draw(@_);
  $self->draw_label(@_)       if $self->option('label');
  $self->draw_description(@_) if $self->option('description');
  $self->draw_part_labels(@_) if $self->option('label') && $self->option('part_labels');
}

sub draw_component {
  my $self = shift;
  $self->SUPER::draw_component(@_);
  $self->draw_translation(@_) if $self->{cds_translation}; # created earlier by calculate_cds()
  $self->draw_sequence(@_)    if $self->option('draw_dna') && $self->dna_fits;
}

# mostly stolen from cds.pm -- draw the protein translation
sub draw_translation {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->bounds(@_);

  my $feature = $self->feature;
  my $strand = $feature->strand;

  my $font    = $self->font;
  my $pixels_per_residue = $self->scale * 3;

  my $y         = $y1 + ($self->height - $font->height)/2;
  my $fontwidth = $font->width;
  my $color     = $self->fontcolor;

  $strand *= -1 if $self->{flip};

  # have to remap feature start and end into pixel coords in order to:
  # 1) correctly align the amino acids with the nucleotide seq
  # 2) correct for the phase offset
  my $start = $self->map_no_trunc($feature->start + $self->{cds_offset});
  my $stop  = $self->map_no_trunc($feature->end   + $self->{cds_offset});

  ($start,$stop) = ($stop,$start) if $stop < $start;  # why does this keep happening?
  my $x_fudge    = $self->{flip} ? 1 : 2;
  my $right      = $self->panel->right;
  my $left       = $self->panel->left;

  my @residues = split '',$self->{cds_translation};
  push @residues,$self->{cds_splice_residue} if $self->{cds_splice_residue};
  for (my $i=0;$i<@residues;$i++) {
    my $x = $strand > 0 ? $start + $i * $pixels_per_residue
                        : $stop  - $i * $pixels_per_residue;
    next unless ($x >= $x1 && $x <= $x2);
    $x -= $fontwidth + 1 if $self->{flip}; # align right when flipped
    last if $x+$fontwidth >= $right;
    last if $x            <= $left;
    $gd->char($font,$x+$x_fudge,$y,$residues[$i],$color);
  }
}

sub draw_sequence {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->bounds(@_);

  my $feature = $self->feature;
  my $strand = $feature->strand;

  my $font            = $self->font;
  my $pixels_per_base = $self->scale;

  my $y         = $y1 + ($self->height - $font->height)/2 - 1;
  my $fontwidth = $font->width;
  my $color     = $self->fontcolor;

  $strand *= -1 if $self->{flip};

  # have to remap feature start and end into pixel coords in order to:
  my $start = $self->map_no_trunc($feature->start);
  my $stop  = $self->map_no_trunc($feature->end);

  ($start,$stop) = ($stop,$start) if $stop < $start;  # why does this keep happening?
  my $x_fudge    = $self->{flip} ? 1 : 2;
  my $right      = $self->panel->right;
  my $left       = $self->panel->left;

  my $seq   = $self->get_seq($self->feature->seq);
  $seq      = $seq->seq if $seq;   # get the dna

  my $canonical = $self->option('canonical_strand');

  my @bases = split '',$seq;
  for (my $i=0;$i<@bases;$i++) {
    my $x = $strand >= 0 ? $start + $i * $pixels_per_base
                         : $stop  - $i * $pixels_per_base;
    next unless ($x >= $x1 && $x <= $x2);
    $x -= $fontwidth + 1 if $self->{flip}; # align right when flipped
    if ($strand >= 0) {
      last if $x + $fontwidth > $right;
    } else {
      next if $x >= $right;
      last if $x < $left;
    }
    my $base = $self->{flip} ? $complement{$bases[$i]} : $bases[$i];
    $base    = $complement{$base} if $canonical && $strand < 0;
    $gd->char($font,$x+$x_fudge,$y,$base,$color);
  }
}

sub min { $_[0] <= $_[1] ? $_[0] : $_[1] }
sub max { $_[0] >= $_[1] ? $_[0] : $_[1] }

sub draw_label {
  my $self = shift;
  my ($gd,$left,$top,$partno,$total_parts) = @_;
  my $label = $self->label or return;
  my $x    = $self->left + $left; # valid for both "top" and "left" because the left-hand side is defined by pad_left
  my $font = $self->labelfont;
  if ($self->label_position eq 'top') {
    $x += $self->pad_left;  # offset to beginning of the drawn part of the feature
    $x = $self->panel->left + 1 if $x <= $self->panel->left;
    $gd->string($font,
		$x,
		$self->top + $top - 1,
		$label,
		$self->fontcolor);
  }
  elsif ($self->label_position eq 'left') {
    $gd->string($font,
		$x,
		$self->{top} + ($self->height - $font->height)/2 + $top,
		$label,
		$self->fontcolor);
  }
}
sub draw_description {
  my $self = shift;
  my ($gd,$left,$top,$partno,$total_parts) = @_;
  my $label = $self->description or return;
  my $x = $self->left + $left;
  $x   += $self->pad_left;  # offset to beginning of drawn part of feature
  $x = $self->panel->left + 1 if $x <= $self->panel->left;
  my $dy= $self->part_labels ? $self->font->height : 0;
  $gd->string($self->descfont,
	      $x,
	      $self->bottom - $self->pad_bottom + $top + $dy,
	      $label,
	      $self->font2color);
}

sub draw_part_labels {
  my $self = shift;
  my ($gd,$left,$top,$partno,$total_parts) = @_;
  return unless $self->{level} == 0;
  my @p = $self->parts or return;
  @p > 1 or return;
  @p = reverse @p if $self->flip;

  my $font  = $self->font;
  my $width = $font->width;
  my $color = $self->fontcolor;

  my $y     = $top + $self->bottom - $self->pad_bottom;
  my $merge_em = $self->part_label_merge;

  my @parts;
  my $previous;

  if ($merge_em) {
    my $current_contig = [];

    for my $part (@p) {
      if (!$previous || $part->feature->start - $previous->feature->end <= 1) {
	push @$current_contig,$part;
      } else {
	push @parts,$current_contig;
	$current_contig = [$part];
      }
      $previous = $part;
    }
    push @parts,$current_contig;
  }

  else {
    @parts = map {[$_]} @p;
  }

  my $last_x;  # avoid overlapping labels
  for (my $i=0; $i<@parts; $i++) {
    my $x1     = $parts[$i][0]->left;
    my $x2     = $parts[$i][-1]->right;

    my $string = $self->part_label($i,scalar @parts);
    my $x    = $left + $x1 + ($x2 - $x1 - $width*length($string))/2;
    my $w    = $width * length($string);
    next if defined $last_x && $self->flip ?  $x + $w > $last_x : $x < $last_x;
    $gd->string($font,
		$x,$y,
		$string,
		$color);
    $last_x = $x + ($self->flip ? 0 : $w);
  }
}

sub part_label {
  my $self = shift;
  my ($part,$total)  = @_;

  local $self->{partno} = $self->feature->strand < 0 ? $total - $part -1 : $part;
  my $label = $self->option('part_labels');
  return unless defined $label;
  return "1"   if $label eq '1 ';
  return $label unless $label eq '1';
  return $self->{partno}+1;
}

sub dna_fits {
  my $self = shift;

  my $pixels_per_base = $self->scale;
  my $font            = $self->font;
  my $font_width      = $font->width;

  return $pixels_per_base >= $font_width;
}

sub protein_fits {
  my $self = shift;
  my $font               = $self->font;

  # return unless $font->height <= $self->height;

  my $font_width         = $font->width;
  my $pixels_per_residue = $self->scale * 3;

  return $pixels_per_residue >= $font_width;
}

sub arrowhead {
  my $self = shift;
  my $image = shift;
  my ($x,$y,$height,$orientation) = @_;

  my $fg = $self->set_pen;
  my $style = $self->option('arrowstyle') || 'regular';

  if ($style eq 'filled') {
    my $poly_pkg = $self->polygon_package;
    my $poly = $poly_pkg->new();
    if ($orientation >= 0) {
      $poly->addPt($x-$height,$y-$height);
      $poly->addPt($x,$y);
      $poly->addPt($x-$height,$y+$height,$y);
    } else {
      $poly->addPt($x+$height,$y-$height);
      $poly->addPt($x,$y);
      $poly->addPt($x+$height,$y+$height,$y);
    }
    $image->filledPolygon($poly,$fg);
  }
  else {
    if ($orientation >= 0) {
      $image->line($x,$y,$x-$height,$y-$height,$fg);
      $image->line($x,$y,$x-$height,$y+$height,$fg);
    } else {
      $image->line($x,$y,$x+$height,$y-$height,$fg);
      $image->line($x,$y,$x+$height,$y+$height,$fg);
    }
  }
}

sub arrow {
  my $self  = shift;
  my $image = shift;
  my ($x1,$x2,$y) = @_;

  my $fg     = $self->set_pen;
  my $height = $self->height/3;

  $image->line($x1,$y,$x2,$y,$fg);
  $self->arrowhead($image,$x2,$y,$height,+1) if $x1 < $x2;
  $self->arrowhead($image,$x2,$y,$height,-1) if $x2 < $x1;
}

sub reversec {
  my $self = shift;
  my $dna  = shift;
  $dna =~ tr/gatcGATC/ctagCTAG/;
  $dna = reverse $dna;
  return $dna;
}

# this gets invoked if the user has requested that the protein translation
# gets drawn using the draw_translation option and protein_fits() returns
# true. It is a rather specialized function and possibly belongs somewhere else,
# but putting it here makes it possible for any feature to display its protein
# translation.
sub calculate_cds {
  my $self = shift;
  my @parts = $self->feature_has_subparts ? $self->parts : $self;

  my $codon_table = $self->option('codontable');
  $codon_table    = 1 unless defined $codon_table;
  require Bio::Tools::CodonTable unless Bio::Tools::CodonTable->can('new');
  my $translate_table = Bio::Tools::CodonTable->new(-id=>$codon_table);

  for (my $i=0; $i < @parts; $i++) {
    my $part    = $parts[$i];
    my $feature = $part->feature;

    my $pos     = $feature->strand >= 0 ? $feature->start : $feature->end;
    my $phase   = eval {$feature->phase};
    next unless defined $phase;
    my $seq     = $feature->seq;
    next unless defined $seq;

    my $strand          = $feature->strand;
    my ($frame,$offset) = frame_and_offset($pos,
					   $strand,
					   -$phase);
    $strand *= -1 if $self->{flip};
    $part->{cds_frame}     = $frame;
    $part->{cds_offset}    = $offset;

    # do in silico splicing in order to find the codon that
    # arises from the splice
    my $protein = $seq->translate(undef,undef,$phase,$codon_table)->seq;
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
      my $left_of_splice  = substr($self->get_seq($feature->seq),     -$next_phase, $next_phase);
      my $right_of_splice = substr($self->get_seq($next_feature->seq),0           , 3-$next_phase);
      $splice_codon = $left_of_splice . $right_of_splice;
      length $splice_codon == 3                      or last BLOCK;
      my $amino_acid = $translate_table->translate($splice_codon);
      $part->{cds_splice_residue} = $amino_acid;
    }
  }
}

# hack around changed feature API
sub get_seq {
  my $self = shift;
  my $seq = shift;
  return $seq if ref $seq && $seq->can('translate');
  require Bio::PrimarySeq unless Bio::PrimarySeq->can('new');
  return Bio::PrimarySeq->new(-seq=>$seq);
}

1;

=head1 NAME

Bio::Graphics::Glyph::generic - The "generic" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This is identical to the "box" glyph except that it will draw the
subparts of features that contain subfeatures.  The subparts are not
connected -- use the "segments" glyph for that.  "Generic" is the
default glyph used when not otherwise specified.

=head2 METHODS

This module overrides the maxdepth() method to return 0 unless the
-maxdepth option is provided explicitly. This means that any module
that inherits from generic will need to override maxdepth() again in
order to draw subfeatures. In general, those implementing
multi-segmented feature glyphs should inherit from
Bio::Graphics::Glyph::segments, which allows for one level of descent.

In addition, the following new methods are implemented:

=over 4

=item labelfont(), descfont(), labelwidth(), descriptionwidth()

Return the font, width for the label or description.

=item label()

Return the glyph label text (printed above the glyph).

=item description()

Return the glyph description text (printed below the glyph).

=item draw_translation()

Draw the protein translation of the feature (assumes that the feature is attached to a DNA sequence).

=item draw_sequence()

Draw the sequence of the feature (either DNA or protein).

=back

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

  -font         Default font                   gdSmallFont

  -label_font   Font used for label	       gdSmallFont

  -desc_font    Font used for description      gdSmallFont

  -connector    Connector type                 0 (false)

  -connector_color
                Connector color                black

  -pad_top      Top padding                    0

  -pad_bottom   Bottom padding                 0

  -label        Whether to draw a label	       0 (false)

  -label_position Where to draw the label      "top" (default) or "left"

  -description  Whether to draw a description  0 (false)

  -strand_arrow Whether to indicate            0 (false)
                 strandedness

  -hilite       Highlight color                undef (no color)

  -draw_dna     If true, draw the dna residues        0 (false)
                 when magnification level
                 allows.

  -canonical_strand If true, draw the dna residues        0 (false)
                 as they appear on the plus strand
                 even if the feature is on the minus
                 strand.

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
L<Bio::Graphics::Glyph::xyplot>,
L<Bio::DB::GFF>,
L<Bio::SeqI>,
L<Bio::SeqFeatureI>,
L<Bio::Das>,
L<GD>

=head1 AUTHOR
p
Allen Day E<lt>day@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
