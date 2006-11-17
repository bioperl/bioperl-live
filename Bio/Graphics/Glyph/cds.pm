package Bio::Graphics::Glyph::cds;

use strict;
use Bio::Graphics::Glyph::segments;
use Bio::Graphics::Util qw(frame_and_offset);
use Bio::Tools::CodonTable;
use base qw(Bio::Graphics::Glyph::segmented_keyglyph Bio::Graphics::Glyph::translation);

my %default_colors = qw(
			frame0f  cornflowerblue
			frame1f  blue
			frame2f  darkblue
			frame0r  magenta
			frame1r  red
			frame2r  darkred
		       );

my %swap_phase = ( 0  => 0,
		   1  => 2,
		   2  => 1,
		   '' => 0);

sub connector   { 0 };
sub description {
  my $self = shift;
  return if $self->level;
  return $self->SUPER::description;
};

sub default_color {
  my ($self,$key) = @_;
  return $self->factory->translate_color($default_colors{$key});
}

sub sixframe {
  my $self = shift;
  return $self->{sixframe} if exists $self->{sixframe};
  my $sixframe = $self->option('sixframe');
  $sixframe    = $self->option('translation') eq '6frame' unless defined $sixframe;
  return $self->{sixframe} = $sixframe;
}

sub maxdepth { 1 };

sub require_subparts {
  my $self = shift;
  my $rs   = $self->option('require_subparts');
  $rs      = $self->feature->type eq 'coding' if !defined $rs;  # shortcut for the "coding" aggregator
  $rs;
}

sub ignore_undef_phase {
  shift->option('ignore_empty_phase');
}

sub ignore_non_cds {
  shift->option('cds_only');
}

sub phase_style {
  shift->option('phase_style') || '012';
}

# figure out (in advance) the color of each component
sub draw {
  my $self = shift;
  my ($gd,$left,$top) = @_;

  my @parts = $self->parts;
  @parts    = $self if !@parts && $self->level == 0 && !$self->require_subparts;

  my $fits = $self->protein_fits;
  my $strand = $self->feature->strand || 1;

  # draw the staff (musically speaking)
  if ($self->level == 0) {
    my ($x1,$y1,$x2,$y2) = $self->bounds($left,$top);
    my $line_count = $self->sixframe ? 6 : 3;
    my $height = ($y2-$y1)/$line_count;
    my $grid  = $self->gridcolor;
    for (0..$line_count-1) {
      my $offset = $y1+$height*$_+1;
      $gd->line($x1,$offset,$x2,$offset,$grid);
      # with three-frame translation, the position of the arrows changes depending on
      # the strand of the feature. With six-frame translation, we draw the first three
      # staff lines with an arrow to the right, and the second three to the left
      my $forward = ($line_count == 6) ? ($_ < 3) : ($strand > 0);
      if ($forward) {
	$gd->line($x2,$offset,$x2-2,$offset-2,$grid);
	$gd->line($x2,$offset,$x2-2,$offset+2,$grid);
      } else {
	$gd->line($x1,$offset,$x1+2,$offset-2,$grid);
	$gd->line($x1,$offset,$x1+2,$offset+2,$grid);
      }
    }
  }

  $self->{cds_part2color} ||= {};
  my $fill   = $self->bgcolor;

  # figure out the colors of each part
  # sort minus strand features backward
  @parts = map { $_->[0] }
  sort { $b->[1] <=> $a->[1] }
  map { [$_, $_->left ] } @parts if $strand < 0;

  my $codon_table = $self->option('codontable');
  $codon_table    = 1 unless defined $codon_table;
  my $translate_table = Bio::Tools::CodonTable->new(-id=>$codon_table);

  my $ignore_undef_phase = $self->ignore_undef_phase;
  my $ignore_non_cds     = $self->ignore_non_cds;
  my $broken_phase       = $self->phase_style eq '021';

  for (my $i=0; $i < @parts; $i++) {
    my $part    = $parts[$i];
    my $feature = $part->feature;

    my $type = $feature->method;
    next if ($self->option('sub_part') && $type ne $self->option('sub_part'));

    next if $ignore_non_cds && lc($type) ne 'cds';

    my $pos     = $feature->strand >= 0 ? $feature->start : $feature->end;
    my $phase   = $feature->can('phase') ? $feature->phase  # bioperl uses "frame" but this is incorrect usage
                 :$feature->can('frame') ? $feature->frame
                 :undef;
    next if $ignore_undef_phase && !defined($phase);
    $phase ||= 0;
    $phase = $swap_phase{$phase} if $broken_phase;
    my $strand  = $feature->strand;
    my ($frame,$offset) = frame_and_offset($pos,
					   $strand,
					   $phase);
    my $suffix = $strand < 0 ? 'r' : 'f';
    my $key = "frame$frame$suffix";
    $self->{cds_frame2color}{$key} ||= $self->color($key) || $self->default_color($key) || $fill;
    $part->{cds_partcolor} = $self->{cds_frame2color}{$key};
    $part->{cds_frame}     = $frame;
    $part->{cds_offset}    = $offset;

    if ($fits && $part->feature->seq) {

      # do in silico splicing in order to find the codon that
      # arises from the splice
      my $seq     = $self->get_seq($part->feature->seq);
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
	my $left_of_splice  = substr($self->get_seq($feature->seq),    -$next_phase, $next_phase);
	my $right_of_splice = substr($self->get_seq($next_feature->seq),0          , 3-$next_phase);
	$splice_codon = $left_of_splice . $right_of_splice;
	length $splice_codon == 3                      or last BLOCK;
	my $amino_acid = $translate_table->translate($splice_codon);
	$part->{cds_splice_residue} = $amino_acid;
      }
    }
  }

  $self->Bio::Graphics::Glyph::generic::draw($gd,$left,$top);
}


# draw the notes on the staff
sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->bounds(@_);

  my $color = $self->{cds_partcolor} or return;
  my $feature   = $self->feature;
  my $frame     = $self->{cds_frame};
  my $linecount = $self->sixframe ? 6 : 3;

  unless ($self->protein_fits) {
    my $height = ($y2-$y1)/$linecount;
    my $offset = $y1 + $height*$frame;
    $offset   += ($y2-$y1)/2 if $self->sixframe && $self->strand < 0;
    $offset   = $y1 + (($y2-$y1) - ($offset-$y1))-$height if $self->{flip}; # ugh. This works, but I don't know why
    $gd->filledRectangle($x1,$offset,$x2,$offset+2,$color);
    return;
  }

  # we get here if there's room to draw the primary sequence
  my $font  = $self->font;
  my $pixels_per_residue = $self->pixels_per_residue;
  my $strand = $feature->strand;
  my $y      = $y1-1;
  my $fontwidth = $font->width;

  $strand *= -1 if $self->{flip};

  # have to remap feature start and end into pixel coords in order to:
  # 1) correctly align the amino acids with the nucleotide seq
  # 2) correct for the phase offset
  my $start = $self->map_no_trunc($feature->start + $self->{cds_offset});
  my $stop  = $self->map_no_trunc($feature->end   + $self->{cds_offset});

  ($start,$stop) = ($stop,$start) if $stop < $start;  # why does this keep happening?
  #  ($start,$stop) = ($stop,$start) if $self->{flip};

  my @residues = split '',$self->{cds_translation};

  push @residues,$self->{cds_splice_residue} if $self->{cds_splice_residue};
  for (my $i=0;$i<@residues;$i++) {
    my $x = $strand > 0 ? $start + $i * $pixels_per_residue
                        : $stop  - $i * $pixels_per_residue;
    next unless ($x >= $x1 && $x <= $x2);
    $x -= $fontwidth + 1 if $self->{flip}; # align right when flipped
    $gd->char($font,$x+1,$y,$residues[$i],$color);
  }
}

sub make_key_feature {
  my $self = shift;
  my @gatc = qw(g a t c);
  my $offset = $self->panel->offset;
  my $scale = 1/$self->scale;  # base pairs/pixel
  my $start = $offset;
  my $stop  = $offset + 100 * $scale;
  my $seq   = join('',map{$gatc[rand 4]} (1..1500));
  my $feature =
    Bio::Graphics::Feature->new(-start=> $start,
				-end  => $stop,
				-seq  => $seq,
				-name => $self->option('key'),
				-strand=> +1,
			       );
  $feature->add_segment(Bio::Graphics::Feature->new(
						    -start=> $start,
						    -end => $start + ($stop - $start)/2,
						    -seq  => $seq,
						    -name => $self->option('key'),
						    -strand=> +1,
						   ),
			Bio::Graphics::Feature->new(
						    -start=> $start + ($stop - $start)/2+1,
						    -end => $stop,
						    -seq  => $seq,
						    -name => $self->option('key'),
						    -phase=> 1,
						    -strand=> +1,
						   ));
  $feature;
}

# never allow our components to bump
sub bump {
  my $self = shift;
  return $self->SUPER::bump(@_) if $self->all_callbacks;
  return 0;
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

  -hilite       Highlight color                undef (no color)

In addition, the cds glyph recognizes the following glyph-specific
options:

  Option      Description                      Default
  ------      -----------                      -------

  -frame0f    Color for first (+) frame        background color

  -frame1f    Color for second (+) frame       background color

  -frame2f    Color for third (+) frame        background color

  -frame0r    Color for first (-) frame        background color

  -frame1r    Color for second (-) frame       background color

  -frame2r    Color for third (-) frame        background color

  -gridcolor  Color for the "staff"            lightslategray

  -translation Number of lines of reading      3frame
               frames to show. One of
               "3frame", or "6frame".
               For 6frame, specify a height
               of at least 30 pixels.

  -sixframe   Draw a six-frame staff           0 (false; usually draws 3 frame)
              This value overrides
              -translation, which essentially
              does the same thing.

  -require_subparts
              Don't draw the reading frame 0   false
              unless it is a feature
              subpart.

  -sub_part   For objects with multiple	       undef
              subpart types, defines which
              is the CDS part.

  -codontable   Codon table to use             1 (see Bio::Tools::CodonTable)

  -phase_style  The way phase is to be
                interpreted. One of            "012"
                "012" or "021"
  -ignore_empty_phase                          false
              Only draw features that have
              their phase defined.

  -cds_only   Only draw features of type       false
              'CDS'

This glyph is more sensitive to the underlying data model than usual,
so there are a few additional options to use to help adapt the glyph
to different environments.

The -require_subparts option is suggested when rendering spliced
transcripts which contain multiple CDS subparts.  Otherwise, the glyph
will hickup when zoomed way down onto an intron between two CDSs (a
phantom reading frame will appear).  For unspliced sequences, do *not*
use -require_subparts.

The -phase_style controls how the value returned by the phase() or
frame() methods is to be interpreted. The official interpretation is
that the phase value indicates the offset into the feature at which
the reading frame starts -- e.g. a phase of "2" means the reading
frame starts after skipping two bases from the beginning of the
feature.  However, many GFF2 format feature files interpret this field
to mean the position reading frame of the first base of the feature --
e.g. a phase of "2" means that the reading frame starts after skipping
just one base from the beginning of the feature. Specify "012" to
interpret the phase field in the correct way, and "021" to interpret
the phase field in the legacy way. The default is "012."

Here is how the option names were chosen:

    * * *                  Base the reading frame starts on
    A B C A B C A B C...
    0 1 2                  PHASE REPRESENTED CORRECTLY
    0 2 1                  PHASE REPRESENTED IN THE LEGACY WAY

Set the -ignore_empty_phase option to true if you wish to skip
subfeatures that do not have a defined phase() or frame(). This is useful
if you are rendering exons that have both translated and untranslated
parts, and you wish to skip the untranslated parts.

Set the -cds_only option to true if you wish to draw the glyph only
for subfeatures of type 'CDS'. This is recommended.

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
