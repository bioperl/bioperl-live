package Bio::Graphics::Glyph::segments;
#$Id$

use strict;
use Bio::Location::Simple;

use constant RAGGED_START_FUZZ => 25;  # will show ragged ends of alignments
                                       # up to this many bp.

use constant DEBUG => 0;

# These are just offsets into an array data structure
use constant TARGET    => 0;
use constant SRC_START => 1;
use constant SRC_END   => 2;
use constant TGT_START => 3;
use constant TGT_END   => 4;

use base qw(Bio::Graphics::Glyph::segmented_keyglyph Bio::Graphics::Glyph::generic);

my %complement = (g=>'c',a=>'t',t=>'a',c=>'g',n=>'n',
		  G=>'C',A=>'T',T=>'A',C=>'G',N=>'N');

sub pad_left {
  my $self = shift;
  return $self->SUPER::pad_left unless $self->level > 0;
  my $ragged = $self->option('ragged_start') 
    ? RAGGED_START_FUZZ 
    : $self->option('ragged_extra');

  return $self->SUPER::pad_left 
    unless $self->draw_target && $ragged && $self->dna_fits;
  my $target = eval {$self->feature->hit} or return $self->SUPER::pad_left;

  return $self->SUPER::pad_left unless $target->start<$target->end && $target->start < $ragged;
  return ($target->start-1) * $self->scale;
}

sub pad_right {
  my $self = shift;
  return $self->SUPER::pad_right unless $self->level > 0;
  my $ragged = $self->option('ragged_start') 
    ? RAGGED_START_FUZZ 
    : $self->option('ragged_extra');
  return $self->SUPER::pad_right 
    unless $self->draw_target && $ragged && $self->dna_fits;
  my $target = eval {$self->feature->hit} or return $self->SUPER::pad_right;
  return $self->SUPER::pad_right unless $target->end < $target->start && $target->start < $ragged;
  return ($target->end-1) * $self->scale;
}

sub draw_target {
  my $self = shift;
  return if $self->option('draw_dna');
  return $self->option('draw_target');
}

sub draw_protein_target {
  my $self = shift;
  return if $self->option('draw_protein');
  return $self->option('draw_protein_target');
  return $self->option('draw_target');
}

sub height {
  my $self = shift;
  my $height = $self->SUPER::height;
  return $height unless $self->draw_target || $self->draw_protein_target;
  if ($self->draw_target) {
    return $height unless $self->dna_fits;
  }
  if ($self->draw_protein_target) {
    return $height unless $self->protein_fits;
  }
  my $fontheight = $self->font->height;
  return $fontheight if $fontheight > $height;
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

sub maxdepth {
  my $self = shift;
  my $md   = $self->Bio::Graphics::Glyph::maxdepth;
  return $md if defined $md;
  return 1;
}

# this was willfully confusing
#sub fontcolor {
#  my $self = shift;
#  return $self->SUPER::fontcolor unless $self->draw_target;# || $self->option('draw_dna');
#  return $self->SUPER::fontcolor unless $self->dna_fits;
#  return $self->bgcolor;
#}

sub draw {
  my $self = shift;

  my $draw_target         = $self->draw_target;
  return $self->SUPER::draw(@_) unless $draw_target;
  return $self->SUPER::draw(@_) unless $self->dna_fits;

  $self->draw_label(@_)       if $self->option('label');
  $self->draw_description(@_) if $self->option('description');
  $self->draw_part_labels(@_) if $self->option('part_labels');

  my $drew_sequence;

  if ($draw_target) {
    return $self->SUPER::draw(@_) unless eval {$self->feature->hit->seq};
    $drew_sequence = $self->draw_multiple_alignment(@_);
  }

  my ($gd,$x,$y) = @_;
  $y  += $self->top + $self->pad_top if $drew_sequence;  # something is wrong - this is a hack/workaround
  my $connector     =  $self->connector;
  $self->draw_connectors($gd,$x,$y)
    if $connector && $connector ne 'none' && $self->level == 0;

}

sub draw_multiple_alignment {
  my $self = shift;
  my $gd   = shift;
  my ($left,$top,$partno,$total_parts) = @_;

  my $flipped              = $self->flip;
  my $ragged_extra         = $self->option('ragged_start') 
                               ? RAGGED_START_FUZZ : $self->option('ragged_extra');
  my $true_target          = $self->option('true_target');
  my $show_mismatch        = $self->option('show_mismatch');
  my $do_realign           = $self->option('realign');

  my $pixels_per_base      = $self->scale;
  my $feature              = $self->feature;
  my $panel                = $self->panel;
  my ($abs_start,$abs_end)     = ($feature->start,$feature->end);
  my ($tgt_start,$tgt_end)     = ($feature->hit->start,$feature->hit->end);
  my ($panel_start,$panel_end) = ($self->panel->start,$self->panel->end);
  my $strand               = $feature->strand;
  my $panel_left           = $self->panel->left;
  my $panel_right          = $self->panel->right;
  my $drew_sequence;

  if ($tgt_start > $tgt_end) { #correct for data problems
    $strand    = -1;
    ($tgt_start,$tgt_end) = ($tgt_end,$tgt_start);
  }

  warn "TGT_START..TGT_END = $tgt_start..$tgt_end" if DEBUG;

  my ($bl,$bt,$br,$bb)     = $self->bounds($left,$top);
  $top = $bt;

  for my $p ($self->parts) {
    my @bounds = $p->bounds($left,$top);
    $self->filled_box($gd,@bounds,$self->bgcolor,$self->bgcolor);
  }

  my @s                     = $self->_subfeat($feature);

  # FIX ME
  # workaround for features in which top level feature does not have a hit but
  # subfeatures do. There is total breakage of encapsulation here because sometimes
  # a chado alignment places the aligned segment in the top-level feature, and sometimes
  # in the child feature.
  unless (@s || $feature->isa('Bio::DB::GFF::Feature')) {
    @s = ($feature);
  }

  my $can_realign = $do_realign && eval { require Bio::Graphics::Browser::Realign; 1 };

  my (@segments,%strands);
  for my $s (@s) {
    my $target = $s->hit;
    my ($src_start,$src_end) = ($s->start,$s->end);
    next unless $src_start <= $panel_end && $src_end >= $panel_start;

    my ($tgt_start,$tgt_end) = ($target->start,$target->end);

    my $strand_bug;
    unless (exists $strands{$target}) {
      my $strand = $feature->strand;
      if ($tgt_start > $tgt_end) { #correct for data problems
	$strand    = -1;
	($tgt_start,$tgt_end) = ($tgt_end,$tgt_start);
	$strand_bug++;
      }
      $strands{$target} = $strand;
    }

    # realign for internal gaps, if requested
    if ($can_realign) {
      warn   "Realigning [$target,$src_start,$src_end,$tgt_start,$tgt_end].\n" if DEBUG;
      my ($sdna,$tdna) = ($s->dna,$target->dna);
      my @result = $self->realign($sdna,$tdna);
      foreach (@result) {
	warn "=========> [$target,@$_]\n" if DEBUG;
	my $a = $strands{$target} >= 0 ? [$target,$_->[0]+$src_start,$_->[1]+$src_start,$_->[2]+$tgt_start,$_->[3]+$tgt_start]
	                               : [$target,$src_end-$_->[1],$src_end-$_->[0],$_->[2]+$tgt_start,$_->[3]+$tgt_start];
	warn "[$target,$_->[0]+$src_start,$_->[1]+$src_start,$tgt_end-$_->[3],$tgt_end-$_->[2]]" if DEBUG;
	warn "=========> [@$a]\n" if DEBUG;
	warn substr($sdna,     $_->[0],$_->[1]-$_->[0]+1),"\n" if DEBUG;
	warn substr($tdna,$_->[2],$_->[3]-$_->[2]+1),"\n"      if DEBUG;
	push @segments,$a;
      }
    }
    else {
      push @segments,[$target,$src_start,$src_end,$tgt_start,$tgt_end];
    }
  }

  # get 'em in the right order so that we don't have to worry about
  # where the beginning and end are.
  @segments = sort {$a->[TGT_START]<=>$b->[TGT_START]} @segments;

  # adjust for ragged (nonaligned) ends
  my ($offset_left,$offset_right) = (0,0);
  if ($ragged_extra && $ragged_extra > 0) {

    # add a little rag to the left end
    $offset_left = $segments[0]->[TGT_START] > $ragged_extra ? $ragged_extra : $segments[0]->[TGT_START]-1;
    if ($strand >= 0) {
      $offset_left     = $segments[0]->[SRC_START]-1 if $segments[0]->[SRC_START] - $offset_left < 1;
      $abs_start                -= $offset_left;
      $tgt_start                -= $offset_left;
      $segments[0]->[SRC_START] -= $offset_left;
      $segments[0]->[TGT_START] -= $offset_left;
    } else {
      $abs_end                  += $offset_left;
      $tgt_start                -= $offset_left;
      $segments[0]->[SRC_END]   += $offset_left;
      $segments[0]->[TGT_START] -= $offset_left;
    }

    # add a little rag to the right end - this is complicated because
    # we don't know what the length of the underlying dna is, so we
    # use the subfeat method to find out
    my $current_end     = $segments[-1]->[TGT_END];
    $offset_right          = length $segments[-1]->[TARGET]->subseq($current_end+1,$current_end+$ragged_extra)->seq;
    if ($strand >= 0) {
      $abs_end                 += $offset_right;
      $tgt_end                 += $offset_left;
      $segments[-1]->[TGT_END] += $offset_right;
      $segments[-1]->[SRC_END] += $offset_right;
    } else {
      $abs_start                 -= $offset_right;
      $tgt_end                   += $offset_left;
      $segments[-1]->[TGT_END]   += $offset_right;
      $segments[-1]->[SRC_START] -= $offset_right;
    }
  }

  # get the DNAs now - a little complicated by the necessity of using
  # the subseq() method
  my $ref_dna = $feature->subseq(1-$offset_left,$feature->length+$offset_right)->seq;
  my $tgt_dna = $feature->hit->subseq(1-$offset_left,$feature->length+$offset_right)->seq;

  # work around changes in the API
  $ref_dna    = $ref_dna->seq if ref $ref_dna and $ref_dna->can('seq');
  $tgt_dna    = $tgt_dna->seq if ref $tgt_dna and $tgt_dna->can('seq');

  $ref_dna    = lc $ref_dna;
  $tgt_dna    = lc $tgt_dna;

  # sanity check.  Let's see if they look like they're lining up
  warn "$feature dna sanity check:\n$ref_dna\n$tgt_dna\n" if DEBUG;

  # now we're all lined up, and we're going to adjust everything to fall within the bounds
  # of the left and right panel coordinates
  my %clip;
  for my $seg (@segments) {

    my $target = $seg->[TARGET];
    warn "preclip [@$seg]\n" if DEBUG;

    # left clipping
    if ( (my $delta = $seg->[SRC_START] - $panel_start) < 0 ) {
      warn "clip left delta = $delta" if DEBUG;
      $seg->[SRC_START] = $panel_start;
      if ($strand >= 0) {
	$seg->[TGT_START] -= $delta;
      }
    }

    # right clipping
    if ( (my $delta = $panel_end - $seg->[SRC_END]) < 0) {
      warn "clip right delta = $delta" if DEBUG;
      $seg->[SRC_END] = $panel_end;
      if ($strand < 0) {
	$seg->[TGT_START] -= $delta;
      }
    }

    my $length = $seg->[SRC_END]-$seg->[SRC_START]+1;
    $seg->[TGT_END] = $seg->[TGT_START]+$length-1;

    warn "Clipping gives [@$seg], tgt_start = $tgt_start\n" if DEBUG;
  }

  # remove segments that got clipped out of existence
  @segments = grep { $_->[SRC_START]<=$_->[SRC_END] } @segments;

  # relativize coordinates
  if ($strand < 0) {
    $ref_dna = $self->reversec($ref_dna);
    $tgt_dna = $self->reversec($tgt_dna);
  }

  for my $seg (@segments) {
    $seg->[SRC_START] -= $abs_start - 1;
    $seg->[SRC_END]   -= $abs_start - 1;
    $seg->[TGT_START] -= $tgt_start - 1;
    $seg->[TGT_END]   -= $tgt_start - 1;

    warn "src segment = $seg->[SRC_START]", "..",$seg->[SRC_END] if DEBUG;
    warn "tgt segment = $seg->[TGT_START]", "..",$seg->[TGT_END] if DEBUG;
    if ($strand < 0) {
      ($seg->[TGT_START],$seg->[TGT_END]) = (length($tgt_dna)-$seg->[TGT_END]+1,length($tgt_dna)-$seg->[TGT_START]+1);
    }
    if (DEBUG) {
      warn "$feature: relativized coordinates = [@$seg]\n";
      warn $self->_subsequence($ref_dna,$seg->[SRC_START],$seg->[SRC_END]),"\n";
      warn $self->_subsequence($tgt_dna,$seg->[TGT_START],$seg->[TGT_END]),"\n";
    }
  }

  # draw
  my $color = $self->fgcolor;
  my $font  = $self->font;
  my $lineheight = $font->height;
  my $fontwidth  = $font->width;

  my $mismatch = $self->factory->translate_color($self->option('mismatch_color') || 'lightgrey');
  my $grey     = $self->factory->translate_color('gray');

  my $base2pixel = 
    $self->flip ?
      sub {
	my ($src,$tgt) = @_;
	my $a = $fontwidth + ($abs_start + $src-$panel_start-1 + $tgt) * $pixels_per_base - 1;    
	$panel_right - $a;
      }
      : sub {
	my ($src,$tgt) = @_;
	$fontwidth/2 + $left + ($abs_start + $src-$panel_start-1 + $tgt) * $pixels_per_base - 1;    
      };

  my ($tgt_last_end,$src_last_end,$leftmost,$rightmost);
  for my $seg (sort {$a->[SRC_START]<=>$b->[SRC_START]} @segments) {
    my $y = $top-1;

    for (my $i=0; $i<$seg->[SRC_END]-$seg->[SRC_START]+1; $i++) {

      my $src_base = $self->_subsequence($ref_dna,$seg->[SRC_START]+$i,$seg->[SRC_START]+$i);
      my $tgt_base = $self->_subsequence($tgt_dna,$seg->[TGT_START]+$i,$seg->[TGT_START]+$i);
#      warn $seg->[TGT_START]+$i,' ',$seg->[TGT_START]+$i if DEBUG;
      my $x = $base2pixel->($seg->[SRC_START],$i);
      $leftmost = $x if !defined $leftmost  || $leftmost  > $x;
      $rightmost= $x if !defined $rightmost || $rightmost < $x;

      next unless $tgt_base && $x >= $panel_left && $x <= $panel_right;

      $self->filled_box($gd,$x-$pixels_per_base/2+2,$y+1,$x+$pixels_per_base/2+1,$y+$lineheight,$mismatch,$mismatch)
	if $show_mismatch && $tgt_base && $src_base ne $tgt_base && $tgt_base !~ /[nN]/;
      $tgt_base = $complement{$tgt_base} if $true_target && $strand < 0;
      $gd->char($font,$x,$y,$tgt_base,$tgt_base =~ /[nN]/ ? $grey : $color);

      $drew_sequence++;
    }

    # indicate the presence of insertions in the target
    if (defined $tgt_last_end) {
      my $delta     = $seg->[TGT_START] - $tgt_last_end;
      my $src_delta = $seg->[SRC_START] - $src_last_end;
      if ($delta > 1 and $src_delta > 0) {  # an insertion in the target relative to the source
	my $gap_left  = $fontwidth + $base2pixel->($src_last_end,0);
	my $gap_right = $base2pixel->($seg->[SRC_START],0);
	($gap_left,$gap_right) = ($gap_right+$fontwidth,$gap_left-$fontwidth) if $self->flip;
	warn "delta=$delta, gap_left=$gap_left, gap_right=$gap_right" if DEBUG;

	if ($delta == $src_delta) {
	  $gap_left  += $pixels_per_base/2-2;
	  $gap_right -= $pixels_per_base/2-2;
	}

	next if $gap_left <= $panel_left || $gap_right >= $panel_right;

	$self->filled_box($gd,$gap_left,$y+1,
			      $gap_right-2,$y+$lineheight,$mismatch,$mismatch) if
				$show_mismatch && $gap_left >= $panel_left && $gap_right <= $panel_right;


	my $gap_distance             = $gap_right - $gap_left + 1;
	my $pixels_per_inserted_base = $gap_distance/($delta-1);

 	if ($pixels_per_inserted_base >= $fontwidth) {  # Squeeze the insertion in
 	  for (my $i = 0; $i<$delta-1; $i++) {
 	    my $x = $gap_left + ($pixels_per_inserted_base-$fontwidth)/2 + $pixels_per_inserted_base * $i;
 	    my $bp = $self->_subsequence($tgt_dna,$tgt_last_end+$i+1,$tgt_last_end+$i+1);
	    next if $x < $panel_left;
 	    $gd->char($font,$x,$y,$bp,$color);
 	  }
	}
	# stick in a blob
	$self->_draw_insertion_point($gd,$gap_left,$gap_right,$y,$y+$lineheight,$mismatch) if $delta > 2;
      }
      # deal with gaps in the alignment
      elsif ( (my $delta = $seg->[SRC_START] - $src_last_end) > 1) {
	for (my $i=0;$i<$delta-1;$i++) {
	  my $x = $base2pixel->($src_last_end,$i+1);
	  next if $x > $panel_right;
	  $self->filled_box($gd,$x-$pixels_per_base/2+2,$y,$x+$pixels_per_base/2+1,$y+$lineheight,$mismatch,$mismatch)
	    if $show_mismatch;
	  $gd->char($font,$x,$y,'-',$color);
	}
	
      }

    }

    $tgt_last_end  = $seg->[TGT_END];
    $src_last_end  = $seg->[SRC_END];
  }

  # additional fixup -- insert dashes at beginnings and ends of the segments, if they don't span the full
  # alignment for some reason - THIS SHOULD NOT BE NECESSARY AND INDICATES THAT THIS WHOLE METHOD NEEDS
  # TO BE REWRITTEN!
  if (defined $leftmost && $leftmost-$bl > $pixels_per_base) {
    $gd->char($font,$_,$top-1,'-',$color) for map {$bl+$_*$pixels_per_base} 0..($leftmost-$bl)/$pixels_per_base-1;
  }
  if (defined $rightmost && $br-$rightmost > $pixels_per_base) {
    $gd->char($font,$_,$top-1,'-',$color) for map {$rightmost+$_*$pixels_per_base} (0..($br-$rightmost)/$pixels_per_base);
  }

  return $drew_sequence;
}

sub _subsequence {
  my $self = shift;
  my ($seq,$start,$end,$strand) = @_;
  my $sub;
  if ((defined $strand && $strand < 0)) {
    my $piece = substr($seq,length($seq)-$end,$end-$start+1);
    $sub = $self->reversec($piece);
  } else {
    $sub = substr($seq,$start-1,$end-$start+1);
  }
  return $self->flip ? $complement{$sub} : $sub;
}

sub realign {
  my $self = shift;
  my ($src,$tgt) = @_;
  return Bio::Graphics::Browser::Realign::align_segs($src,$tgt);
}

# Override _subfeat() method to make it appear that a top-level feature that
# has no subfeatures appears as a feature that has a single subfeature.
# Otherwise at high mags gaps will be drawn as components rather than
# as connectors.  Because of differing representations of split features
# in Bio::DB::GFF::Feature and Bio::SeqFeature::Generic, there is
# some breakage of encapsulation here.
sub _subfeat {
  my $self    = shift;
  my $feature = shift;
  my @subfeat  = $self->SUPER::_subfeat($feature);
  return @subfeat if @subfeat;
  if ($self->level == 0 && !@subfeat && !$self->feature_has_subparts) {
    return $self->feature;
  } else {
    return;
  }
}

# draw the classic "i-beam" icon to indicate that an insertion fits between
# two bases
# sub _draw_insertion_point {
#   my $self = shift;
#   my ($gd,$x,$y,$color) = @_;
#   my $top    = $y;
#   $x--;
#   my $bottom = $y + $self->font->height - 4;
#   $gd->line($x,$top+2, $x,$bottom-2,$color);
#   $gd->setPixel($x+1,  $top+1,$color);
#   $gd->setPixel($x+$_, $top,$color) for (2..3);
#   $gd->setPixel($x-1,  $top+1,$color);
#   $gd->setPixel($x-$_, $top,$color) for (2..3);

#   $gd->setPixel($x+1,  $bottom-1,$color);
#   $gd->setPixel($x+$_, $bottom,  $color) for (2..3);
#   $gd->setPixel($x-1,  $bottom-1,$color);
#   $gd->setPixel($x-$_, $bottom,  $color) for (2..3);
# }

# don't like that -- try drawing carets
sub _draw_insertion_point {
   my $self = shift;
   my ($gd,$left,$right,$top,$bottom,$color) = @_;

   my $poly = GD::Polygon->new();
   $poly->addPt($left-3,$top+1);
   $poly->addPt($right+2,$top+1);
   $poly->addPt(($left+$right)/2-1,$top+3);
   $gd->filledPolygon($poly,$color);

   $poly = GD::Polygon->new();
   $poly->addPt($left-3,$bottom);
   $poly->addPt($right+2,$bottom);
   $poly->addPt(($left+$right)/2-1,$bottom-2);
   $gd->filledPolygon($poly,$color);
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

=head2 METHODS

This module overrides the maxdepth() method to return 1 unless
explicitly specified by the -maxdepth option. This means that modules
inheriting from segments will only be presented with one level of
subfeatures. Override the maxdepth() method to get more levels.

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

In addition, the following glyph-specific options are recognized:

  -draw_dna     If true, draw the dna residues        0 (false)
                 when magnification level
                 allows.

  -draw_target  If true, draw the dna residues        0 (false)
                 of the TARGET sequence when
                 magnification level allows.
                 See "Displaying Alignments".

  -draw_protein_target  If true, draw the protein residues        0 (false)
                 of the TARGET sequence when
                 magnification level allows.
                 See "Displaying Alignments".

  -ragged_extra When combined with -draw_target,      0 (false)
                draw extra bases beyond the end
                of the alignment. The value is
                the maximum number of extra
                bases.
                See "Displaying Alignments".

  -ragged_start  Deprecated option.  Use
                 -ragged_extra instead

  -show_mismatch When combined with -draw_target,     0 (false)
                 highlights mismatched bases in
                 the mismatch color.  
                 See "Displaying Alignments".

  -mismatch_color The mismatch color to use           'lightgrey'

  -true_target   Show the target DNA in its native    0 (false)
                 (plus strand) orientation, even if
                 the alignment is to the minus strand.
                 See "Displaying Alignments".

  -realign       Attempt to realign sequences at      0 (false)
                 high mag to account for indels.
                 See "Displaying Alignments".

If the -draw_dna flag is set to a true value, then when the
magnification is high enough, the underlying DNA sequence will be
shown.  This option is mutually exclusive with -draw_target. See
Bio::Graphics::Glyph::generic for more details.

The -draw_target, -ragged_extra, and -show_mismatch options only work
with seqfeatures that implement the hit() method
(Bio::SeqFeature::SimilarityPair). -draw_target will cause the DNA of
the hit sequence to be displayed when the magnification is high enough
to allow individual bases to be drawn. The -ragged_extra option will
cause the alignment to be extended at the extreme ends by the
indicated number of bases, and is useful for looking for polyAs and
cloning sites at the ends of ESTs and cDNAs. -show_mismatch will cause
mismatched bases to be highlighted in with the color indicated by
-mismatch_color (default lightgray).

At high magnifications, minus strand matches will automatically be
shown as their reverse complement (so that the match has the same
sequence as the plus strand of the source dna).  If you prefer to see
the actual sequence of the target as it appears on the minus strand,
then set -true_target to true.

Note that -true_target has the opposite meaning from
-canonical_strand, which is used in conjunction with -draw_dna to draw
minus strand features as if they appear on the plus strand.

=head2 Displaying Alignments

When the B<-draw_target> option is true, this glyph can be used to
display nucleotide alignments such as BLAST, FASTA or BLAT
similarities.  At high magnification, this glyph will attempt to show
how the sequence of the source (query) DNA matches the sequence of the
target (the hit).  For this to work, the feature must implement the
hit() method, and both the source and the target DNA must be
available.  If you pass the glyph a series of
Bio::SeqFeature::SimilarityPair objects, then these criteria will be
satisified.

Without additional help, this glyph cannot display gapped alignments
correctly.  To display gapped alignments, you can use the
Bio::Graphics::Brower::Realign module, which is part of the Generic
Genome Browser package (http://www.gmod.org).  If you wish to install
the Realign module and not the rest of the package, here is the
recipe:

  cd Generic-Genome-Browser-1.XX
  perl Makefile.PL DO_XS=1
  make
  make install_site

If possible, build the gbrowse package with the DO_XS=1 option.  This
compiles a C-based DP algorithm that both gbrowse and gbrowse_details
will use if they can.  If DO_XS is not set, then the scripts will use
a Perl-based version of the algorithm that is 10-100 times slower.

The display of alignments can be tweaked using the -ragged_extra,
-show_mismatch, -true_target, and -realign options.  See the options
section for further details.

There is also a B<-draw_protein_target> option, which is designed for
protein to nucleotide alignments. It draws the target sequence every
third base pair and is supposed to align correctly with the forward
and reverse translation glyphs. This option is experimental at the
moment, and may not work correctly, to use with care.

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
