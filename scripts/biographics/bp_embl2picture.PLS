#!/usr/bin/perl

# file: embl2picture.pl
# This is code example 6 in the Graphics-HOWTO
# Author: Lincoln Stein

use strict;
use lib "$ENV{HOME}/projects/bioperl-live";
use Bio::Graphics;
use Bio::SeqIO;

=head1 NAME

bp_embl2picture.pl - Render a Genbank/EMBL sequence file graphically as a png image

=head1 SYNOPSIS

  % bp_embl2picture.pl factor7.embl | display -

=head1 DESCRIPTION

Render a GenBank/EMBL entry into drawable form.  Return as a GIF or
PNG image on standard output.

The file must be in embl, genbank, or another SeqIO- recognized format.
Only the first entry will be rendered.

=head1 SEE ALSO

L<Bio::Graphics>, the BioGraphics HOWTO.

=head1 AUTHOR

Lincoln Stein, lstein@cshl.edu

Copyright (c) 2004 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

use constant USAGE =><<END;
Usage: $0 <file>
   Render a GenBank/EMBL entry into drawable form.
   Return as a GIF or PNG image on standard output.

   File must be in embl, genbank, or another SeqIO-
   recognized format.  Only the first entry will be 
   rendered.

Example to try:
   embl2picture.pl factor7.embl | display -

END

my $file = shift                       or die USAGE;
my $io = Bio::SeqIO->new(-file=>$file) or die USAGE;
my $seq = $io->next_seq                or die USAGE;
my $wholeseq = Bio::SeqFeature::Generic->new(-start=>1,-end=>$seq->length,
					     -display_name=>$seq->display_name);

my @features = $seq->all_SeqFeatures;

# sort features by their primary tags
my %sorted_features;
for my $f (@features) {
  my $tag = $f->primary_tag;
  push @{$sorted_features{$tag}},$f;
}

my $panel = Bio::Graphics::Panel->new(
				      -length    => $seq->length,
				      -key_style => 'between',
				      -width     => 800,
				      -pad_left  => 10,
				      -pad_right => 10,
				      );
$panel->add_track($wholeseq,
		  -glyph => 'arrow',
		  -bump => 0,
		  -double=>1,
		  -tick => 2);

$panel->add_track($wholeseq,
		  -glyph  => 'generic',
		  -bgcolor => 'blue',
		  -label  => 1,
		 );

# special cases
if ($sorted_features{CDS}) {
  $panel->add_track($sorted_features{CDS},
		    -glyph      => 'transcript2',
		    -bgcolor    => 'orange',
		    -fgcolor    => 'black',
		    -font2color => 'red',
		    -key        => 'CDS',
		    -bump       =>  +1,
		    -height     =>  12,
		    -label      => \&gene_label,
		    -description=> \&gene_description,
		   );
  delete $sorted_features{'CDS'};
}

if ($sorted_features{tRNA}) {
  $panel->add_track($sorted_features{tRNA},
		    -glyph     =>  'transcript2',
		    -bgcolor   =>  'red',
		    -fgcolor   =>  'black',
		    -font2color => 'red',
		    -key       => 'tRNAs',
		    -bump      =>  +1,
		    -height    =>  12,
		    -label     => \&gene_label,
		   );
  delete $sorted_features{tRNA};
}

# general case
my @colors = qw(cyan orange blue purple green chartreuse magenta yellow aqua);
my $idx    = 0;
for my $tag (sort keys %sorted_features) {
  my $features = $sorted_features{$tag};
  $panel->add_track($features,
		    -glyph    =>  'generic',
		    -bgcolor  =>  $colors[$idx++ % @colors],
		    -fgcolor  => 'black',
		    -font2color => 'red',
		    -key      => "${tag}s",
		    -bump     => +1,
		    -height   => 8,
		    -description => \&generic_description
		   );
}

print $panel->png;
exit 0;

sub gene_label {
  my $feature = shift;
  my @notes;
  foreach (qw(product gene)) {
    next unless $feature->has_tag($_);
    @notes = $feature->each_tag_value($_);
    last;
  }
  $notes[0];
}

sub gene_description {
  my $feature = shift;
  my @notes;
  foreach (qw(note)) {
    next unless $feature->has_tag($_);
    @notes = $feature->each_tag_value($_);
    last;
  }
  return unless @notes;
  substr($notes[0],30) = '...' if length $notes[0] > 30;
  $notes[0];
}

sub generic_description {
  my $feature = shift;
  my $description;
  foreach ($feature->all_tags) {
    my @values = $feature->each_tag_value($_);
    $description .= $_ eq 'note' ? "@values" : "$_=@values; ";
  }
  $description =~ s/; $//; # get rid of last
  $description;
}
