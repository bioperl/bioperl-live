#!/usr/bin/perl

#Thu Aug 19 22:29:23 EDT 2004
#Simon Ilyushchenko - demonstrating several new glyphs for gbrowse, part 1.

use strict;

use Bio::Graphics::Panel;

=head1 NAME

bp_glyphs1-demo.pl - First demo of Bio::Graphics glyphs

=head1 SYNOPSIS

  % bp_glyphs2-demo.pl | display -

=head1 DESCRIPTION

Generates a PNG image of some of the more esoteric Bio::Graphics glyphs.

=head1 SEE ALSO

L<Bio::Graphics>, the BioGraphics HOWTO.

=head1 AUTHOR

Simon Ilyushchenko

Copyright (c) 2004 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

my $ftr = 'Bio::Graphics::Feature';
my $segment = $ftr->new(-start=>1,-end=>1000,-name=>'ZK154',-type=>'clone');
my $subseg1 = $ftr->new(-start=>100,-end=>600,-name=>'saw teeth');

my $panel = Bio::Graphics::Panel->new(
                    -grid => 1,
                    -segment => $segment,
                    -key_style => 'bottom');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-interval => 20,
			-width => 20,
			-glyph => 'saw_teeth');

$subseg1->name('frequent saw teeth');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-interval => 0,
			-width => 10,
			-fgcolor => 'red',
			-glyph => 'saw_teeth');

$subseg1->name('dashed line');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-dash_size => 10,
			-space_size => 5,
			-glyph => 'dashed_line');

$subseg1->name('thick colored dashed line with shear');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-fgcolor => 'red',
			-dash_size => 20,
			-space_size => 5,
			-space_color => 'blue',
			-shear => 'yes',
			-linewidth => 2,
			-glyph => 'dashed_line');

$subseg1->name('three letters');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-width => 20,
			-interval => 10,
			-pad_top => 30,
			-glyph => 'three_letters');

$subseg1->name('flag');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-width => 20,
			-text => "ori",
			-height => 30,
			-glyph => 'flag');

$subseg1->name('dumbbell - square ');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-shape_sive => 20,
			-end_shape => "square",
			-height => 20,
			-glyph => 'dumbbell');

$subseg1->name('dumbbell - diamond');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-shape_size => 20,
			-end_shape => "diamond",
			-fgcolor => 'orange',
			-height => 20,
			-glyph => 'dumbbell');

$subseg1->name('dumbbell - tree');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-shape_size => 20,
			-end_shape => "tree",
			-fgcolor => 'green',
			-height => 20,
			-glyph => 'dumbbell');

$subseg1->name('dumbbell - clover');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-shape_size => 20,
			-end_shape => "clover",
			-fgcolor => 'pink',
			-height => 20,
			-glyph => 'dumbbell');


$subseg1->name('dumbbell - star with text');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-shape_size => 40,
			-end_shape => "star",
			-fgcolor => 'red',
			-height => 40,
			-caption => 'Back in USSR',
			-glyph => 'dumbbell');


$subseg1->name('dumbbell - bubble text');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-shape_size => 20,
			-end_shape => "bubble",
			-bubble_text => 'CpG',
			-fgcolor => 'red',
			-height => 20,
			-glyph => 'dumbbell');







open OUT,">glyphs1.png" or die "Couldn't open glyphs1.png for writing: $!";
print OUT $panel->gd->png;
close OUT;


