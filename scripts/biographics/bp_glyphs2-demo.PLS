#!/usr/bin/perl

#Wed Sep  1 18:54:18 EDT 2004
#Simon Ilyushchenko - demonstrating several new glyphs for gbrowse, part 2.

use strict;

use Bio::Graphics::Panel;

=head1 NAME

bp_glyphs2-demo.pl - Second demo of Bio::Graphics glyphs

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
my $segment = $ftr->new(-start=>1,-end=>400,-name=>'ZK154',-type=>'clone');
my $subseg1 = $ftr->new(-start=>100,-end=>300,-name=>'glyphs 2');

my $panel = Bio::Graphics::Panel->new(
                    -grid => 1,
                    -segment => $segment,
                    -key_style => 'bottom');

$subseg1->name('dumbbell - arrows with arc ');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-width => 50,
			-height => 30,
			-arc => 1,
			-shape_size => 20,
			-fgcolor => 'crimson',
			-end_shape => "arrow",
			-glyph => 'dumbbell');

$subseg1->name('dumbbell - wave ');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-shape_size => 20,
			-end_shape => "wave",
			-height => 20,
			-fgcolor => 'green',
			-glyph => 'dumbbell');

$subseg1->name('two bolts');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-shape_size => 40,
			-bolt_color => 'violet',
			-height => 20,
			-glyph => 'two_bolts');

$subseg1->name('wave');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-height => 10,
			-circle => 1,
			-glyph => 'wave');

$subseg1->name('broken line');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-height => 20,
			-glyph => 'broken_line');

$subseg1->name('tic_tac_toe');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-height => 20,
			-glyph => 'tic_tac_toe');

$subseg1->name('text_in_box');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-height => 20,
			-text_bgcolor => 'yellow',
			-glyph => 'text_in_box');


$subseg1->name('christmas arrow');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-height => 20,
			-fgcolor => 'steelblue',
			-glyph => 'christmas_arrow');


$subseg1->name('pentagram');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-height => 20,
			-glyph => 'pentagram');


$subseg1->name('weighted arrow');

$panel->add_track(segments=>[$subseg1],
            -label => 1,
			-height => 20,
			-fgcolor => 'sienna',
			-glyph => 'weighted_arrow');



open OUT,">glyphs2.png" or die "Couldn't open glyphs2.png for writing: $!";
print OUT $panel->gd->png;
close OUT;


