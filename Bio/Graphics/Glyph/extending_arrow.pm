package Bio::Graphics::Glyph::extending_arrow;

use strict;
use vars '@ISA';
use Bio::Graphics::Glyph::anchored_arrow;
@ISA = 'Bio::Graphics::Glyph::anchored_arrow';

=head1 NAME

Bio::Graphics::Glyph::extending_arrow -- The "extending arrow" glyph

=head1 SYNOPSIS

This is deprecated.  Use L<Bio::Graphics::Glyph::anchored_arrow>
instead.

=head1 DESCRIPTION

This glyph was designed to show a segment that goes beyond the panel.
If the segment is contained within the panel, a vertical base is
shown.  Otherwise, an arrow is shown.

Also see the arrow glyph.

=head2 OPTIONS

See L<Bio::Graphics::Glyph::anchored_arrow>.  This glyph has been
deprecated.

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

Originally by Shengqiang Shu.  Temporarily deprecated by Lincoln
Stein.

Copyright (c) 2001 Berkeley Drosophila Genome Project

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
