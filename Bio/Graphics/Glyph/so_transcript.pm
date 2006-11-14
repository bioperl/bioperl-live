package Bio::Graphics::Glyph::so_transcript;

# $Id$

use strict;
use base qw(Bio::Graphics::Glyph::processed_transcript);

1;


__END__

=head1 NAME

Bio::Graphics::Glyph::so_transcript - The sequence ontology transcript glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This is a sequence-ontology compatible glyph, which works hand-in-hand
with the so_transcript aggregator in BioPerl.

This glyph is identical to "processed_transcript," which is described
in detail in L<Bio::Graphics::Glyph::processed_transcript>.

=head1 BUGS

Please report them.

=head1 SEE ALSO

L<Bio::Graphics::Panel>,
L<Bio::Graphics::Glyph>,
L<Bio::Graphics::Glyph::processed_transcript>,
L<Bio::DB::GFF::Aggregators::so_transcript>,
L<GD>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>

Copyright (c) 2005 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
