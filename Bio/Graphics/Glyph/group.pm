package Bio::Graphics::Glyph::group;

use strict;
use vars '@ISA';
use Bio::Graphics::Glyph::segmented_keyglyph;
@ISA = 'Bio::Graphics::Glyph::segmented_keyglyph';

# group sets connector to 'dashed'
sub connector {
  my $self = shift;
  my $super = $self->SUPER::connector(@_);
  return $super if $self->all_callbacks;
  return 'dashed' unless defined($super) && ($super eq 'none' or !$super);
}

# we don't label group (yet)
sub label { 0 }

sub new {
  my $self = shift->SUPER::new(@_);
  # reset our parts to level zero
  foreach (@{$self->{parts}}) {
    $_->{level} = 0;
  }
  $self;
}

# don't allow simple bumping in groups -- it looks terrible...
sub bump {
  my $bump = shift->SUPER::bump(@_);
  return unless defined $bump;
  return 1  if $bump >  1;
  return -1 if $bump < -1;
  return $bump;
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::group - The "group" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph is used internally by Bio::Graphics::Panel for laying out
groups of glyphs that move in concert.  It should not be used
explicitly.

=head1 BUGS

Please report them.

=head1 SEE ALSO

L<Ace::Sequence>, L<Ace::Sequence::Feature>, L<Bio::Graphics::Panel>,
L<Bio::Graphics::Track>, L<Bio::Graphics::Glyph::anchored_arrow>,
L<Bio::Graphics::Glyph::arrow>,
L<Bio::Graphics::Glyph::box>,
L<Bio::Graphics::Glyph::primers>,
L<Bio::Graphics::Glyph::segments>,
L<Bio::Graphics::Glyph::toomany>,
L<Bio::Graphics::Glyph::transcript>,

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
