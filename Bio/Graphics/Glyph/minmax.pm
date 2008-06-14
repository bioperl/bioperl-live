package Bio::Graphics::Glyph::minmax;
# $Id$

use strict;
use base qw(Bio::Graphics::Glyph::segments);

sub min_score {
  shift->option('min_score');
}

sub max_score {
  shift->option('max_score');
}

sub minmax {
  my $self = shift;
  my $parts = shift;

  # figure out the colors
  my $max_score = $self->max_score;
  my $min_score = $self->min_score;

  my $do_min = !defined $min_score;
  my $do_max = !defined $max_score;

  if ($do_min or $do_max) {
    my $first = $parts->[0];
    for my $part (@$parts) {
      my $s = eval { $part->feature->score } || $part;
      next unless defined $s;
      $max_score = $s if $do_max && (!defined $max_score or $s > $max_score);
      $min_score = $s if $do_min && (!defined $min_score or $s < $min_score);
    }
  }

  ($min_score,$max_score);
}

sub midpoint {
    my $self    = shift;
    my $default = shift;

    my $pivot = $self->bicolor_pivot;

    if ($pivot eq 'zero') {
	return 0;
    } elsif ($pivot eq 'mean') {
	return eval {$self->series_mean} || $self->SUPER::midpoint($default);
    } elsif  ($pivot =~ /^[\d.eE+-]+$/){
	return $pivot;
    } else {
	$self->SUPER::midpoint($default);
    }
}

sub bicolor_pivot {
    my $self = shift;
    return $self->option('bicolor_pivot');
}

sub pos_color {
    my $self = shift;
    return $self->color('pos_color') || $self->bgcolor;
}

sub neg_color {
    my $self = shift;
    return $self->color('neg_color') || $self->bgcolor;
}



1;

__END__

=head1 NAME

Bio::Graphics::Glyph::minmax - The minmax glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph is a common base class for
L<Bio::Graphics::Glyph::graded_segments> and
L<Bio::Graphics::Glyph::xyplot>.  It adds an internal method named
minmax() for calculating the upper and lower boundaries of scored
features, and is not intended for end users.

=head1 BUGS

Please report them.

=head1 SEE ALSO

L<Bio::Graphics::Panel>,
L<Bio::Graphics::Track>,
L<Bio::Graphics::Glyph::graded_segments>,
L<Bio::Graphics::Glyph::xyplot>,

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>

Copyright (c) 2003 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

