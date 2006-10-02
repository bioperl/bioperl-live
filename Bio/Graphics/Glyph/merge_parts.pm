package Bio::Graphics::Glyph::merge_parts;

use strict;
use base qw(Bio::Graphics::Glyph);

sub merge_parts {
    my ($self,@parts)  = @_;
    
    # This is the largest gap across which adjacent segments will be merged
    my $max_gap = $self->max_gap;

    my $last_part;

    my @sorted_parts = sort {$a->start <=> $b->start} @parts;

    for my $part (@sorted_parts) {
        if ($last_part) {
            my $start  = $part->start;
            my $end    = $part->stop;
            my $score  = $part->score;
            my $pstart = $last_part->start;
            my $pend   = $last_part->stop;
            my $pscore = $last_part->score || 0;
            my $len    = 1 + abs($end - $start);
            my $plen   = 1 + abs($pend - $pstart);

            # weighted average score
            my $new_score = (($score*$len)+($pscore*$plen))/($len+$plen);

            # don't merge if there is a gap > than the allowed size
            my $gap   = abs($start - $pend);
            my $total = abs($end - $pstart);

	    my $last_f = $last_part->feature;
            if ($gap > $max_gap) {
                $last_part = $part;
                next;
            }

            $part->{start}    = $pstart;
            $part->{score}    = $new_score;
            my ($left,$right) = $self->map_pt($pstart,$end+1);
            $part->{left}     = $left;
            $part->{width}    = ($right - $left) + 1;

            # flag the left feature for removal
            $last_part->{remove} = 1;
        }

        $last_part = $part;

    }

    @parts =  grep {!defined $_->{remove}} @parts;

    return @parts;
}

sub max_gap {
    my $self = shift;
    $self->panel->{max_gap} ||= $self->option('max_gap');
    return $self->panel->{max_gap} || $self->calculate_max_gap;
}

sub calculate_max_gap {
    my $self = shift;
    my $segment_length = $self->panel->length;

    # allow more aggressive merging for larger segments
    # by exponentially increasing max_gap
    my $max_gap = ($segment_length/10000)*($segment_length/500);

    $self->panel->{max_gap} = $max_gap;

    return $max_gap;
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::merge_parts - a base class which suppors semantic zooming of scored alignment features

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This is a base class for
Bio::Graphics::Glyph::graded_segments, 
Bio::Graphics::Glyph::heterogeneous_segments
and Bio::Graphics::Glyph::merged_alignment.
It adds internal methods to support semantic zooming of scored
alignment features. It is not intended for end users.

=head1 BUGS

Please report them.

=head1 SEE ALSO

L<Bio::Graphics::Panel>,
L<Bio::Graphics::Track>,
L<Bio::Graphics::Glyph::graded_segments>
L<Bio::Graphics::Glyph::heterogeneous_segments>
L<Bio::Graphics::Glyph::merged_alignment>

=head1 AUTHOR

Sheldon McKay E<lt>mckays@cshl.eduE<gt>

Copyright (c) 2005 Cold Spring Harbor Laboratory

    This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
