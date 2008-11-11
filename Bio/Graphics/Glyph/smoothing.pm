package Bio::Graphics::Glyph::smoothing;

use strict;

use constant SMOOTHING  => 'mean';

sub get_smoothing {
  my $self = shift;
  return 'none' if $self->smooth_window == 1;
  return $self->option('smoothing') or SMOOTHING;
}

sub smooth_window {
  my $self    = shift;

  my $smooth_window = $self->option('smoothing_window') 
                    || $self->option('smoothing window'); # drat!
  return $smooth_window if defined $smooth_window; 

  my $start = $self->smooth_start;
  my $end   = $self->smooth_end;

  $smooth_window = int (($end - $start)/(2*$self->width));
  $smooth_window = 1 unless $smooth_window > 2;
  return $smooth_window;
}

sub smooth_start {
  my $self = shift;
  my ($start) = sort {$b<=>$a} ($self->feature->start,$self->panel->start);
  return $start;
}

sub smooth_end {
  my $self = shift;
  my ($end) = sort {$a<=>$b} ($self->feature->end,$self->panel->end);
  return $end;
}

1;

