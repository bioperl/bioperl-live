package Range;

use strict;
use Carp;

use RangeI;

use vars qw(@ISA);

@ISA = qw(RangeI);

sub new {
  my $thingy = shift;
  my $package = ref($thingy) || $thingy;
  my $self = bless {}, $package;
  my $usageMessage = "Specify exactly two of -start, -end, -length";
  my %args = @_;
  $self->strand($args{-strand} || 0);
  
  if($args{-start} && $args{-end} && $args{-length}) {
    confess $usageMessage;
  }
  
  if($args{-start}) {
    $self->start($args{-start});
    if($args{-end}) {
      $self->end($args{-end});
    } elsif($args{-length}) {
      $self->end($self->start()+$args{-length}-1);
    } else {
      confess $usageMessage;
    }
  } elsif($args{-end} && $args{-length}) {
    $self->end($args{-end});
    $self->start($self->end() - $args{-length} + 1);
  } else {
    confess $usageMessage;
  }
  return $self;
}

sub start {
  my $self = shift;
  @_ ? $self->{start} = shift
     : $self->{start};
}

sub end {
  my $self = shift;
  @_ ? $self->{end} = shift
     : $self->{end};
}

sub strand {
  my $self = shift;
  if(@_) {
    my $val = shift;
    $val =~ tr/+/1/;
    $val =~ tr/-/-1/;
    $val =~ tr/./0/;
    if($val == -1 || $val == 0 || $val == 1 ) {
      $self->{strand} = $val;
    }
  }  
  return $self->{strand};
}

sub length {
  my $self = shift;
  return $self->end() - self->start() + 1;
}

sub toString {
  my $self = shift;
  return  "(" . $self->start() . "," . $self->end() . ")";
}
1;
