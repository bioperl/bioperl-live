package Bio::DB::GFF::Adaptor::memory_iterator;
use strict;
# $Id$
# this module needs to be cleaned up and documented

#use constant STH         => 0;
#use constant CALLBACK    => 1;
#use constant CACHE       => 2;

*next_seq = \&next_feature;

sub new {
  my $class = shift;
  my ($data,$callback) = @_;
  my $pos = 0;
  return bless {data     => $data,
		pos      => $pos,
		callback => $callback,
                cache    => []},$class;
  #return bless [$sth,$callback,[]],$class;
}

sub next_feature {
  my $self = shift;
  return shift @{$self->{cache}} if @{$self->{cache}};
  my $data = $self->{data} or return;
  my $next_feature_pos = $self->{pos}; 
  my $callback = $self->{callback};

  my $features;
  while (1) {
    #if (my $feature = $self->{data}[$next_feature_pos]) {
    if ($next_feature_pos < @{$self->{data}}){
      my $feature = $self->{data}[$next_feature_pos];
      $features = $callback->(@{$feature});
      $self->{pos}++;
      last if $features;
    } else {
      undef $self->{pos};
      undef $self->{data};
      $features = $callback->();
      #last;
      return;
    }
  }
  #return $features;
  $self->{cache} = $features or return;
  shift @{$self->{cache}}; 
}

1;
