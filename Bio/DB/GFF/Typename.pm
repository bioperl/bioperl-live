package Bio::DB::GFF::Typename;

use strict;
use overload 
  '""'     => 'asString',
  fallback => 1;

# cut down on the number of equivalent objects we have to create
my %OBJECT_CACHE;

sub new    {
  my $package = shift;
  my ($method,$source) = @_;
  return $OBJECT_CACHE{"$method:$source"} ||= bless { method => $method,
						      source => $source
						    },$package;
}

sub method {
  my $self = shift;
  my $d = $self->{method};
  $self->{method} = shift if @_;
  $d;
}

sub source {
  my $self = shift;
  my $d = $self->{source};
  $self->{source} = shift if @_;
  $d;
}

sub asString { join ':',@{$_[0]}{qw(method source)} }

sub clone {
  my $self = shift;
  return bless {%$self},ref $self;
}

1;
