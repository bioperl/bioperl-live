package Bio::DB::GFF::Typename;

use overload '""' => 'asString';

sub new    { bless {method=>$_[1],source=>$_[2]},$_[0] }
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
sub asString { "$_[0]->{method}:$_[0]->{source}" }
sub clone { 
  my $self = shift;
  return bless {%$self},ref $self;
}

1;
