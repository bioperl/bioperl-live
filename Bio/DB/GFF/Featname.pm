package Bio::DB::GFF::Featname;
use strict;

use overload '""' => 'asString';

sub new    { bless {class=>$_[1],name=>$_[2]},$_[0] }
sub id     {
  my $self = shift;
  return join ':',@{$self}{qw(class name)};
}
sub name   { shift->{name} }
sub class  { shift->{class} }
sub asString { shift->name }

1;
