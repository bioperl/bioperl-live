package Bio::DB::GFF::Typename;

use overload '""' => 'asString';

sub new    { bless {method=>$_[1],source=>$_[2]},$_[0] }
sub method { shift->{method} }
sub source { shift->{source} }
sub asString { "$_[0]->{method}:$_[0]->{source}" }

1;
