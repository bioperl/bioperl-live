package Bio::DB::GFF::Homol;
use strict;

use Bio::DB::GFF::Segment;
use vars qw(@ISA);
@ISA = 'Bio::DB::GFF::Segment';

sub name     { shift->refseq }
sub asString { shift->name }
sub id       {
  my $self = shift;
  return "$self->{class}:$self->{name}";
}

1;
