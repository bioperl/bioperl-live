package Bio::DB::SeqFeature::Store::DBI::Iterator;


=head1 NAME

Bio::DB::SeqFeature::Store::DBI::Iterator - utility methods for creating and iterating over SeqFeature records

=cut

sub new {
  my $class          = shift;
  my ($sth,$store)   = @_;
  return bless {sth   => $sth,
		store => $store
	       },ref($class) || $class;
}

sub next_seq {
  my $self  = shift;
  my $sth   = $self->{sth}   or return;
  my $store = $self->{store} or return;
  my $obj   = $store->_sth2obj($sth);
  if (!$obj) {
    $self->{sth}->finish;
    undef $self->{sth};
    undef $self->{store};
    return;
  }
  return $obj;
}

1;
