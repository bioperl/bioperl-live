package Bio::DB::SeqFeature::Store::Cacher;

use strict;

use base 'Bio::DB::SeqFeature::Store';
use Tie::Cache::LRU ();
use Carp 'croak';

use constant CACHE_SIZE => 500;  # keep up to 1000 objects cached in memory
use constant DEBUG => 0;

sub new {
  no strict 'refs';

  my $adaptor = shift->SUPER::new(@_);

  # what we're doing here is kinda' weird. We insert ourselves into
  # the ISA chain of the adaptor so that our caching fetch() and store()
  # methods are used preferentially.
  my $my_class    = __PACKAGE__;
  my $child_class = ref($adaptor);
  my $glob   = "$child_class\:\:ISA";
  my $isa    = *$glob{ARRAY};
  my %isa    = map {$_=>1} @$isa;
  unshift @$isa,$my_class unless $isa{$my_class};

  my %cache;
  tie %cache,'Tie::Cache::LRU',CACHE_SIZE or croak "Couldn't tie: $!";
  $adaptor->{cache} = \%cache;
  return $adaptor;
}

# look in cache before fetching from database
sub fetch {
  my $self       = shift;
  my $primary_id = shift or croak "usage: fetch(\$primary_id)";
  warn "Cache hit on $primary_id" if DEBUG && exists $self->{cache}{$primary_id};
  return $self->{cache}{$primary_id} ||= $self->SUPER::fetch($primary_id);
}

# when storing objects, keep copies in cache
sub store {
  my $self       = shift;
  my $result = $self->SUPER::store(@_);
  for my $obj (@_) {
    defined (my $id     = eval {$obj->primary_id}) or next;
    $self->{cache}{$id} = $obj;
  }
  $result;
}

# given a statement handler that is expected to return rows of (id,object)
# unthaw each object and return a list of 'em
sub _thaw {
  my $self = shift;
  my ($obj,$primary_id) = @_;
  warn "Cache hit on $primary_id" if DEBUG && exists $self->{cache}{$primary_id};
  return $self->{cache}{$primary_id} ||= $self->SUPER::_thaw($obj,$primary_id);
}

1;
