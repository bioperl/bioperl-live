package Bio::DB::SeqFeature::LazyFeature;

use strict;
use Carp 'croak';
use base 'Bio::Graphics::Feature';
use base 'Bio::DB::SeqFeature::NormalizedFeatureI';

*sub_SeqFeature = \&segments;
*add_SeqFeature = \&add_segment;

sub new {
  my $class = shift;
  my %args  = @_;
  my $self  = $class->SUPER::new(@_);
  if (my $db = $args{-store}) {
    $db->store($self);
  }
  $self;
}

sub object_store {
  my $self = shift;
  my $d = $self->{store};
  $self->{store} = shift if @_;
  $d;
}

sub add_segment {
  my $self = shift;
  $self->_add_segment(1,@_);
}

sub add_location {
  my $self = shift;
  $self->_add_segment(0,@_);
}

# This adds subfeatures. It has the property of converting the
# provided features into an object like itself and storing them
# into the database. If the feature already has a primary id and
# an object_store() method, then it is not stored into the database,
# but its primary id is reused.
sub _add_segment {
  my $self       = shift;
  my $normalized = shift;

  my $type = $self->{subtype} || $self->{type};
  my $ref   = $self->seq_id;
  my $name  = $self->name;
  my $class = $self->class;
  my $store = $self->object_store;

  my @segments;

  for my $seg (@_) {

    if (UNIVERSAL::isa($seg,ref $self)) {

      if (!$normalized) {  # make sure the object has no lazy behavior
	$seg->primary_id(undef);
	$seg->object_store(undef);
      }
      push @segments,$seg;
    }

    elsif (ref($seg) eq 'ARRAY') {
      my ($start,$stop) = @{$seg};
      next unless defined $start && defined $stop;  # fixes an obscure bug somewhere above us
      my $strand = $self->{strand};

      if ($start > $stop) {
	($start,$stop) = ($stop,$start);
	$strand = -1;
      }
      push @segments,$self->new(-start  => $start,
				-stop   => $stop,
				-strand => $strand,
				-ref    => $ref,
				-type   => $type,
			        -name   => $name,
			        -class  => $class,
				-store  => undef,
			       );
    }


    elsif (UNIVERSAL::isa($seg,'Bio::SeqFeatureI')) {
      my $f = $self->new(-start  => $seg->start,
			 -end    => $seg->end,
			 -strand => $seg->strand,
			 -seq_id => $seg->seq_id,
			 -name   => $seg->display_name,
			 -primary_tag => $seg->primary_tag,
			 -source_tag  => $seg->source,
			 -score       => eval {$seg->score} || undef,
			 -store       => undef,
			);
      for my $tag ($seg->get_all_tags) {
	my @values = $seg->get_tag_values($tag);
	$f->{attributes}{$tag} = \@values;
      }
      push @segments,$f;
    }

    else {
      croak "$seg is neither a Bio::SeqFeatureI object nor an arrayref";
    }
  }

  return unless @segments;

  my $min_start = $self->start ||  999_999_999_999;
  my $max_stop  = $self->end   || -999_999_999_999;

  my $mystore = $self->object_store;

  if ($normalized) {
    my @need_loading = grep {!defined $_->primary_id || $_->object_store ne $mystore} @segments;
    @need_loading && ($mystore->store(@need_loading) or croak "Couldn't store one or more subseqfeatures");
  }

  for my $seg (@segments) {
    $min_start     = $seg->start if $seg->start < $min_start;
    $max_stop      = $seg->end   if $seg->end   > $max_stop;
    my $id_or_seg  = $normalized ? $seg->primary_id : $seg;
    defined $id_or_seg or croak "No primary ID when there should be";
    push @{$self->{segments}},$id_or_seg;
  }


  # adjust our boundaries, etc.
  $self->start($min_start) if $min_start < $self->start;
  $self->end($max_stop)    if $max_stop  > $self->end;
  $self->{ref}        ||= $segments[0]->seq_id;
  $self->{strand}     ||= $segments[0]->strand;

  $self->update; # write us back to disk
}

sub update {
  my $self = shift;
  $self->object_store->store($self);
}

# segments can be either normalized IDs or ordinary feature objects
sub segments {
  my $self = shift;
  my $s     = $self->{segments} or return;
  my $store = $self->object_store;
  my (@ordinary,@ids);
  for (@$s) {
    if (ref ($_)) {
      push @ordinary,$_;
    } else {
      push @ids,$_;
    }
  }
  my @r = (@ordinary,$store->fetch_many(\@ids));
  return @r;
}

sub load_id {
  return shift->attributes('load_id');
}

sub primary_id {
  my $self = shift;
  my $d    = $self->{primary_id};
  $self->{primary_id} = shift if @_;
  $d;
}

1;
