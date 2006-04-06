package Bio::DB::SeqFeature::LazyTableFeature;

use strict;
use Carp 'croak';
use Bio::DB::SeqFeature::Store;
use base 'Bio::Graphics::Feature';
use base 'Bio::DB::SeqFeature::NormalizedTableFeatureI';
use overload '""' => \&as_string;

my $USE_OVERLOADED_NAMES     = 1;

# some of this is my fault and some of it is changing bioperl API
*merged_segments = *segments = *get_all_SeqFeatures = *sub_SeqFeature = \&get_SeqFeatures;

##### CLASS METHODS ####

sub overloaded_names {
  my $class = shift;
  my $d     = $USE_OVERLOADED_NAMES;
  $USE_OVERLOADED_NAMES = shift if @_;
  $d;
}

##### CONSTRUCTOR ####

sub new {
  my $class = shift;
  my %args  = @_;
  my $db      = $args{-store};
  my $index = exists $args{-index} ? $args{-index} : 1;
  my $self  = $class->SUPER::new(@_);
  $self->object_store($db);

  if ($db) {
    if ($index) {
      $db->store($self); # this will set the primary_id
    } else {
      $db->store_noindex($self); # this will set the primary_id
    }
  }
  $self;
}

##### INSTANCE METHODS ####

sub object_store {
  my $self = shift;
  my $d = $self->{store};
  $self->{store} = shift if @_;
  $d;
}

sub add_SeqFeature {
  my $self = shift;
  $self->_add_segment(1,@_);
}

sub add_segment {
  my $self = shift;
  $self->_add_segment(0,@_);
}

sub seq {
  my $self = shift;
  if (my $store = $self->object_store) {
    return $store->fetch_sequence(@_);
  } else {
    return $self->SUPER::seq(@_);
  }
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
  my $index_subfeatures_policy = $store->index_subfeatures;
  $normalized &&= $store->_can_store_subFeatures;

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
			       );
    }


    elsif (UNIVERSAL::isa($seg,'Bio::SeqFeatureI')) {
      my $score = $seg->score if $seg->can('score');
      my $f = $self->new(-start  => $seg->start,
			 -end    => $seg->end,
			 -strand => $seg->strand,
			 -seq_id => $seg->seq_id,
			 -name   => $seg->display_name,
			 -primary_tag => $seg->primary_tag,
			 -source_tag  => $seg->source,
			 -score       => $score,
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
  my $method  = $index_subfeatures_policy ? 'store' : 'store_noindex';

  if ($normalized) {  # parent/child data is going to be stored in the database

    my @need_loading = grep {!defined $_->primary_id || $_->object_store ne $mystore} @segments;
    if (@need_loading) {
      my $result;
      if ($index_subfeatures_policy) {
	$result = $mystore->store(@need_loading);
      } else {
	$result = $mystore->store_noindex(@need_loading);
      }
      $result or croak "Couldn't store one or more subseqfeatures";
    }
  }

  for my $seg (@segments) {
    $min_start     = $seg->start if $seg->start < $min_start;
    $max_stop      = $seg->end   if $seg->end   > $max_stop;
  }

  # adjust our boundaries, etc.
  $self->start($min_start) if $min_start < $self->start;
  $self->end($max_stop)    if $max_stop  > $self->end;
  $self->{ref}        ||= $segments[0]->seq_id;
  $self->{strand}     ||= $segments[0]->strand;

  # write our children out
  if ($normalized) {
    $mystore->add_SeqFeature($self,@segments);
  } else {
    push @{$self->{segments}},@segments;
  }
  $self->update if $self->primary_id; # write us back to disk
}

sub update {
  my $self = shift;
  $self->object_store->update($self);
}

# segments can be stored directly in the object (legacy behavior)
# or stored in the database
# an optional list of types can be used to specify which types to return
sub get_SeqFeatures {
  my $self         = shift;
  my @types        = @_;

  my @inline_segs  = exists $self->{segments} ? @{$self->{segments}} : ();
  my $store        = $self->object_store;
  return @inline_segs unless $store && $store->_can_store_subFeatures;

  my @db_segs;

  if (!@types || $store->subfeatures_are_indexed) {
    @db_segs = $store->get_SeqFeatures($self,@types);
  } else {
    @db_segs     = grep {$_->type_match(@types)} $store->get_SeqFeatures($self);
  }
  my @segs         = (@inline_segs,@db_segs);
  return @segs;
}

sub denormalized_segments {
  my $self = shift;
  return exists $self->{segments} ? @{$self->{segments}} : ();
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

sub target {
  my $self    = shift;
  my @targets = $self->attributes('Target');
  my @result;
  for my $t (@targets) {
    my ($seqid,$start,$end,$strand) = split /\s+/,$t;
    $strand ||= +1;
    push @result,Bio::DB::SeqFeature::Segment->new($self->object_store,
						   $seqid,
						   $start,
						   $end,
						   $strand);
  }
  return wantarray ? @result : $result[0];
}

sub as_string {
  my $self = shift;
  return overload::StrVal($self) unless $self->overloaded_names;
  my $name  = $self->display_name || $self->load_id || "id=".$self->primary_id;
  my $method = $self->primary_tag;
  my $source= $self->source_tag;
  my $type  = $source ? "$method:$source" : $method;
  return "$type($name)";
}

sub type_match {
  my $self = shift;
  my @types = @_;
  my $method = $self->primary_tag;
  my $source = $self->source_tag;
  for my $t (@types) {
    my ($m,$s) = split /:/,$t;
    return 1 if $method eq $m && (!defined $s || $source eq $s);
  }
  return;
}

1;
