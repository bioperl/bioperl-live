package Bio::DB::SeqFeature::LazyFeature;

use strict;
use Carp 'croak';
use base 'Bio::Graphics::FeatureBase';
use base 'Bio::DB::SeqFeature::NormalizedFeatureI';
use overload '""' => \&as_string;

my $USE_OVERLOADED_NAMES     = 1;

# some of this is my fault and some of it is changing bioperl API
*get_all_SeqFeatures = *sub_SeqFeature = *merged_segments = \&segments;

##### CLASS METHODS ####

sub new {
  my $class = shift;
  my %args  = @_;
  my $db      = $args{-store};
  my $index = exists $args{-index} ? $args{-index} : 1;
  my $self  = $class->SUPER::new(@_);

  if ($db) {
    if ($index) {
      $db->store($self); # this will set the primary_id
    } else {
      $db->store_noindex($self); # this will set the primary_id
    }
    $self->object_store($db);
  }
  $self;
}

sub overloaded_names {
  my $class = shift;
  my $d     = $USE_OVERLOADED_NAMES;
  $USE_OVERLOADED_NAMES = shift if @_;
  $d;
}

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

# This adds subfeatures. It has the property of converting the
# provided features into an object like itself and storing them
# into the database. If the feature already has a primary id and
# an object_store() method, then it is not stored into the database,
# but its primary id is reused.
sub _add_segment {
  my $self       = shift;
  my $normalized = shift;

  my @segments   = $self->_create_subfeatures($normalized,@_);

  my $min_start = $self->start ||  999_999_999_999;
  my $max_stop  = $self->end   || -999_999_999_999;

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
  $self->{ref}           ||= $segments[0]->seq_id;
  $self->{strand}        ||= $segments[0]->strand;

  $self->update if $self->primary_id; # write us back to disk
}

sub _create_subfeatures {
  my $self = shift;
  my $normalized = shift;

  my $type = $self->{subtype} || $self->{type};
  my $ref   = $self->seq_id;
  my $name  = $self->name;
  my $class = $self->class;
  my $store = $self->object_store or return;

  my $index_subfeatures_policy = $store->index_subfeatures;

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

  if ($normalized) {  # parent/child data is going to be stored in the database

    my @need_loading = grep {!defined $_->primary_id || $_->object_store ne $store} @segments;
    if (@need_loading) {
      my $result;
      if ($index_subfeatures_policy) {
	$result = $store->store(@need_loading);
      } else {
	$result = $store->store_noindex(@need_loading);
      }
      $result or croak "Couldn't store one or more subseqfeatures";
    }
  }

  return @segments;
}

sub update {
  my $self = shift;
  $self->object_store->store($self);
}

# segments can be either normalized IDs or ordinary feature objects
sub get_SeqFeatures {
  my $self = shift;
  my @types        = @_;

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
  my @r = grep {$_->type_match(@types)} (@ordinary,$store->fetch_many(\@ids));
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

sub segments { shift->get_SeqFeatures(@_) }

1;
