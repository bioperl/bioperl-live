package Bio::DB::SeqFeature::NormalizedFeature;


=head1 NAME

Bio::DB::SeqFeature::NormalizedFeature -- Normalized feature for use with Bio::DB::SeqFeature::Store

=head1 SYNOPSIS

 use Bio::DB::SeqFeature::Store;
 # Open the sequence database
 my $db      = Bio::DB::SeqFeature::Store->new( -adaptor => 'DBI::mysql',
                                                -dsn     => 'dbi:mysql:test');
 my ($feature)   = $db->get_features_by_name('ZK909');
 my @subfeatures = $feature->get_SeqFeatures();
 my @exons_only  = $feature->get_SeqFeatures('exon');

 # create a new object
 $db->seqfeature_class('Bio::DB::SeqFeature::NormalizedFeature');
 my $new = $db->new_feature(-primary_tag=>'gene',
                            -seq_id     => 'chr3',
                            -start      => 10000,
                            -end        => 11000);

 # add a new exon
 $feature->add_SeqFeature($db->new_feature(-primary_tag=>'exon',
                                           -seq_id     => 'chr3',
                                           -start      => 5000,
                                           -end        => 5551));

=head1 DESCRIPTION

The Bio::DB::SeqFeature::NormalizedFeature object is an alternative
representation of SeqFeatures for use with Bio::DB::SeqFeature::Store
database system. It is identical to Bio::DB::SeqFeature, except that
instead of storing feature/subfeature relationships in a database
table, the information is stored in the object itself. This actually
makes the objects somewhat inconvenient to work with from SQL, but
does speed up access somewhat.

To use this class, pass the name of the class to the
Bio::DB::SeqFeature::Store object's seqfeature_class() method. After
this, $db-E<gt>new_feature() will create objects of type
Bio::DB::SeqFeature::NormalizedFeature. If you are using the GFF3
loader, pass Bio::DB::SeqFeature::Store::GFF3Loader-E<gt>new() the
-seqfeature_class argument:

  use Bio::DB::SeqFeature::Store::GFF3Loader;

  my $store  = connect_to_db_somehow();
  my $loader = Bio::DB::SeqFeature::Store::GFF3Loader->new(
                -store=>$db,
                -seqfeature_class => 'Bio::DB::SeqFeature::NormalizedFeature'
               );

=cut

use strict;
use Carp 'croak';
use base 'Bio::SeqFeature::Lite';
use base 'Bio::DB::SeqFeature::NormalizedFeatureI';
use overload '""' => \&as_string,
              eq  => \&eq,
              ne  => \&ne,
              fallback => 1;

use vars '$AUTOLOAD';

my $USE_OVERLOADED_NAMES     = 1;

# some of this is my fault and some of it is changing bioperl API
*get_all_SeqFeatures = *sub_SeqFeature = *merged_segments = \&segments;

##### CLASS METHODS ####

=head2 new

 Title   : new
 Usage   : $feature = Bio::DB::SeqFeature::NormalizedFeature->new(@args)
 Function: create a new feature
 Returns : the new seqfeature
 Args    : see below
 Status  : public

This method creates and, if possible stores into a database, a new
Bio::DB::SeqFeature::NormalizedFeature object using the specialized
Bio::DB::SeqFeature class.

The arguments are the same to Bio::SeqFeature::Generic-E<gt>new() and
Bio::Graphics::Feature-E<gt>new(). The most important difference is the
B<-store> option, which if present creates the object in a
Bio::DB::SeqFeature::Store database, and he B<-index> option, which
controls whether the feature will be indexed for retrieval (default is
true). Ordinarily, you would only want to turn indexing on when
creating top level features, and off only when storing
subfeatures. The default is on.

Arguments are as follows:

  -seq_id       the reference sequence
  -start        the start position of the feature
  -end          the stop position of the feature
  -display_name the feature name (returned by seqname)
  -primary_tag  the feature type (returned by primary_tag)
  -source       the source tag
  -score        the feature score (for GFF compatibility)
  -desc         a description of the feature
  -segments     a list of subfeatures (see Bio::Graphics::Feature)
  -subtype      the type to use when creating subfeatures
  -strand       the strand of the feature (one of -1, 0 or +1)
  -phase        the phase of the feature (0..2)
  -url          a URL to link to when rendered with Bio::Graphics
  -attributes   a hashref of tag value attributes, in which the key is the tag
                  and the value is an array reference of values
  -store        a previously-opened Bio::DB::SeqFeature::Store object
  -index        index this feature if true

Aliases:

  -id           an alias for -display_name
  -seqname      an alias for -display_name
  -display_id   an alias for -display_name
  -name         an alias for -display_name
  -stop         an alias for end
  -type         an alias for primary_tag

=cut

sub new {
  my $class = shift;
  my %args  = @_;
  my $db      = $args{-store} || $args{-factory};
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

=head2 Bio::SeqFeatureI methods

The following Bio::SeqFeatureI methods are supported:

 seq_id(), start(), end(), strand(), get_SeqFeatures(),
 display_name(), primary_tag(), source_tag(), seq(),
 location(), primary_id(), overlaps(), contains(), equals(),
 intersection(), union(), has_tag(), remove_tag(),
 add_tag_value(), get_tag_values(), get_all_tags()

Some methods that do not make sense in the context of a genome
annotation database system, such as attach_seq(), are not supported.

Please see L<Bio::SeqFeatureI> for more details.

=cut

sub seq {
  my $self = shift;

  require Bio::PrimarySeq unless Bio::PrimarySeq->can('new');

  my ($start,$end) = ($self->start,$self->end);
  if ($self->strand < 0) {
    ($start,$end) = ($end,$start);
  }

  if (my $store = $self->object_store) {
    return Bio::PrimarySeq->new(-seq => $store->fetch_sequence($self->seq_id,$start,$end) || '',
				-id  => $self->display_name);
  } else {
      return $self->SUPER::seq($self->seq_id,$start,$end);
  }
}

sub subseq {
  my $self = shift;
  my ($newstart,$newstop) = @_;
  my $store = $self->object_store or return;
  my ($start,$stop) = ($self->start+$newstart-1,$self->end+$newstop-1);
  if ($self->strand < 0) {
    ($start,$stop) = ($stop,$start);
  }
  my $seq = $store->fetch_sequence($self->seq_id,$start,$stop);
  return Bio::PrimarySeq->new($seq);
}

=head2 add_SeqFeature

 Title   : add_SeqFeature
 Usage   : $flag = $feature->add_SeqFeature(@features)
 Function: Add subfeatures to the feature
 Returns : true if successful
 Args    : list of Bio::SeqFeatureI objects
 Status  : public

Add one or more subfeatures to the feature. For best results,
subfeatures should be of the same class as the parent feature
(i.e. don't try mixing Bio::DB::SeqFeature::NormalizedFeature with
other feature types).

An alias for this method is add_segment().

=cut

sub add_SeqFeature {
  my $self = shift;
  $self->_add_segment(1,@_);
}

=head2 update

 Title   : update
 Usage   : $flag = $feature->update()
 Function: Update feature in the database
 Returns : true if successful
 Args    : none
 Status  : public

After changing any fields in the feature, call update() to write it to
the database. This is not needed for add_SeqFeature() as update() is
invoked automatically.

=cut

sub update {
  my $self = shift;
  my $store = $self->object_store or return;
  $store->store($self);
}

=head2 get_SeqFeatures

 Title   : get_SeqFeature
 Usage   : @subfeatures = $feature->get_SeqFeatures([@types])
 Function: return subfeatures of this feature
 Returns : list of subfeatures
 Args    : list of subfeature primary_tags (optional)
 Status  : public

This method extends the Bio::SeqFeatureI get_SeqFeatures() slightly by
allowing you to pass a list of primary_tags, in which case only
subfeatures whose primary_tag is contained on the list will be
returned. Without any types passed all subfeatures are returned.

=cut


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
  for my $r (@r) {
    eval {$r->object_store($store) };
  }
  return @r;
}

=head2 object_store

 Title   : object_store
 Usage   : $store = $feature->object_store([$new_store])
 Function: get or set the database handle
 Returns : current database handle
 Args    : new database handle (optional)
 Status  : public

This method will get or set the Bio::DB::SeqFeature::Store object that
is associated with the feature. After changing the store, you should
probably unset the feature's primary_id() and call update() to ensure
that the object is written into the database as a new feature.

=cut

sub object_store {
  my $self = shift;
  my $d = $self->{store};
  $self->{store} = shift if @_;
  $d;
}


=head2 overloaded_names

 Title   : overloaded_names
 Usage   : $overload = $feature->overloaded_names([$new_overload])
 Function: get or set overloading of object strings
 Returns : current flag
 Args    : new flag (optional)
 Status  : public

For convenience, when objects of this class are stringified, they are
represented in the form "primary_tag(display_name)". To turn this
feature off, call overloaded_names() with a false value. You can
invoke this on an individual feature object or on the class:

  Bio::DB::SeqFeature::NormalizedFeature->overloaded_names(0);

=cut


sub overloaded_names {
  my $class = shift;
  my $d     = $USE_OVERLOADED_NAMES;
  $USE_OVERLOADED_NAMES = shift if @_;
  $d;
}

=head2 segment

 Title   : segment
 Usage   : $segment = $feature->segment
 Function: return a Segment object corresponding to feature
 Returns : a Bio::DB::SeqFeature::Segment
 Args    : none
 Status  : public

This turns the feature into a Bio::DB::SeqFeature::Segment object,
which you can then use to query for overlapping features. See
L<Bio::DB::SeqFeature::Segment>.

=cut

sub segment  {
  my $self = shift;
  return Bio::DB::SeqFeature::Segment->new($self);
}

### instance methods

=head2 AUTOLOADED methods

 @subfeatures = $feature->Exon;

If you use an unknown method that begins with a capital letter, then
the feature autogenerates a call to get_SeqFeatures() using the
lower-cased method name as the primary_tag. In other words
$feature-E<gt>Exon is equivalent to:

 @subfeature s= $feature->get_SeqFeatures('exon')

If you use an unknown method that begins with Tag_(tagname),
Att_(tagname) Is_(tagname), then it will be the same as calling the
each_tag_value() method with the tagname. In a list context, these
autogenerated procedures return the list of results. In scalar
context, they return the first item in the list!!

=cut


sub AUTOLOAD {
  my($pack,$func_name) = $AUTOLOAD=~/(.+)::([^:]+)$/;
  my $sub = $AUTOLOAD;
  my $self = $_[0];

  # ignore DESTROY calls
  return if $func_name eq 'DESTROY';

  # call attributes if func_name begins with "Tag_" or "Att_":
  if ($func_name =~ /^(Tag|Att|Is)_(\w+)/) {
    my @result = $self->each_tag_value($2);
    return wantarray ? @result : $result[0];
  }

  # fetch subfeatures if func_name has an initial cap
  if ($func_name =~ /^[A-Z]/) {
    return $self->get_SeqFeatures(lc $func_name);
  }

  # error message of last resort
  $self->throw(qq(Can't locate object method "$func_name" via package "$pack"));
}#'


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
  my $store      = $self->object_store;

  my @segments   = $self->_create_subfeatures($normalized,@_);

  # fix boundaries
  $self->_fix_boundaries(\@segments);

  # freakish fixing of our non-standard Target attribute
  $self->_fix_target(\@segments);

  for my $seg (@segments) {
    my $id  = $normalized ? $seg->primary_id : $seg;
    defined $id or $self->throw("No primary ID when there should be");
    push @{$self->{segments}},$id;
  };

  $self->update if $self->primary_id; # write us back to disk
}

sub _fix_boundaries {
  my $self       = shift;
  my $segments   = shift;
  my $normalized = shift;

  my $min_start = $self->start ||  999_999_999_999;
  my $max_stop  = $self->end   || -999_999_999_999;

  for my $seg (@$segments) {
    $min_start     = $seg->start if $seg->start < $min_start;
    $max_stop      = $seg->end   if $seg->end   > $max_stop;
  }

  # adjust our boundaries, etc.
  $self->start($min_start) if $min_start < $self->start;
  $self->end($max_stop)    if $max_stop  > $self->end;
  $self->{ref}           ||= $segments->[0]->seq_id;
  $self->{strand}        ||= $segments->[0]->strand;
}

sub _fix_target {
  my $self = shift;
  my $segs = shift;
  my $normalized = shift;  # ignored for now

  # freakish fixing of our non-standard Target attribute
  if (my $t = ($self->attributes('Target'))[0]) {
    my ($seqid,$tstart,$tend,$strand) = split /\s+/,$t;
    if (defined $tstart && defined $tend) {
	my $min_tstart = $tstart;
	my $max_tend   = $tend;
	for my $seg (@$segs) {
	    my $st = ($seg->attributes('Target'))[0] or next;
	    (undef,$tstart,$tend) = split /\s+/,$st;
	    next unless defined $tstart && defined $tend;
	    $min_tstart     = $tstart if $tstart < $min_tstart;
	    $max_tend       = $tend   if $tend   > $max_tend;
	}
	if ($min_tstart < $tstart or $max_tend > $tend) {
	    $self->{attributes}{Target}[0] = join ' ',($seqid,$min_tstart,$max_tend,$strand||'');
	}
    }
  }
}

# undo the load_id and Target hacks on the way out
sub format_attributes {
  my $self   = shift;
  my $parent      = shift;
  my $fallback_id = shift;

  my $load_id   = $self->load_id || '';

  my $targobj = ($self->attributes('Target'))[0];
  # was getting an 'Use of uninitialized value with split' here, changed to cooperate -cjf 7/10/07
  my ($target)  = $targobj ? split /\s+/,($self->attributes('Target'))[0] : ('');
  my @tags = $self->all_tags;
  my @result;
  for my $t (@tags) {
    my @values = $self->each_tag_value($t);

    # This line prevents Alias from showing up if it matches the load id, but this is not good
    # @values = grep {$_ ne $load_id && $_ ne $target} @values if $t eq 'Alias';

    # these are hacks, which we don't want to appear in the file
    next if $t eq 'load_id';
    next if $t eq 'parent_id';

    foreach (@values) { s/\s+$// } # get rid of trailing whitespace
    push @result,join '=',$self->escape($t),join(',', map {$self->escape($_)} @values) if @values;
  }
  my $id         = $self->primary_id || $fallback_id;
  my $parent_id;
  if (@$parent) {
      $parent_id  = join (',',map {$self->escape($_)} @$parent);
  }
  my $name = $self->display_name;
  unshift @result,"ID=".$self->escape($id)                       if defined $id;
  unshift @result,"Parent=".$parent_id                           if defined $parent_id;
  unshift @result,"Name=".$self->escape($name)                   if defined $name;
  return join ';',@result;
}

sub _create_subfeatures {
  my $self = shift;
  my $normalized = shift;

  my $type = $self->{subtype} || $self->{type};
  my $ref   = $self->seq_id;
  my $name  = $self->name;
  my $class = $self->class;
  my $store = $self->object_store;
  my $source = $self->source;

  if ($normalized) {
    $store or $self->throw("Feature must be associated with a Bio::DB::SeqFeature::Store database before attempting to add subfeatures to a normalized object");
  }

  my $index_subfeatures_policy = eval{$store->index_subfeatures};

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
				-source => $source,
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
			 -source => $source,
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

  if ($normalized && $store) {  # parent/child data is going to be stored in the database

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

=head2 load_id

 Title   : load_id
 Usage   : $id = $feature->load_id
 Function: get the GFF3 load ID
 Returns : the GFF3 load ID (string)
 Args    : none
 Status  : public

For features that were originally loaded by the GFF3 loader, this
method returns the GFF3 load ID. This method may not be supported in
future versions of the module.

=cut

sub load_id {
  return (shift->attributes('load_id'))[0];
}


=head2 notes

 Title   : notes
 Usage   : @notes = $feature->notes
 Function: get contents of the GFF3 Note tag
 Returns : List of GFF3 Note tags
 Args    : none
 Status  : public

For features that were originally loaded by the GFF3 loader, this
method returns the contents of the Note tag as a list. This is a
convenience for Bio::Graphics, which looks for notes() when it
constructs a default description line.

=cut

sub notes {
  return shift->attributes('Note');
}

=head2 primary_id

 Title   : primary_id
 Usage   : $id = $feature->primary_id([$new_id])
 Function: get/set the feature's database ID
 Returns : the current primary ID
 Args    : none
 Status  : public

This method gets or sets the primary ID of the feature in the
underlying Bio::DB::SeqFeature::Store database. If you change this
field and then call update(), it will have the effect of making a copy
of the feature in the database under a new ID.

=cut

sub primary_id {
  my $self = shift;
  my $d    = $self->{primary_id};
  $self->{primary_id} = shift if @_;
  $d;
}

=head2 target

 Title   : target
 Usage   : $segment = $feature->target
 Function: return the segment correspondent to the "Target" attribute
 Returns : a Bio::DB::SeqFeature::Segment object
 Args    : none
 Status  : public

For features that are aligned with others via the GFF3 Target
attribute, this returns a segment corresponding to the aligned
region. The CIGAR gap string is not yet supported.

=cut

sub target {
  my $self    = shift;
  my @targets = $self->attributes('Target');
  my @result;
  for my $t (@targets) {
    my ($seqid,$start,$end,$strand) = split /\s+/,$t;
    $strand ||= '';
    $strand = $strand eq '+' ? 1
              : $strand eq '-' ? -1
	      : 0;
    push @result,Bio::DB::SeqFeature::Segment->new($self->object_store,
						   $seqid,
						   $start,
						   $end,
						   $strand);
  }
  return wantarray ? @result : $result[0];
}

=head2 Internal methods

=over 4

=item $feature-E<gt>as_string()

Internal method used to implement overloaded stringification.

=item $boolean = $feature-E<gt>type_match(@list_of_types)

Internal method that will return true if the feature's primary_tag and
source_tag match any of the list of types (in primary_tag:source_tag
format) provided.

=back

=cut

sub as_string {
  my $self = shift;
  return overload::StrVal($self) unless $self->overloaded_names;
  my $name   = $self->display_name || $self->load_id;
  $name    ||= "id=".$self->primary_id if $self->primary_id;
  $name    ||= "<unnamed>";
  my $method = $self->primary_tag;
  my $source= $self->source_tag;
  my $type  = $source ? "$method:$source" : $method;
  return "$type($name)";
}

sub eq {
  my $self = shift;
  my $b    = shift;
  my $store1 = $self->object_store;
  my $store2 = eval {$b->object_store} || '';
  return $store1 eq $store2 && $self->primary_id eq $b->primary_id;
}

sub ne {
  my $self = shift;
  return !$self->eq(shift);
}

# completely case insensitive
sub type_match {
  my $self = shift;
  my @types = @_;
  my $method = lc $self->primary_tag;
  my $source = lc $self->source_tag;
  for my $t (@types) {
    my ($m,$s) = map {lc $_} split /:/,$t;
    return 1 if $method eq $m && (!defined $s || $source eq $s);
  }
  return;
}

sub segments { shift->get_SeqFeatures(@_) }

1;


__END__

=head1 BUGS

This is an early version, so there are certainly some bugs. Please
use the BioPerl bug tracking system to report bugs.

=head1 SEE ALSO

L<bioperl>,
L<Bio::DB::SeqFeature>,
L<Bio::DB::SeqFeature::Store>,
L<Bio::DB::SeqFeature::Segment>,
L<Bio::DB::SeqFeature::GFF3Loader>,
L<Bio::DB::SeqFeature::Store::DBI::mysql>,
L<Bio::DB::SeqFeature::Store::bdb>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2006 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut


