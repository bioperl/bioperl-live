package Bio::DB::SeqFeature;


=head1 NAME

Bio::DB::SeqFeature -- Normalized feature for use with Bio::DB::SeqFeature::Store

=head1 SYNOPSIS

 use Bio::DB::SeqFeature::Store;

 # Open the sequence database
 my $db      = Bio::DB::SeqFeature::Store->new( -adaptor => 'DBI::mysql',
                                                -dsn     => 'dbi:mysql:test');
 my ($feature)   = $db->get_features_by_name('ZK909');
 my @subfeatures = $feature->get_SeqFeatures();
 my @exons_only  = $feature->get_SeqFeatures('exon');

 # create a new object
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

The Bio::DB::SeqFeature object is the default SeqFeature class stored
in Bio::DB::SeqFeature databases. It implements both the
Bio::DB::SeqFeature::NormalizedFeatureI and
Bio::DB::SeqFeature::NormalizedTableFeatureI interfaces, which means that its
subfeatures, if any, are stored in the database in a normalized
fashion, and that the parent/child hierarchy of features and
subfeatures are also stored in the database as set of tuples. This
provides efficiencies in both storage and retrieval speed.

Typically you will not create Bio::DB::SeqFeature directly, but will
ask the database to do so on your behalf, as described in
L<Bio::DB::SeqFeature::Store>.

=cut

# just like Bio::DB::SeqFeature::NormalizedFeature except that the parent/child 
# relationships are stored in a table in the Bio::DB::SeqFeature::Store

use strict;
use Carp 'croak';
use Bio::DB::SeqFeature::Store;
use base qw(Bio::DB::SeqFeature::NormalizedFeature 
				Bio::DB::SeqFeature::NormalizedTableFeatureI);

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
Bio::DB::SeqFeature::Store database, and the B<-index> option, which
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

sub add_segment {
  my $self = shift;
  $self->_add_segment(0,@_);
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

=head2 add_SeqFeature

 Title   : add_SeqFeature
 Usage   : $flag = $feature->add_SeqFeature(@features)
 Function: Add subfeatures to the feature
 Returns : true if successful
 Args    : list of Bio::SeqFeatureI objects
 Status  : public

Add one or more subfeatures to the feature. For best results,
subfeatures should be of the same class as the parent feature
(i.e. do not try mixing Bio::DB::SeqFeature::NormalizedFeature with
other feature types).

An alias for this method is add_segment().

=cut

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

=head2 object_store

 Title   : object_store
 Usage   : $store = $feature->object_store([$new_store])
 Function: get or set the database handle
 Returns : current database handle
 Args    : new database handle (optional)
 Status  : public

This method will get or set the Bio::DB::SeqFeature::Store object that
is associated with the feature. After changing the store, you should
probably unset the primary_id() of the feature and call update() to ensure
that the object is written into the database as a new feature.

=cut

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

=head2 AUTOLOADED methods

 @subfeatures = $feature->Exon;

If you use an unknown method that begins with a capital letter, then
the feature autogenerates a call to get_SeqFeatures() using the
lower-cased method name as the primary_tag. In other words
$feature-E<gt>Exon is equivalent to:

 @subfeature s= $feature->get_SeqFeatures('exon')

=cut

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

=head2 primary_id

 Title   : primary_id
 Usage   : $id = $feature->primary_id([$new_id])
 Function: get/set the database ID of the feature
 Returns : the current primary ID
 Args    : none
 Status  : public

This method gets or sets the primary ID of the feature in the
underlying Bio::DB::SeqFeature::Store database. If you change this
field and then call update(), it will have the effect of making a copy
of the feature in the database under a new ID.

=cut

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

=head2 Internal methods

=over 4

=item $feature-E<gt>as_string()

Internal method used to implement overloaded stringification.

=item $boolean = $feature-E<gt>type_match(@list_of_types)

Internal method that will return true if the primary_tag of the feature and
source_tag match any of the list of types (in primary_tag:source_tag
format) provided.

=back

=cut

# This adds subfeatures. It has the property of converting the
# provided features into an object like itself and storing them
# into the database. If the feature already has a primary id and
# an object_store() method, then it is not stored into the database,
# but its primary id is reused.
sub _add_segment {
  my $self       = shift;
  my $normalized = shift;

  my $store      = $self->object_store;
  my $store_parentage = eval{$store->can_store_parentage};

  return         $self->SUPER::_add_segment($normalized,@_)
    unless $normalized && $store_parentage;

  my @segments   = $self->_create_subfeatures($normalized,@_);

  my $pos = "@{$self}{'start','stop','ref','strand'}";

  # fix boundaries
  $self->_fix_boundaries(\@segments,1);

  # freakish fixing of our non-standard Target attribute
  $self->_fix_target(\@segments);

  # write our children out
  if ($normalized) {
    $store->add_SeqFeature($self,@segments);
  } else {
    push @{$self->{segments}},@segments;
  }

  # write us back to disk
  $self->update if $self->primary_id && $pos ne "@{$self}{'start','stop','ref','strand'}"; 
}

# segments can be stored directly in the object (legacy behavior)
# or stored in the database
# an optional list of types can be used to specify which types to return
sub get_SeqFeatures {
  my $self         = shift;
  my @types        = @_;

  my @inline_segs  = exists $self->{segments} ? @{$self->{segments}} : ();
  @inline_segs     = grep {$_->type_match(@types)} @inline_segs if @types;
  my $store        = $self->object_store;

  my @db_segs;

  if ($store && $store->can_store_parentage) {
    if (!@types || $store->subfeature_types_are_indexed) {
      @db_segs = $store->fetch_SeqFeatures($self,@types);
    } else {
      @db_segs     = grep {$_->type_match(@types)} $store->fetch_SeqFeatures($self);
    }
  }

  my @segs         = (@inline_segs,@db_segs);
  foreach (@segs) {
    eval {$_->object_store($store)};
  }
  return @segs;
}

sub denormalized_segments {
  my $self = shift;
  return exists $self->{segments} ? @{$self->{segments}} : ();
}

sub denormalized_segment_count {
  my $self = shift;
  return 0 unless exists $self->{segments};
  return scalar @{$self->{segments}};
}

# for Bio::LocationI compatibility
sub is_remote { return }

# for Bio::LocationI compatibility
sub location_type { return 'EXACT' }

# for Bio::DB::GFF compatibility

sub feature_id {shift->primary_id}

1;


__END__

=head1 BUGS

This is an early version, so there are certainly some bugs. Please
use the BioPerl bug tracking system to report bugs.

=head1 SEE ALSO

L<bioperl>,
L<Bio::DB::SeqFeature::Store>,
L<Bio::DB::SeqFeature::Segment>,
L<Bio::DB::SeqFeature::NormalizedFeature>,
L<Bio::DB::SeqFeature::GFF3Loader>,
L<Bio::DB::SeqFeature::Store::DBI::mysql>,
L<Bio::DB::SeqFeature::Store::bdb>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2006 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut
