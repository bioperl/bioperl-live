package Bio::DB::SeqFeature;

# $Id$

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
Bio::DB::SeqFeature::TableFeatureI interfaces, which means that its
subfeatures, if any, are stored in the database in a normalized
fashion, and that the parent/child hierarchy of features and
subfeatures are also stored in the database as set of tuples. This
provides efficiencies in both storage and retrieval speed.

Typically you will not create Bio::DB::SeqFeature directly, but will
ask the database to do so on your behalf, as described in
L<Bio::DB::SeqFeature::Store>.

=cut

# just like Bio::DB::SeqFeature::NormalizedFeature except that the parent/child relationships are
# stored in a table in the Bio::DB::SeqFeature::Store

use strict;
use Carp 'croak';
use Bio::DB::SeqFeature::Store;
use base qw(Bio::DB::SeqFeature::NormalizedFeature Bio::DB::SeqFeature::NormalizedTableFeatureI);

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

The arguments are the same to Bio::SeqFeature::Generic->new() and
Bio::Graphics::Feature->new(). The most important difference is the
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

  my $store      = $self->object_store;
  return         $self->SUPER::_add_segment($normalized,@_)
    unless $normalized && eval{$store->can_store_parentage};

  my @segments   = $self->_create_subfeatures($normalized,@_);

  my $min_start = $self->start ||  999_999_999_999;
  my $max_stop  = $self->end   || -999_999_999_999;

  for my $seg (@segments) {
    $min_start     = $seg->start if $seg->start < $min_start;
    $max_stop      = $seg->end   if $seg->end   > $max_stop;
  }

  # adjust our boundaries, etc.
  $self->start($min_start) if $min_start < $self->start;
  $self->end($max_stop)    if $max_stop  > $self->end;
  $self->{ref}        ||= $segments[0]->seq_id;
  $self->{strand}     ||= $segments[0]->strand;

  my $pos = "@{$self}{'start','end','ref','strand'}";

  # write our children out
  if ($normalized) {
    $store->add_SeqFeature($self,@segments);
  } else {
    push @{$self->{segments}},@segments;
  }

  # write us back to disk
  $self->update if $self->primary_id && $pos ne "@{$self}{'start','end','ref','strand'}"; 
}

# segments can be stored directly in the object (legacy behavior)
# or stored in the database
# an optional list of types can be used to specify which types to return
sub get_SeqFeatures {
  my $self         = shift;
  my @types        = @_;

  my @inline_segs  = exists $self->{segments} ? @{$self->{segments}} : ();
  my $store        = $self->object_store;
  return @inline_segs unless $store && $store->can_store_parentage;

  my @db_segs;

  if (!@types || $store->subfeatures_are_indexed) {
    @db_segs = $store->fetch_SeqFeatures($self,@types);
  } else {
    @db_segs     = grep {$_->type_match(@types)} $store->fetch_SeqFeatures($self);
  }
  my @segs         = (@inline_segs,@db_segs);
  return @segs;
}

sub denormalized_segments {
  my $self = shift;
  return exists $self->{segments} ? @{$self->{segments}} : ();
}

1;


__END__

=head1 BUGS

This is an early version, so there are certainly some bugs. Please
use the BioPerl bug tracking system to report bugs.

=head1 SEE ALSO

L<bioperl>,
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


