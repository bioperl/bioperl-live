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
 my $new = Bio::DB::SeqFeature->new(-primary_tag=>'gene',
                                    -seq_id     => 'chr3',
                                    -start      => 10000,
                                    -end        => 11000,
                                    -store      => $db);

 # add a new exon
 $feature->add_SeqFeature(Bio::SeqFeature::Generic->new(-primary_tag=>'exon',
                                                        -seq_id     => 'chr3',
                                                        -start      => 5000,
                                                        -end        => 5551));

=head1 DESCRIPTION

The Bio::DB::SeqFeature object 

=cut

# just like Bio::DB::SeqFeature::NormalizedFeature except that the parent/child relationships are
# stored in a table in the Bio::DB::SeqFeature::Store

use strict;
use Carp 'croak';
use Bio::DB::SeqFeature::Store;
use base qw(Bio::DB::SeqFeature::NormalizedFeature Bio::DB::SeqFeature::NormalizedTableFeatureI);

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


