package Bio::DB::GFF::Adaptor::memory::feature_serializer;


=head1 NAME

Bio::DB::GFF::Adaptor::memory::feature_serializer - utility methods for serializing and deserializing GFF features

=cut

use strict;

require Exporter;
use vars qw(@EXPORT @EXPORT_OK @hash2array_map);
use base qw(Exporter);
@EXPORT_OK = qw(feature2string string2feature @hash2array_map);
@EXPORT = @EXPORT_OK;

@hash2array_map = qw(ref start stop source method score strand phase gclass gname tstart tstop feature_id group_id bin);

sub feature2string {
  my $feature = shift;
  local $^W = 0;
  my @a = @{$feature}{@hash2array_map};
  push @a,map {join "\0",@$_} @{$feature->{attributes}} if $feature->{attributes};
  return join $;,@a;
}

sub string2feature {
  my $string  = shift;
  my (@attributes,%feature);
  (@feature{@hash2array_map},@attributes) = split $;,$string;
  $feature{attributes} = [map {[split "\0",$_]} @attributes];
  undef $feature{group_id};
  \%feature;
}

1;
