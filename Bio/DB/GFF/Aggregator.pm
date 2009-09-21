=head1 NAME

Bio::DB::GFF::Aggregator -- Aggregate GFF groups into composite features

=head1 SYNOPSIS

 use Bio::DB::GFF;

 my $agg1 = Bio::DB::GFF::Aggregator->new(-method       => 'cistron',
                                          -main_method  => 'locus',
                                          -sub_parts    => ['allele','variant']
                                         );

 my $agg2 = Bio::DB::GFF::Aggregator->new(-method       => 'splice_group',
                                          -sub_parts    => 'transcript');

 my $db      = Bio::DB::GFF->new( -adaptor    => 'dbi:mysql',
			          -aggregator => [$agg1,$agg2],
                                  -dsn        => 'dbi:mysql:elegans42',
				 );


=head1 DESCRIPTION

Bio::DB::GFF::Aggregator is used to aggregate GFF groups into
composite features.  Each composite feature has a "main part", the
top-level feature, and a series of zero or more subparts, retrieved
with the sub_SeqFeature() method.  The aggregator class is designed to
be subclassable, allowing a variety of GFF feature types to be
supported.

The base Bio::DB::GFF::Aggregator class is generic, and can be used to
create specific instances to be passed to the -aggregator argument of
Bio::DB::GFF-E<gt>new() call.  The various subclasses of
Bio::DB::GFF::Aggregator are tuned for specific common feature types
such as clones, gapped alignments and transcripts.

Instances of Bio::DB::GFF::Aggregator have three attributes:

=over 3

=item *

method

This is the GFF method field of the composite feature as a whole.  For
example, "transcript" may be used for a composite feature created by
aggregating individual intron, exon and UTR features.

=item *

main method

Sometimes GFF groups are organized hierarchically, with one feature
logically containing another.  For example, in the C. elegans schema,
methods of type "Sequence:curated" correspond to regions covered by
curated genes.  There can be zero or one main methods.

=item *

subparts

This is a list of one or more methods that correspond to the component
features of the aggregates.  For example, in the C. elegans database,
the subparts of transcript are "intron", "exon" and "CDS".

=back

Aggregators have two main methods that can be overridden in
subclasses:

=over 4

=item *

disaggregate()

This method is called by the Adaptor object prior to fetching a list
of features.  The method is passed an associative array containing the
[method,source] pairs that the user has requested, and it returns a
list of raw features that it would like the adaptor to fetch.

=item *

aggregate()

This method is called by the Adaptor object after it has fetched 
features.  The method is passed a list of raw features and is expected 
to add its composite features to the list.

=back

The disaggregate() and aggregate() methods provided by the base
Aggregator class should be sufficient for many applications.  In this
case, it suffices for subclasses to override the following methods:

=over 4

=item *

method()

Return the default method for the composite feature as a whole.

=item *

main_name()

Return the default main method name.

=item *

part_names()

Return a list of subpart method names.

=back

Provided that method() and part_names() are overridden (and optionally
main_name() as well), then the bare name of the aggregator subclass
can be passed to the -aggregator of Bio::DB::GFF-E<gt>new().  For example,
this is a small subclass that will aggregate features of type "allele"
and "polymorphism" into an aggregate named "mutant":

  package Bio::DB::GFF::Aggregator::mutant;

  use strict;
  use Bio::DB::GFF::Aggregator;

  use base qw(Bio::DB::GFF::Aggregator);

  sub method { 'mutant' }

  sub part_names {
    return qw(allele polymorphism);
  }

  1;

Once installed, this aggregator can be passed to Bio::DB::GFF-E<gt>new()
by name like so:

 my $db      = Bio::DB::GFF->new( -adaptor    => 'dbi:mysql',
			          -aggregator => 'mutant',
                                  -dsn        => 'dbi:mysql:elegans42',
				 );

=head1 API

The remainder of this document describes the public and private
methods implemented by this module.

=cut

package Bio::DB::GFF::Aggregator;

use strict;
use Bio::DB::GFF::Util::Rearrange;  # for rearrange()
use Bio::DB::GFF::Feature;

use base qw(Bio::Root::Root);

my $ALWAYS_TRUE   = sub { 1 };

=head2 new

 Title   : new
 Usage   : $a = Bio::DB::GFF::Aggregator->new(@args)
 Function: create a new aggregator
 Returns : a Bio::DB::GFF::Aggregator object
 Args    : see below
 Status  : Public

This is the constructor for Bio::DB::GFF::Aggregator.  Named arguments 
are as follows:

  -method           the method for the composite feature

  -main_method      the top-level raw feature, if any

  -sub_parts        the list of raw features that will form the subparts
		    of the composite feature (array reference or scalar)

=cut

sub new {
  my $class = shift;
  my ($method,$main,$sub_parts,$whole_object) = rearrange(['METHOD',
							   ['MAIN_PART','MAIN_METHOD'],
							   ['SUB_METHODS','SUB_PARTS'],
							   'WHOLE_OBJECT'
							  ],@_);
  return bless {
		method      => $method,
		main_method => $main,
		sub_parts   => $sub_parts,
		require_whole_object => $whole_object,
	       },$class;
}

=head2 disaggregate

 Title   : disaggregate
 Usage   : $a->disaggregate($types,$factory)
 Function: disaggregate type list into components
 Returns : a true value if this aggregator should be called to reaggregate
 Args    : see below
 Status  : Public

This method is called to disaggregate a list of types into the set of
low-level features to be retrieved from the GFF database.  The list of
types is passed as an array reference containing a series of
[method,source] pairs.  This method synthesizes a new set of
[method,source] pairs, and appends them to the list of requested
types, changing the list in situ.

Arguments:

  $types           reference to an array of [method,source] pairs

  $factory         reference to the Adaptor object that is calling
		   this method

Note that the API allows disaggregate() to remove types from the type
list.  This feature is probably not desirable and may be deprecated in 
the future.

=cut

# this is called at the beginning to turn the pseudo-type 
# into its component feature types
sub disaggregate {
  my $self  = shift;
  my $types = shift;
  my $factory = shift;

  my $sub_features = $factory->parse_types($self->get_part_names);
  my $main_feature = $factory->parse_types($self->get_main_name);

  if (@$types) {
    my (@synthetic_types,@unchanged);
    foreach (@$types) {
      my ($method,$source) = @$_;
      if (lc $method eq lc $self->get_method) { # e.g. "transcript"
	push @synthetic_types,map { [$_->[0],$_->[1] || $source] } @$sub_features,@$main_feature;
      }
      else {
	push @unchanged,$_;
      }
    }
    # remember what we're searching for
    $self->components(\@synthetic_types);
    $self->passthru(\@unchanged);
    @$types = (@unchanged,@synthetic_types);
  }

  # we get here when no search types are listed
  else {
    my @stypes = map { [$_->[0],$_->[1]] }  @$sub_features,@$main_feature;
    $self->components(\@stypes);
    $self->passthru(undef);
  }

  return $self->component_count > 0;
}


=head2 aggregate

 Title   : aggregate
 Usage   : $features = $a->aggregate($features,$factory)
 Function: aggregate a feature list into composite features
 Returns : an array reference containing modified features
 Args    : see below
 Status  : Public

This method is called to aggregate a list of raw GFF features into the
set of composite features.  The method is called an array reference to
a set of Bio::DB::GFF::Feature objects.  It runs through the list,
creating new composite features when appropriate.  The method result
is an array reference containing the composite features.

Arguments:

  $features        reference to an array of Bio::DB::GFF::Feature objects

  $factory         reference to the Adaptor object that is calling
		   this method

NOTE: The reason that the function result contains the raw features as
well as the aggregated ones is to allow queries like this one:

  @features =  $segment->features('exon','transcript:curated');

Assuming that "transcript" is the name of an aggregated feature and
that "exon" is one of its components, we do not want the transcript
aggregator to remove features of type "exon" because the user asked
for them explicitly.

=cut

sub aggregate {
  my $self = shift;
  my $features = shift;
  my $factory  = shift;

  my $main_method = $self->get_main_name;
  my $matchsub    = $self->match_sub($factory) or return;
  my $strictmatch = $self->strict_match();
  my $passthru    = $self->passthru_sub($factory);

  my (%aggregates,@result);
  for my $feature (@$features) {

    if ($feature->group && $matchsub->($feature)) {
      my $key = $strictmatch->{lc $feature->method,lc $feature->source} 
          ? join ($;,$feature->group,$feature->refseq,$feature->source)
          : join ($;,$feature->group,$feature->refseq);
      if ($main_method && lc $feature->method eq lc $main_method) {
	$aggregates{$key}{base} ||= $feature->clone;
      } else {
	push @{$aggregates{$key}{subparts}},$feature;
      }
      push @result,$feature if $passthru && $passthru->($feature);

    } else {
      push @result,$feature;
    }
  }

  # aggregate components
  my $pseudo_method        = $self->get_method;
  my $require_whole_object = $self->require_whole_object;
  foreach (keys %aggregates) {
    if ($require_whole_object && $self->components) {
      next unless $aggregates{$_}{base}; # && $aggregates{$_}{subparts};
    }
    my $base = $aggregates{$_}{base};
    unless ($base) { # no base, so create one
      my $first = $aggregates{$_}{subparts}[0];
      $base = $first->clone;     # to inherit parent coordinate system, etc
      $base->score(undef);
      $base->phase(undef);
    }
    $base->method($pseudo_method);
    $base->add_subfeature($_) foreach @{$aggregates{$_}{subparts}};
    $base->adjust_bounds;
    $base->compound(1);  # set the compound flag
    push @result,$base;
  }
  @$features = @result;
}


=head2 method

 Title   : method
 Usage   : $string = $a->method
 Function: get the method type for the composite feature
 Returns : a string
 Args    : none
 Status  : Protected

This method is called to get the method to be assigned to the
composite feature once it is aggregated.  It is called if the user did
not explicitly supply a -method argument when the aggregator was
created.

This is the method that should be overridden in aggregator subclasses.

=cut

# default method - override in subclasses
sub method {
  my $self = shift;
  $self->{method};
}

=head2 main_name

 Title   : main_name
 Usage   : $string = $a->main_name
 Function: get the method type for the "main" component of the feature
 Returns : a string
 Args    : none
 Status  : Protected

This method is called to get the method of the "main component" of the
composite feature.  It is called if the user did not explicitly supply
a -main-method argument when the aggregator was created.

This is the method that should be overridden in aggregator subclasses.

=cut

# no default main method
sub main_name {
  my $self = shift;
  return;
}

=head2 part_names

 Title   : part_names
 Usage   : @methods = $a->part_names
 Function: get the methods for the non-main various components of the feature
 Returns : a list of strings
 Args    : none
 Status  : Protected

This method is called to get the list of methods of the "main component" of the
composite feature.  It is called if the user did not explicitly supply
a -main-method argument when the aggregator was created.

This is the method that should be overridden in aggregator subclasses.

=cut

# no default part names
sub part_names {
  my $self = shift;
  return;
}

=head2 require_whole_object

 Title   : require_whole_object
 Usage   : $bool = $a->require_whole_object
 Function: see below
 Returns : a boolean flag
 Args    : none
 Status  : Internal

This method returns true if the aggregator should refuse to aggregate
an object unless both its main part and its subparts are present.

=cut

sub require_whole_object {
  my $self = shift;
  my $d    = $self->{require_whole_object};
  $self->{require_whole_object} = shift if @_;
  $d;
}

=head2 match_sub

 Title   : match_sub
 Usage   : $coderef = $a->match_sub($factory)
 Function: generate a code reference that will match desired features
 Returns : a code reference
 Args    : see below
 Status  : Internal

This method is used internally to generate a code sub that will
quickly filter out the raw features that we're interested in
aggregating.  The returned sub accepts a Feature and returns true if
we should aggregate it, false otherwise.

=cut

#' make emacs happy

sub match_sub {
  my $self    = shift;
  my $factory = shift;
  my $types_to_aggregate = $self->components() or return;  # saved from disaggregate call
  return unless @$types_to_aggregate;
  return $factory->make_match_sub($types_to_aggregate);
}

=head2 strict_match

 Title   : strict_match
 Usage   : $strict = $a->strict_match
 Function: generate a hashref that indicates which subfeatures
           need to be tested strictly for matching sources before
           aggregating
 Returns : a hash ref
 Status  : Internal

=cut

sub strict_match {
  my $self = shift;
  my $types_to_aggregate = $self->components();
  my %strict;
  for my $t (@$types_to_aggregate) {
    $strict{lc $t->[0],lc $t->[1]}++ if defined $t->[1];
  }
  \%strict;
}

sub passthru_sub {
  my $self    = shift;
  my $factory = shift;
  my $passthru = $self->passthru() or return;
  return unless @$passthru;
  return $factory->make_match_sub($passthru);
}

=head2 components

 Title   : components
 Usage   : @array= $a->components([$components])
 Function: get/set stored list of parsed raw feature types
 Returns : an array in list context, an array ref in scalar context
 Args    : new arrayref of feature types
 Status  : Internal

This method is used internally to remember the parsed list of raw
features that we will aggregate.  The need for this subroutine is
seen when a user requests a composite feature of type
"clone:cosmid".  This generates a list of components in which the
source is appended to the method, like "clone_left_end:cosmid" and
"clone_right_end:cosmid".  components() stores this information for
later use.

=cut

sub components {
  my $self = shift;
  my $d = $self->{components};
  $self->{components} = shift if @_;
  return unless ref $d;
  return wantarray ? @$d : $d;
}

sub component_count {
  my @c = shift->components;
  scalar @c;
}

sub passthru {
  my $self = shift;
  my $d = $self->{passthru};
  $self->{passthru} = shift if @_;
  return unless ref $d;
  return wantarray ? @$d : $d;
}

sub clone {
  my $self = shift;
  my %new = %{$self};
  return bless \%new,ref($self);
}

=head2 get_part_names

 Title   : get_part_names
 Usage   : @array = $a->get_part_names
 Function: get list of sub-parts for this type of feature
 Returns : an array
 Args    : none
 Status  : Internal

This method is used internally to fetch the list of feature types that
form the components of the composite feature.  Type names in the
format "method:source" are recognized, as are "method" and
Bio::DB::GFF::Typename objects as well.  It checks instance variables
first, and if not defined calls the part_names() method.

=cut

sub get_part_names {
  my $self = shift;
  if ($self->{sub_parts}) {
    return ref $self->{sub_parts} ? @{$self->{sub_parts}} : $self->{sub_parts};
  } else {
    return $self->part_names;
  }
}

=head2 get_main_name

 Title   : get_main_name
 Usage   : $string = $a->get_main_name
 Function: get the "main" method type for this feature
 Returns : a string
 Args    : none
 Status  : Internal

This method is used internally to fetch the type of the "main part" of
the feature.  It checks instance variables first, and if not defined
calls the main_name() method.

=cut

sub get_main_name {
  my $self = shift;
  return $self->{main_method} if defined $self->{main_method};
  return $self->main_name;
}

=head2 get_method

 Title   : get_method
 Usage   : $string = $a->get_method
 Function: get the method type for the composite feature
 Returns : a string
 Args    : none
 Status  : Internal

This method is used internally to fetch the type of the method that
will be assigned to the composite feature once it is synthesized.

=cut

sub get_method {
  my $self = shift;
  return $self->{method} if defined $self->{method};
  return $self->method;
}

1;

=head1 BUGS

None known yet.

=head1 SEE ALSO

L<Bio::DB::GFF>,
L<Bio::DB::GFF::Aggregator::alignment>,
L<Bio::DB::GFF::Aggregator::clone>,
L<Bio::DB::GFF::Aggregator::coding>,
L<Bio::DB::GFF::Aggregator::match>,
L<Bio::DB::GFF::Aggregator::processed_transcript>,
L<Bio::DB::GFF::Aggregator::transcript>,
L<Bio::DB::GFF::Aggregator::none>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

