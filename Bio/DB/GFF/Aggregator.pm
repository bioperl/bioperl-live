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
Bio::DB::GFF->new() call.  The various subclasses of
Bio::DB::GFF::Aggregator are tuned for specific common feature types
such as clones, gapped alignments and transcripts.

Instances of Bio::DB::GFF::Aggregator have three attributes:

=over 4

=item method

This is the GFF method field of the composite feature as a whole.  For
example, "transcript" may be used for a composite feature created by
aggregating individual intron, exon and UTR features.

=item main method

Sometimes GFF groups are organized hierarchically, with one feature
logically containing another.  For example, in the C. elegans schema,
methods of type "Sequence:curated" correspond to regions covered by
curated genes.  There can be zero or one main methods.

=item subparts

This is a list of one or more methods that correspond to the component
features of the aggregates.  For example, in the C. elegans database,
the subparts of transcript are "intron", "exon" and "CDS".

=back

Aggregators have two main methods that can be overridden in
subclasses:

=over 4

=item disaggregate()

This method is called by the Adaptor object prior to fetching a list
of features.  The method is passed an associative array containing the
[method,source] pairs that the user has requested, and it returns a
list of raw features that it would like the adaptor to fetch.

=item aggregate()

This method is called by the Adaptor object after it has fetched 
features.  The method is passed a list of raw features and is expected 
to add its composite features to the list.

=back

The disaggregate() and aggregate() methods provided by the base
Aggregator class should be sufficient for many applications.  In this
case, it suffices for subclasses to override the following methods:

=over 4

=item method()

Return the default method for the composite feature as a whole.

=item main_name()

Return the default main method name.

=item part_names()

Return a list of subpart method names.

=back

Provided that method() and part_names() are overridden (and optionally
main_name() as well), then the bare name of the aggregator subclass
can be passed to the -aggregator of Bio::DB::GFF->new().  For example,
this is a small subclass that will aggregate features of type "allele"
and "polymorphism" into an aggregate named "mutant":

  package Bio::DB::GFF::Aggregator::mutant;

  use strict;
  use Bio::DB::GFF::Aggregator;

  use vars '@ISA';
  @ISA = 'Bio::DB::GFF::Aggregator';

  sub method { 'mutant' }

  sub part_names {
    return qw(allele polymorphism);
  }

  1;

Once installed, this aggregator can be passed to Bio::DB::GFF->new()
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
use vars qw($VERSION @ISA);

$VERSION = '0.10';
@ISA = qw(Bio::Root::RootI);

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
  my ($method,$main,$sub_parts) = rearrange([qw(METHOD MAIN_PART SUB_PARTS)],@_);
  $method ||= $main;
  $main   ||= $method;
  return bless {
		method      => $method,
		main_method => $main,
		sub_parts   => $sub_parts,
	       },$class;
}

=head2 disaggregate

 Title   : disaggregate
 Usage   : $arrayref = $a->disaggregate($types,$factory)
 Function: disaggregate type list into components
 Returns : an array reference
 Args    : see below
 Status  : Public

This method is called to disaggregate a list of types into the set of
low-level features to be retrieved from the GFF database.  The list of
types is passed as an array reference containing a series of
[method,source] pairs.  This method synthesizes a new set of
[method,source] pairs, and appends them to the list of requested
types, return the new list as an array reference.

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

  my $sub_features = Bio::DB::GFF->parse_types($self->get_part_names);
  my $main_feature = Bio::DB::GFF->parse_types($self->get_main_name);

  my (@synthetic_types,@unchanged);
  foreach (@$types) {
    my ($method,$source) = @$_;
    if (lc($method) eq $self->method) { # e.g. "transcript"
      push @synthetic_types,map { [$_->[0],$_->[1] || $source] } @$sub_features,@$main_feature;
    }
    else {
      push @unchanged,$_;
    }
  }

  # remember what we're searching for
  $self->components(\@synthetic_types);
  @$types = (@unchanged,@synthetic_types);
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
creating new composite features when appropriate, and appending them to 
the list.  The method result is a new array reference containing all
the raw features plus the new composite ones.

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

  my $main_method = $self->main_name;
  my $matchsub    = $self->match_sub($factory) or return;

  my %aggregates;
  for my $feature (@$features) {
    next unless $feature->group;
    next unless $matchsub->($feature);
    if ($main_method && lc $feature->method eq lc $main_method) {
      $aggregates{$feature->group}{base} ||= $feature->clone;
    } else {
      push @{$aggregates{$feature->group}{subparts}},$feature;
    }
  }

  # aggregate components
  my @result;
  my $pseudo_method = $self->get_method;
  foreach (keys %aggregates) {
    next unless exists $aggregates{$_}{base};
    next unless exists $aggregates{$_}{subparts};
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
    push @result,$base;
  }

  \@result;
}


sub match_sub {
  my $self    = shift;
  my $factory = shift;

  my $types_to_aggregate = $self->components();  # saved from disaggregate call
  return unless @$types_to_aggregate;
  return $factory->make_match_sub($types_to_aggregate);
}

sub components {
  my $self = shift;
  my $d = $self->{components};
  $self->{components} = shift if @_;
  return unless ref $d;
  return wantarray ? @$d : $d;
}

sub get_part_names {
  my $self = shift;
  if ($self->{sub_parts}) {
    return ref $self->{sub_parts} ? @{$self->{sub_parts}} : $self->{sub_parts};
  } else {
    return $self->part_names;
  }
}

sub get_main_name {
  my $self = shift;
  return $self->{main_method} if defined $self->{main_method};
  return $self->main_name;
}

sub get_method {
  my $self = shift;
  return $self->{method} if defined $self->{method};
  return $self->method;
}


# no default method
sub method {
  return;
}

# no default main method
sub main_name {
  return;
}

# no default part names
sub part_names {
  my $self = shift;
  return;
}

1;
