package Bio::DB::GFF;

use Bio::DB::GFF::RelSegment;
use Bio::DB::GFF::Feature;

our $VERSION = '0.25';

use strict;
use Carp 'croak';
use Bio::DB::GFF::Util::Rearrange;

# features() is the pseudonym for overlapping_features()
*features = \&overlapping_features;

sub new {
  my $class   = shift;
  my ($adaptor,$aggregators,$args) = rearrange([
						[qw(ADAPTOR FACTORY)],
						[qw(AGGREGATOR AGGREGATORS)]
						],@_);

  $adaptor    ||= 'dbi::mysqlopt';
  $class = "Bio::DB::GFF::Adaptor::\L${adaptor}\E";
  eval "require $class";
  croak "Unable to load $adaptor adaptor: $@" if $@;

  my $self = $class->new($args);

  # handle the aggregators.
  # aggregators are responsible for creating complex multi-part features
  # from the GFF "group" field.  If none are provided, then we provide a
  # list of the two used in WormBase.
  # Each aggregator can be a scalar or a ref.  In the former case
  # it is treated as a class name to call new() on.  In the latter
  # the aggreator is treated as a ready made object.
  $aggregators = $self->default_aggregators unless defined $aggregators;
  my @a = ref($aggregators) ? @$aggregators : $aggregators;
  my @aggregators;
  for my $a (@a) {
    $self->add_aggregator($a);
  }
  $self;
}

sub error {
  my $self = shift;
  my $g = $self->{error};
  $self->{error} = shift if @_;
  $g;
}

sub debug {
  my $self = shift;
  my $g = $self->{debug};
  $self->{debug} = shift if @_;
  $g;
}

sub aggregators {
  my $self = shift;
  return unless $self->{aggregators};
  return @{$self->{aggregators}};
}

sub add_aggregator {
  my $self       = shift;
  my $aggregator = shift;
  my $list = $self->{aggregators} ||= [];
  if (ref $aggregator) { # an object
    push @$list,$a;
  } else {
    my $class = "Bio::DB::GFF::Aggregator::\L${aggregator}\E";
    eval "require $class";
    croak "Unable to load $aggregator aggregator: $@" if $@;
    push @$list,$class->new();
  }
}

sub default_aggregators {
  my $self = shift;
  return ['transcript','clone','alignment'];
}

sub segment {
  my $self = shift;
  # (see Ace::Sequence::DBI::Segment for all the arguments)
  return $_[0] =~ /^-/ ? Bio::DB::GFF::RelSegment->new(-factory => $self,@_)
                       : Bio::DB::GFF::RelSegment->new($self,@_);
}

# This call is responsible for turning a line of GFF into a
# feature object.
# The $parent argument is a Bio::DB::GFF::Segment object and is used
# to establish the coordinate system for the new feature.
# The $group_hash argument is an hash ref that holds previously-
# generated group objects.
# Other arguments are taken right out of the GFF table.
sub make_feature {
  my $self = shift;
  my ($parent,$group_hash,
      $start,$stop,$method,$source,
      $score,$strand,$phase,
      $group_class,$group_name,
      $tstart,$tstop) = @_;

  my $group;  # undefined
  if (defined $group_class && defined $group_name) {
    $group = $group_hash->{$group_class,$group_name,$tstart,$tstop} 
      ||= $self->make_object($group_class,$group_name,$tstart,$tstop);
  }

  return Bio::DB::GFF::Feature->new_feature($parent,$start,$stop,
					    $method,$source,
					    $score,$strand,$phase,
					    $group);
}

# call to return the DNA string for the indicated region
# real work is done by get_dna()
sub dna {
  my $self = shift;
  my ($id,$start,$stop) = rearrange([
				     [qw(NAME ID REF REFSEQ)],
				     qw(START),
				     [qw(STOP END)],
				     ],@_);
  return unless defined $start && defined $stop;
  $self->get_dna($id,$start,$stop);
}

# call to return the features that overlap the named region
# real work is done by get_features
sub overlapping_features {
  my $self = shift;
  my ($refseq,$start,$stop,$types,$parent,$automerge) = rearrange([
								   [qw(REF REFSEQ)],
								   qw(START),
								   [qw(STOP END)],
								   [qw(TYPE TYPES)],
								   qw(PARENT),
								   [qw(MERGE AUTOMERGE)],
								  ],@_);
  return unless defined $start && defined $stop;
  $automerge = 1 unless defined $automerge;
  $self->_features(0,$refseq,$start,$stop,$types,$parent,$automerge);
}


# The same, except that it only returns features that are completely contained within the
# range (much faster usually)
sub contained_features {
  my $self = shift;
  my ($refseq,$start,$stop,$types,$parent,$automerge) = rearrange([
								   [qw(REF REFSEQ)],
								   qw(START),
								   [qw(STOP END)],
								   [qw(TYPE TYPES)],
								   qw(PARENT),
								   [qw(MERGE AUTOMERGE)],
								  ],@_);
  return unless defined $start && defined $stop;
  $automerge = 1 unless defined $automerge;
  $self->_features(1,$refseq,$start,$stop,$types,$parent,$automerge);
}

sub types {
  my $self = shift;
  my ($refseq,$start,$stop,$enumerate) = rearrange ([
						     [qw(REF REFSEQ)],
						     qw(START),
						     [qw(STOP END)],
						     [qw(ENUMERATE COUNT)],
						     ],@_);
  $self->get_types($refseq,$start,$stop,$enumerate);
}

sub _features {
  my $self = shift;
  my ($range_query,$refseq,$start,$stop,$types,$parent,$automerge) = @_;


  $types = $self->parse_types($types);  # parse out list of types
  my $aggregated_types = $types;         # keep a copy


  # allow the aggregators to operate on the original
  if ($automerge) {
    for my $a ($self->aggregators) {
      $aggregated_types = $a->disaggregate($aggregated_types,$self);
    }
  }

  my (%groups);  # cache groups so that we don't create them unecessarily
  my $features = [];

  my $callback = sub { push @$features,$self->make_feature($parent,\%groups,@_) } if $parent;
  $self->get_features($range_query,$refseq,$start,$stop,$aggregated_types,$callback) ;

  if ($automerge) {
    warn "aggregating...\n" if $self->debug;
    my @aggregated;
    foreach my $a (reverse $self->aggregators) {  # last aggregator gets first shot
      $features = $a->aggregate($features,$self);
    }
  }

  warn "filtering...\n" if $self->debug;

  # remove anything from the features list that was not specifically requested.
  my $match = $self->make_match_sub($types);
  return grep { $match->($_) } @$features;
}

# turn feature types in the format "method:source" into a list of [method,source] refs
sub parse_types {
  my $self  = shift;
  return [] if !@_ or !defined($_[0]);

  my @types = ref($_[0]) ? @{$_[0]} : @_;
  my @type_list = map { [split(':',$_,2)] } @types;
  return \@type_list;
}

# a subroutine that matches features indicated by list of types
sub make_match_sub {
  my $self = shift;
  my $types = shift;

  return sub { 1 } unless ref $types && @$types;

  my @expr;
  for my $type (@$types) {
    my ($method,$source) = @$type;
    $method ||= '.*';
    $source ||= '.*';
    push @expr,"$method:$source";
  }
  my $expr = join '|',@expr;
  return $self->{match_subs}{$expr} if $self->{match_subs}{$expr};

  my $sub =<<END;
sub {
  my \$feature = shift;
  return \$feature->type =~ /^$expr\$/i;
}
END
  my $compiled_sub = eval $sub;
  croak $@ if $@;
  return $self->{match_subs}{$expr} = $compiled_sub;
}

# abstract call to turn a feature into an object, given its class and name
sub make_object {
  my $self = shift;
  my ($class,$name,$start,$stop) = @_;
  return Bio::DB::GFF::Homol->new($self,$name,$class,$start,$stop) if defined $start;
  return Bio::DB::GFF::Featname->new($class,$name);
}

# given a sequence class and name, return its coordinates in format (reference,start,stop,strand)
sub abscoords {
  my $self = shift;
  my ($class,$name) = @_;
  if (!defined($name)) {
    ($name,$class) = ($class,'Sequence');
  }
  $self->get_abscoords($class,$name);
}

# THESE ARE THE ROUTINES THAT NEED TO BE OVERRIDDEN IN SUBCLASSES

sub get_dna {
  my $self = shift;
  my ($id,$start,$stop) = @_;
  croak "get_dna() must be implemented by an adaptor";
}

sub get_features{
  my $self = shift;
  my ($isrange,$refseq,$start,$stop,$types,$callback) = @_;
  croak "get_features() must be implemented by an adaptor";
}


sub get_abscoords {
  my $self = shift;
  my ($class,$name) = @_;
  croak "get_abscoords() must be implemented by an adaptor";
}

sub get_types {
  my $self = shift;
  my ($refseq,$start,$stop,$count) = @_;
  croak "get_types() must be implemented by an adaptor";
}

1;
