package Bio::DB::GFF::Feature;

use strict;

use Bio::DB::GFF::Util::Rearrange;
use Bio::DB::GFF::RelSegment;
use Bio::DB::GFF::Featname;
use Bio::DB::GFF::Typename;
use Bio::DB::GFF::Homol;
use Bio::Root::RootI;

use vars qw($VERSION @ISA $AUTOLOAD);
@ISA = qw(Bio::DB::GFF::RelSegment Bio::Root::RootI);

$VERSION = '0.20';

*segments = \&sub_SeqFeature;
*name     = \&group;

sub new_feature {
  my $class   = shift;
  my ($parent,
      $start,$stop,
      $method,$source,$score,
      $fstrand,$phase,
      $group) = @_;

  ($start,$stop) = ($stop,$start) if defined($fstrand) and $fstrand eq '-';
  my $self = bless {
		    %$parent,
		    start  => $start,
		    stop   => $stop,
		    type   => Bio::DB::GFF::Typename->new($method,$source),
		    fstrand => $fstrand,
		    score  => $score,
		    phase  => $phase,
		    group  => $group,
		   },$class;
}

sub method {
  my $self = shift;
  my $d = $self->{type}->method;
  $self->{type}->method(shift) if @_;
  $d;
}
sub source {
  my $self = shift;
  my $d = $self->{type}->source;
  $self->{type}->source(shift) if @_;
  $d;
}
sub score  {
  my $self = shift;
  my $d    = $self->{score};
  $self->{score} = shift if @_;
  $d;
}
sub phase  { 
  my $self = shift;
  my $d    = $self->{phase};
  $self->{phase} = shift if @_;
  $d;
}
sub group  {
  my $self = shift;
  my $d    = $self->{group};
  $self->{group} = shift if @_;
  $d;
}
sub info   { shift->group(@_) }
sub target { shift->group(@_) }

sub clone {
  my $self = shift;
  my $clone = $self->SUPER::clone;

  if (ref(my $t = $clone->type)) {
    my $type = $t->can('clone') ? $t->clone : bless {%$t},ref $t;
    $clone->type($type);
  }

  if (ref(my $g = $clone->group)) {
    my $group = $g->can('clone') ? $g->clone : bless {%$g},ref $g;
    $clone->group($group);
  }

  $clone;
}


sub type   {
  my $self = shift;
  my $d = $self->{type};
  $self->{type} = shift if @_;
  $d;
}

sub add_subfeature {
  my $self    = shift;
  my $feature = shift;
  my $type = $feature->method;
  my $subfeat = $self->{subfeatures}{lc $type} ||= [];
  push @{$subfeat},$feature;
}

# adjust a feature so that its boundaries are synched with its subparts' boundaries.
# this works recursively, so subfeatures can contain other features
sub adjust_bounds {
  my $self = shift;
  my $g = $self->{group};

  if (my $subfeat = $self->{subfeatures}) {
    for my $list (values %$subfeat) {
      for my $feat (@$list) {

	# fix up our bounds to hold largest subfeature
	my($start,$stop,$strand) = $feat->adjust_bounds;
	$self->{fstrand} = $strand unless defined $self->{fstrand};
	if ($self->strand >= 0) {
	  $self->{start} = $start if !defined($self->{start}) || $start < $self->{start};
	  $self->{stop}  = $stop  if !defined($self->{stop})  || $stop  > $self->{stop};
	} else {
	  $self->{start} = $start if !defined($self->{start}) || $start > $self->{start};
	  $self->{stop}  = $stop  if !defined($self->{stop})  || $stop  < $self->{stop};
	}

	# fix up endpoints of targets too (for homologies only)
	my $h = $feat->group;
	next unless $h && $h->isa('Bio::DB::GFF::Homol');
	next unless $g && $g->isa('Bio::DB::GFF::Homol');
	($start,$stop) = ($h->{start},$h->{stop});
	if ($h->strand >= 0) {
	  $g->{start} = $start if !defined($g->{start}) || $start < $g->{start};
	  $g->{stop}  = $stop  if !defined($g->{stop})  || $stop  > $g->{stop};
	} else {
	  $g->{start} = $start if !defined($g->{start}) || $start > $g->{start};
	  $g->{stop}  = $stop  if !defined($g->{stop})  || $stop  < $g->{stop};
	}
      }
    }
  }

  ($self->{start},$self->{stop},$self->strand);
}

sub merged_segments {
  my $self = shift;
  my $type = shift;
  $type ||= '';    # prevent uninitialized variable warnings

  return @{$self->{merged_segs}{$type}} if exists $self->{merged_segs}{$type};

  my @segs = sort {$a->start <=> $b->start} $self->sub_SeqFeature($type);
  # attempt to merge overlapping segments
  my @merged;
  for my $s (@segs) {
    my $previous = $merged[-1];
    if ($previous && $previous->stop+1 >= $s->start) {
      $previous->{stop} = $s->{stop};
      # fix up the target too
      my $g = $previous->{group};
      if ( ref($g) &&  $g->isa('Bio::DB::GFF::Homol')) {
	my $cg = $s->{group};
	$g->{stop} = $cg->{stop};
      }
    } else {
      my $copy = $s->clone;
      push @merged,$copy;
    }
  }
  $self->{merged_segs}{$type} = \@merged;
  return @merged;
}

sub sub_types {
  my $self = shift;
  my $subfeat = $self->{subfeatures} or return;
  return keys %$subfeat;
}

sub sub_SeqFeature {
  my $self = shift;
  my $type = shift;
  my $subfeat = $self->{subfeatures} or return;
  $self->sort_features;
  if ($type) {
    my $features = $subfeat->{lc $type} or return;
    return @{$features};
  } else {
    return map {@{$_}} values %{$subfeat};
  }
}

sub strand {
  my $self = shift;
  return 0 unless $self->{fstrand};
  return 0 unless defined $self->{start};
  return $self->{start} < $self->{stop} ? '+1' : '-1';
}

# sort features
sub sort_features {
  my $self = shift;
  return if $self->{sorted}++;
  my $strand = $self->strand or return;
  my $subfeat = $self->{subfeatures} or return;
  for my $type (keys %$subfeat) {
    $subfeat->{$type} = [sort {$a->start<=>$b->start} @{$subfeat->{$type}}] if $strand > 0;
    $subfeat->{$type} = [sort {$b->start<=>$a->start} @{$subfeat->{$type}}] if $strand < 0;
  }
}

sub asString {
  my $self = shift;
  my $type = $self->method;
  my $id   = $self->group || 'unidentified';
  return join '/',$id,$type,$self->SUPER::asString;
}

sub AUTOLOAD {
  my($pack,$func_name) = $AUTOLOAD=~/(.+)::([^:]+)$/;
  my $sub = $AUTOLOAD;
  my $self = $_[0];

  # ignore DESTROY calls
  return if $func_name eq 'DESTROY';

  # fetch subfeatures if func_name has an initial cap
  return $self->sub_SeqFeature($func_name) if $func_name =~ /^[A-Z]/;

  # error message of last resort
  $self->throw(qq(Can\'t locate object method "$func_name" via package "$pack"));
}

1;

