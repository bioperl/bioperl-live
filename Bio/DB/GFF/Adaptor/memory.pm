package Bio::DB::GFF::Adaptor::memory;
use strict;

use Bio::DB::GFF;
use Bio::DB::GFF::Util::Rearrange; # for rearrange()
use vars qw($VERSION @ISA);

@ISA =  qw(Bio::DB::GFF);
$VERSION = '0.01';

sub new {
  my $class = shift;
  my ($file) = rearrange([
			  [qw(FILE DIRECTORY)]
			 ],@_);

  # fill in object
  my $self = bless{ data => [] },$class;
  $self->load($file) if $file;
  return $self;
}

sub load_gff_line {
  my $self = shift;
  my $feature_hash  = shift;
  push @{$self->{data}},$feature_hash;
}

sub get_abscoords {
  my $self = shift;
  my ($name,$class,$refseq) = @_;
  my %refs;

  # Find all features that have the requested name and class.
  # Sort them by reference point.
  for my $feature (@{$self->{data}}) {
    next unless $feature->{gname} eq $name;
    next unless $feature->{gclass} eq $class;
    push @{$refs{$feature->{ref}}},$feature;
  }

  # find out how many reference points we recovered
  if (! %refs) {
    $self->error("$name not found in database");
    return;
  } elsif (keys %refs > 1) {
    $self->error("$name has more than one reference sequence in database");
    return;
  }

  # compute min and max
  my ($ref) = keys %refs;
  my @found = @{$refs{$ref}};
  my ($strand,$start,$stop);
  foreach (@found) {
    $strand ||= $_->{strand};
    $start  = $_->{start} if !defined($start) || $start > $_->{start};
    $stop   = $_->{stop}  if !defined($stop)  || $stop  < $_->{stop};
  }
  return ($ref,$class,$start,$stop,$strand);
}

sub get_features{
  my $self = shift;
  my ($search,$options,$callback) = @_;
  my @found_features;

  for my $feature (@{$self->{data}}) {
    my $feature_start = $feature->{start};
    my $feature_stop  = $feature->{stop};
    my $feature_ref   = $feature->{ref};
    next unless $feature_ref eq $search->{refseq};

    my $rangetype = $search->{rangetype};
    if ($rangetype eq 'overlap') {
      next unless $feature_stop >= $search->{start} && $feature_start <= $search->{stop};
    } elsif ($rangetype eq 'contains') {
      next unless $feature_start >= $search->{start} && $feature_stop <= $search->{stop};
    } elsif ($rangetype eq 'contained_in') {
      next unless $feature_start <= $search->{start} && $feature_stop >= $search->{stop};
    } else {
      next unless $feature_start == $search->{start} && $feature_stop == $search->{stop};
    }

    my $feature_source = $feature->{source};
    my $feature_method = $feature->{method};
    foreach (@{$search->{types}}) {
      my ($search_method,$search_source) = @$_;
      next if $search_method ne $feature_method;
      next if defined($search_source) && $search_source ne $feature_source;
    }

    # if we get here, then we have a feature that meets the criteria.
    # If we were asked to sort by group, then we just push onto an array
    # of found features and continue.  Otherwise we call the callback
    # immediately.
    if ($options->{sort_by_group}) {
      push @found_features,$feature;
      next;
    } else {
      $callback->($feature_ref,
		  $feature_start,
		  $feature_stop,
		  $feature_source,
		  $feature_method,
		  $feature->{score},
		  $feature->{strand},
		  $feature->{phase},
		  $feature->{gclass},
		  $feature->{gname},
		  $feature->{tstart},
		  $feature->{tstop}
		 );
    }
  }

  for my $feature (sort
		   {"$a->{gclass}:$a->{gname}" cmp "$b->{gclass}:$b->{gname}"
		  } @found_features) {  # only true if the sort by group option was specified
    $callback->(
		@{$feature}{qw(ref start stop source method score strand phase gclass gname tstart tstop)}
	       );
  }
}

sub setup_load { }
sub finish_load { }

1;

