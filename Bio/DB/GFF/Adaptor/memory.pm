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
  for my $feature (@{$self->{data}}) {
    next unless $feature->{gname} eq $name;
    next unless $feature->{gclass} eq $class;
    return ($feature->{ref},$feature->{start},$feature->{stop},$feature->{strand});
  }
}

sub get_features{
  my $self = shift;
  my ($rangetype,$srcseq,$class,$start,$stop,$types,$sparse,$callback,$automerge) = @_;
  for my $feature (@{$self->{data}}) {
    my $feature_start = $feature->{start};
    my $feature_stop  = $feature->{stop};
    my $feature_ref   = $feature->{ref};
    next unless $feature_ref eq $srcseq;

    if ($rangetype eq 'overlap') {
      next unless $feature_stop >= $start && $feature_start <= $stop;
    } elsif ($rangetype eq 'contains') {
      next unless $feature_start >= $start && $feature_stop <= $stop;
    } elsif ($rangetype eq 'contained_in') {
      next unless $feature_start <= $start && $feature_stop >= $stop;
    } else {
      next unless $feature_start == $start && $feature_stop == $stop;
    }

    my $feature_source = $feature->{source};
    my $feature_method = $feature->{method};
    foreach (@$types) {
      my ($search_method,$search_source) = @$_;
      next if $search_method ne $feature_method;
      next if defined($search_source) && $search_source ne $feature_source;
    }

    # if we get here, then we have a feature that meets the criteria
    # we pass its information to the callback
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

sub setup_load { }
sub finish_load { }

1;

