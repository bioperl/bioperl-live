package Bio::DB::GFF::Adaptor::memory;
use strict;

use Bio::DB::GFF;
use Bio::DB::GFF::Util::Rearrange; # for rearrange()
use vars qw($VERSION @ISA);

use constant MAX_SEGMENT => 100_000_000;  # the largest a segment can get

@ISA =  qw(Bio::DB::GFF);
$VERSION = '0.01';

sub new {
  my $class = shift ;
  my ($file) = rearrange([
			  [qw(FILE DIRECTORY)]
			 ],@_);

  # fill in object
  my $self = bless{ data => [] },$class;
  $self->load($file) if $file;
  return $self;
}

sub insert_sequence {
  my $self = shift;
  my($id,$offset,$seq) = @_;
  $self->{dna}{$id} .= $seq;
}

sub get_dna {
  my $self = shift;
  my ($id,$start,$stop,$class) = @_;
  my $reversed = 0;
  if ($start > $stop) {
    $reversed++;
    ($start,$stop) = ($stop,$start);
  }
  my $dna = substr($self->{dna}{$id},$start-1,$stop-$start+1);
  if ($reversed) {
    $dna =~ tr/gatcGATC/ctagCTAG/;
    $dna = reverse $dna;
  }

  $dna;
}

sub load_gff_line {
  my $self = shift;
  my $feature_hash  = shift;
  $feature_hash->{strand} = '+' if $feature_hash->{strand} eq '.'; 
  $feature_hash->{phase} = '+' if $feature_hash->{phase} eq '.';
  push @{$self->{data}},$feature_hash;
}

sub get_abscoords {
  my $self = shift;
  my ($name,$class,$refseq) = @_;
  my %refs;

  # Find all features that have the requested name and class.
  #for my $type (@$typelist) {
  #	my ($method,$source) = @$type;
  #	if defined $method && length $method {
  #	  next unless $feature_method ;
  #	}
  #} 
  # Sort them by reference point.
  for my $feature (@{$self->{data}}) {
    next unless $feature->{gname} eq $name;
    next unless $feature->{gclass} eq $class;
    push @{$refs{$feature->{ref}}},$feature;
  }

  # find out how many reference points we recovered
  my @found_segments;
  foreach my $ref (keys %refs) {
    next if defined($refseq) and $ref ne $refseq;
    my @found = @{$refs{$ref}};
    my ($strand,$start,$stop);
    foreach (@found) {
      $strand ||= $_->{strand};
      $strand = '+' if $strand eq '.'; 
      $start  = $_->{start} if !defined($start) || $start > $_->{start};
      $stop   = $_->{stop}  if !defined($stop)  || $stop  < $_->{stop};
    }
    push @found_segments,[$ref,$class,$start,$stop,$strand];
  }
  return \@found_segments;
}


#sub get_features{
#  my $self = shift;
#  my ($search,$options,$callback) = @_;
#  my @found_features;
#  my (%result,%obj);

#  for my $feature (@{$self->{data}}) {
#    my $feature_start = $feature->{start};
#    my $feature_stop  = $feature->{stop};
#    my $feature_ref   = $feature->{ref};
#    next unless $feature_ref eq $search->{refseq};

#    my $rangetype = $search->{rangetype};
#    if ($rangetype eq 'overlap') {
#      next unless $feature_stop >= $search->{start} && $feature_start <= $search->{stop};
#    } elsif ($rangetype eq 'contains') {
#      next unless $feature_start >= $search->{start} && $feature_stop <= $search->{stop};
#    } elsif ($rangetype eq 'contained_in') {
#      next unless $feature_start <= $search->{start} && $feature_stop >= $search->{stop};
#    } else {
#      next unless $feature_start == $search->{start} && $feature_stop == $search->{stop};
#    }

#    my $feature_source = $feature->{source};
#    my $feature_method = $feature->{method};
    
#    foreach (@{$search->{types}}) {
#      my ($search_method,$search_source) = @$_;
#      next if $search_method ne $feature_method;
#      next if defined($search_source) && $search_source ne $feature_source;
#    }

    # if we get here, then we have a feature that meets the criteria.
    # If we were asked to sort by group, then we just push onto an array
    # of found features and continue.  Otherwise we call the callback
    # immediately.
#    if ($options->{sort_by_group}) {
#      push @found_features,$feature;
#      next;
#    } else {
#      $callback->($feature_ref,
#		  $feature_start,
#		  $feature_stop,
#		  $feature_source,
#		  $feature_method,
#		  $feature->{score},
#		  $feature->{strand},
#		  $feature->{phase},
#		  $feature->{gclass},
#		  $feature->{gname},
#		  $feature->{tstart},
#		  $feature->{tstop}
#		 );
#    }
#  }

#  for my $feature (sort
#		   {"$a->{gclass}:$a->{gname}" cmp "$b->{gclass}:$b->{gname}"
#		  } @found_features) {  # only true if the sort by group option was specified
 #   $callback->(
#		@{$feature}{qw(ref start stop source method score strand phase gclass gname tstart tstop)}
#	       );
#  }
#}
 

sub get_features{
  my $self = shift;
  my ($rangetype,$refseq,$class,$start,$stop,$types,$sparse,$callback,$order_by_group) = @_;
  my @found_features;
  my (%result,%obj);

  for my $feature (@{$self->{data}}) {
    my $feature_start = $feature->{start};
    my $feature_stop  = $feature->{stop};
    my $feature_ref   = $feature->{ref};
    next unless $feature_ref eq $refseq;


     if (defined $start or defined $stop) {
      $start = 0               unless defined($start);
      $stop  = MAX_SEGMENT     unless defined($stop);
    
      if ($rangetype eq 'overlaps') {
	next unless $feature_stop >= $start && $feature_start <= $stop;
      } elsif ($rangetype eq 'contains') {
	next unless $feature_start >= $start && $feature_stop <= $stop;
      } elsif ($rangetype eq 'contained_in') {
	next unless $feature_start <= $start && $feature_stop >= $stop;
      } else {
	next unless $feature_start == $start && $feature_stop == $stop;
      }

    }
    
    my $feature_source = $feature->{source};
    my $feature_method = $feature->{method};

    if (defined $types && @$types){
      next unless _matching_typelist($feature_method,$feature_source,$types);
    } 

    # if we get here, then we have a feature that meets the criteria.
    # If we were asked to sort by group, then we just push onto an array
    # of found features and continue.  Otherwise we call the callback
    # immediately.
    if ($order_by_group) {
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


sub get_types {
  my $self = shift;
  my ($srcseq,$class,$start,$stop,$want_count,$typelist) = @_;
  my(%result,%obj);

  for my $feature (@{$self->{data}}) {
    my $feature_start = $feature->{start};
    my $feature_stop  = $feature->{stop};
    my $feature_ref   = $feature->{ref};
    my $feature_class = $feature->{class};
    my $feature_method = $feature->{method};
    my $feature_source = $feature->{source};
   
    if (defined $srcseq){
      next unless $feature_ref eq $srcseq ;
    }
    
    if (defined $class){ 
      next unless $feature_class eq $class ;
    }
    
     # the requested range should OVERLAP the retrieved features
     if (defined $start or defined $stop) {
      $start = 1           unless defined $start;
      $stop  = MAX_SEGMENT unless defined $stop;
      next unless $feature_stop >= $start && $feature_start <= $stop;
    }
    
    if (defined $typelist && @$typelist){
      next unless _matching_typelist($feature_method,$feature_source,$typelist);
    }

    my $type = Bio::DB::GFF::Typename->new($feature_method,$feature_source);
    $result{$type}++;
    $obj{$type} = $type;

  }   #end features loop
  
  return $want_count ? %result : values %obj;
 
}

sub _matching_typelist{
  
  my ($feature_method,$feature_source,$typelist) = @_; 

  foreach (@$typelist) {
	 my ($search_method,$search_source) = @$_;
	 next if $search_method ne $feature_method;
	 next if defined($search_source) && $search_source ne $feature_source;
	 return 1;
  }
 
  return 0;
}

sub do_initialize { 1; }
sub setup_load { }
sub finish_load { 1; }

1;

