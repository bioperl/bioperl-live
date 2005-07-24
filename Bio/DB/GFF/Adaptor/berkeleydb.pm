package Bio::DB::GFF::Adaptor::berkeleydb;

use strict;

use Data::Dumper;
use Bio::DB::GFF::Util::Rearrange; # for rearrange()
use Bio::DB::GFF::Util::Binning;

use BerkeleyDB;
use Storable qw(freeze thaw);

# this is the smallest bin (1 K)
use constant MIN_BIN    => 1000;
# this is the largest that any reference sequence can be (100 megabases)
use constant MAX_BIN    => 100_000_000;

#We have to define a limit because Berkeleydb sorts in lexicografic order,
#so all the numbers have to have the same length.	
use constant MAX_NUM_LENGTH => length(MAX_BIN);

use base 'Bio::DB::GFF::Adaptor::memory';

sub setup_load {
  my $self = shift;
	$self->{id} = 1;
	
	my $filename = $ENV{SERVER_SOFTWARE} ? "/tmp/db" :"/tmp/db1" ;
	unlink $filename ;
	$self->{db} = new BerkeleyDB::Btree
				-Filename => $filename,
				-Flags    => DB_CREATE,
				-Property => DB_DUP
		or die "Cannot open file $filename: $! $BerkeleyDB::Error\n" ;		

	$filename = $ENV{SERVER_SOFTWARE} ? "/tmp/recno" :"/tmp/recno1" ;
	$self->{iddb} = new BerkeleyDB::Recno
				-Filename => $filename,
				-Flags    => DB_CREATE,
		or die "Cannot open file $filename: $! $BerkeleyDB::Error\n" ;		
}

sub load_gff_line {

  my ($self, $feat) = @_;

  #Continue if this is the '##sequence-region' header line.
  return if lc($feat->{source}) eq "reference";

  $feat->{strand} = '' if $feat->{strand} && $feat->{strand} eq '.';
  $feat->{phase} = ''  if $feat->{phase}  && $feat->{phase} eq '.';

  my $start = $feat->{start};
  my $stop = $feat->{stop};

  my $bin =  bin($feat->{start},$feat->{stop},MIN_BIN);
  $feat->{bin} = $bin;

  $feat->{feature_id} = $self->{id};
  $bin = $self->normalizeNumber($bin);

  $self->{db}->db_put("__class__". $feat->{gclass}, $self->{id});
  $self->{db}->db_put("__name__".(lc $feat->{gname}), $self->{id}); 
  $self->{db}->db_put("__bin__".$bin, $self->{id});

  for my $attr (@{$feat->{attributes}}) {
	my ($attr_name,$attr_value) = @$attr;
	$self->{db}->db_put("__attr__".$attr_name."__".$attr_value, $self->{id});
  }

  #warn "Storing start $start, stop $stop, bin $bin, id ".$self->{id}."\n";
  $self->{iddb}->db_put($self->{id}, freeze($feat));

  $self->{id}++;
}


# given sequence name, return (reference,start,stop,strand)
sub get_abscoords {
  my $self = shift;
  my ($name,$class,$refseq) = @_;
  my %refs;
  my $regexp;

  if ($name =~ /[*?]/) {  # uh oh regexp time
    $name = quotemeta($name);
    $name =~ s/\\\*/.*/g;
    $name =~ s/\\\?/.?/g;
    $regexp++;
  }
  # Find all features that have the requested name and class.
  # Sort them by reference point.
  my @features = @{$self->retrieve_features(-table => "class", key => $class)};

  foreach my $feature (@features){
    my $no_match_class_name;
    my $empty_class_name;
    if (defined $feature->{gname}){
      my $matches = $regexp ? $feature->{gname} =~ /$name/i : $feature->{gname} eq $name;
      $no_match_class_name = !$matches;  # to accomodate Shuly's interesting logic
    }

    else{
      $empty_class_name = 1;
    }

    if ($no_match_class_name || $empty_class_name){
      my $feature_attributes = $feature->{attributes};
      my $attributes = {Alias => $name};
      if (!$self->_matching_attributes($feature_attributes,$attributes)){
		next;
      }
    }

    push @{$refs{$feature->{ref}}},$feature;
  }

  # find out how many reference points we recovered
  if (! %refs) {
    $self->error("$name not found in database");
    return;
  }

  # compute min and max
  my ($ref) = keys %refs;
  my @found = @{$refs{$ref}};
  my ($strand,$start,$stop);

  my @found_segments;
  foreach my $ref (keys %refs) {
    next if defined($refseq) and $ref ne $refseq;
    my @found = @{$refs{$ref}};
    my ($strand,$start,$stop,$name);
    foreach (@found) {
      $strand ||= $_->{strand};
      $strand = '+' if $strand && $strand eq '.'; 
      $start  = $_->{start} if !defined($start) || $start > $_->{start};
      $stop   = $_->{stop}  if !defined($stop)  || $stop  < $_->{stop};
      $name ||= $_->{gname};
    }
    push @found_segments,[$ref,$class,$start,$stop,$strand,$name];

  }

  return \@found_segments;
}


# Low level implementation of fetching a named feature.
# GFF annotations are named using the group class and name fields.
# May return zero, one, or several Bio::DB::GFF::Feature objects.

=head2 _feature_by_name

 Title   : _feature_by_name
 Usage   : $db->get_features_by_name($class,$name,$callback)
 Function: get a list of features by name and class
 Returns : count of number of features retrieved
 Args    : name of feature, class of feature, and a callback
 Status  : protected

This method is used internally.  The callback arguments are those used
by make_feature().

=cut

sub _feature_by_name {
  my $self = shift;
  my ($class,$name,$location,$callback) = @_;
  $callback || $self->throw('must provide a callback argument');

  #use Devel::StackTrace;
  #warn Devel::StackTrace->new->as_string;

  my $count = 0;
  my $id    = -1;
  my ($use_regexp, $use_glob);

  if ($name =~ /[*?]/) {  # uh oh regexp time
	
	#If there is only one trailing *, do a range search
	if ($name =~ /^([^\*]+)\*$/)
	{
	  $name = $1;
	  $use_glob++;
	}
	else
	{
	  $name = quotemeta($name);
	  $name =~ s/\\\*/.*/g;
	  $name =~ s/\\\?/.?/g;
	  $use_regexp++;
	}
  }
  my @features;
  if ($use_glob)
  {
	my $callback = sub {my $feat = shift; $feat->{gname} =~ /^$name/i};
	@features = @{$self->retrieve_features_range
	  (-table => "name", -start => $name, -do_while => $callback)};
  }
  elsif ($use_regexp)
  {
	my $filter = sub {my $feat = shift; $feat->{gname} =~ /$name/i};
    @features = @{$self->filter_features(-table => "name", -filter => $filter)};
  }
  else
  {
    @features = @{$self->retrieve_features(-table => "name", -key => lc $name)};
  }
  foreach my $feature (@features){
    $id++;
	#next unless ($regexp && $feature->{gname} =~ /$name/i);
    next unless $feature->{gclass} eq $class;

    if ($location) {
      next if $location->[0] ne $feature->{ref};
      next if $location->[1] && $location->[1] > $feature->{stop};
      next if $location->[2] && $location->[2] < $feature->{start};
    }
    $count++;
    $callback->(@{$feature}{qw(
			       ref
			       start
			       stop
			       source
			       method
			       score
			       strand
			       phase
			       gclass
			       gname
			       tstart
			       tstop
				   feature_id
			      )},0
	       );
  }
  return $count;
}

#sub get_feature_by_attribute{
sub _feature_by_attribute{
  my $self = shift;
  my ($attributes,$callback) = @_;
  $callback || $self->throw('must provide a callback argument');
  my $count = 0;
  my $feature_id = -1;
  my $feature_group_id = undef;

  #there could be more than one set of attributes......
  while (my ($key, $value) = each %$attributes) {
	
	my @features = @{$self->retrieve_features
	  (-table => "attr", -key => $key."__".$value)};

    for my $feature (@features) {

	  $feature_id++;
	  for my $attr (@{$feature->{attributes}}) {
			 $callback->($feature->{ref},
			  $feature->{start},
			  $feature->{stop},
			  $feature->{source},
			  $feature->{method},
			  $feature->{score},
			  $feature->{strand},
			  $feature->{phase},
			  $feature->{gclass},
			  $feature->{gname},
		  $feature->{tstart},
		  $feature->{tstop},
			  $feature_id,
		  $feature_group_id);
		 $count++;
		  }
      }
  }

}

sub search_notes {
  my $self = shift;
  my ($search_string,$limit) = @_;
  my (@results, @matches);
  my @words = map {quotemeta($_)} $search_string =~ /(\w+)/g;
  my $search = join '|',@words;

  my $filter = sub {
	my $feature = shift;
    next unless defined $feature->{gclass} && defined $feature->{gname}; # ignore NULL objects
    next unless $feature->{attributes};
    my @attributes = @{$feature->{attributes}};
    my @values     = map {$_->[1]} @attributes;
    my $value      = "@values";
    my $matches    = 0;
    for my $w (@words) {
      my @hits = $value =~ /($w)/ig;
      $matches += @hits;
    }
	
	if ($matches)
	{
	  push @matches, $matches;
	  return -1 if @matches >= $limit;
	}
	
	return $matches;
  };

  my @features = @{$self->filter_features(-table => "attrib__Note", -filter => $filter)};

  for (my $i=0; $i<scalar @matches; $i++)
  {
	my $feature = $features[$i];
	my $matches = $matches[$i];
	
	my $relevance = 10 * $matches;
	my $featname = Bio::DB::GFF::Featname->new($feature->{gclass}=>$feature->{gname});
	my $note;
	$note   = join ' ',map {$_->[1]} grep {$_->[0] eq 'Note'}                @{$feature->{attributes}};
	$note  .= join ' ',grep /$search/,map {$_->[1]} grep {$_->[0] ne 'Note'} @{$feature->{attributes}};
	push @results,[$featname,$note,$relevance];
  }	

  return @results;
}

sub _get_features_by_search_options
{
#The $data argument is not used and is preserved for superclass compatibility
  my ($self, $data,$search,$options) = @_;
  my $count = 0;

  my ($rangetype,$refseq,$class,$start,$stop,$types,$sparse,$order_by_group,$attributes) = 
    (@{$search}{qw(rangetype refseq refclass start stop types)},
    @{$options}{qw(sparse sort_by_group ATTRIBUTES)}) ;
	
  $start = 0               unless defined($start);
  $stop  = MAX_BIN unless defined($stop);
	
  my $cursor = $self->{db}->db_cursor;
  my %found;

  my $bin =  bin($start,$stop,MIN_BIN);  
  $bin = $self->normalizeNumber($bin);

  my @features;
  my $filter = sub {
	my $feature = shift;
	
	my $feature_start = $feature->{start};
	my $feature_stop  = $feature->{stop};
	my $feature_id  = $feature->{feature_id};

	return 0 if $found{$feature_id};
	
	if ($rangetype eq 'overlaps') {
	  return 0 unless $feature_stop >= $start && $feature_start <= $stop;
	} elsif ($rangetype eq 'contains') {
	  return 0 unless $feature_start >= $start && $feature_stop <= $stop;
	} elsif ($rangetype eq 'contained_in') {
	  return 0 unless $feature_start <= $start && $feature_stop >= $stop;
	} else {
	  return 0 unless $feature_start == $start && $feature_stop == $stop;
	}
	
	$found{$feature_id} = 1;
	return 1;
  };

  my $tier = MAX_BIN;
  while ($tier >= MIN_BIN) {
	my ($tier_start,$tier_stop) = (bin_bot($tier,$start),bin_top($tier,$stop));
	warn "Using $tier_start $tier_stop\n";
	if ($tier_start == $tier_stop) {
	  push @features, @{$self->retrieve_features(
		-table => "bin", -key => $tier_start, -filter => $filter)};
	} else {
	  my $callback = sub {my $feat = shift; $feat->{bin} <= $tier_stop};
	  push @features, @{$self->retrieve_features_range
		(-table => "bin", -start => $tier_start, -do_while => $callback, -filter => $filter
		)};
	}
	
	$tier /= 10;
  }

  return \@features;
}

sub retrieve_features
{
  my ($self) = shift;

  my ($table, $key, $filter) = rearrange(['TABLE','KEY','FILTER'],@_);

  my $frozen;
  my @ids = $self->{db}->get_dup("__".$table."__".$key);
  my @result;

  foreach my $id (@ids)
  {
	$self->{iddb}->db_get($id, $frozen);
	my $feat = thaw($frozen);
	next if $filter && !$filter->($feat);
	push @result, $feat;
  }
  return \@result;
}

sub retrieve_features_range
{
  my ($self) = shift;
  my ($table, $start, $do_while, $filter) = rearrange(['TABLE','START','DO_WHILE', 'FILTER'],@_);

  my @result;
  my ($id, $key, $frozen);

  my $cursor = $self->{db}->db_cursor;

  $key = "__".$table."__".$start;  
  my $status = $cursor->c_get($key, $id, DB_SET_RANGE);
  return \@result if $status;

  $self->{iddb}->db_get($id, $frozen);
  my $feature = thaw($frozen);

  while ($do_while->($feature))
  {
	unless ($filter)
	{
	  push @result, $feature;
	}
	else
	{
	  my $filter_result = $filter->($feature);
	  push @result, $feature if $filter_result;
	  last if $filter_result == -1;
	}
	
	$status = $cursor->c_get($key, $id, DB_NEXT);
	last if $status;
	$self->{iddb}->db_get($id, $frozen);
	$feature = thaw($frozen);
  }

  return \@result;
}

sub filter_features
{
  my ($self) = shift;

  my ($table, $filter) = rearrange(['TABLE','FILTER'],@_);

  my @result;
  my ($id, $frozen);

  my $cursor = $self->{iddb}->db_cursor;

  my $status = $cursor->c_get($id, $frozen, DB_NEXT);
  return \@result if $status;

  $self->{iddb}->db_get($id, $frozen);
  my $feature = thaw($frozen);

  while (!$status)
  {
	my $filter_result = $filter->($feature);
	push @result, $feature if $filter_result;
	last if $filter_result == -1;
	$status = $cursor->c_get($id, $frozen, DB_NEXT);
	$feature = thaw($frozen) if $frozen;
  }

  return \@result;
}


sub _basic_features_by_id{
  my $self = shift;
  my ($ids) = @_;

  $ids = [$ids] unless ref $ids =~ /ARRAY/;

  my @result;
  for my $feature_id (@$ids){
	  my $frozen;
	  $self->{iddb}->db_get($feature_id, $frozen);  
	  push @result, thaw($frozen);
  }

  return wantarray() ? @result : $result[0];
}

sub normalizeNumber
{
  my ($self, $num) = @_;
  while ((length $num) < MAX_NUM_LENGTH)
  {
	$num = "0".$num;
  }
  return $num;
}

1;

