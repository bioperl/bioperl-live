package Bio::DB::GFF::Adaptor::berkeleydb;

use strict;

use Data::Dumper;
use Bio::DB::GFF::Util::Rearrange; # for rearrange()
use Bio::DB::GFF::Util::Binning;
use Bio::DB::Fasta;

use DB_File;
use File::Path 'mkpath';
use Storable qw(freeze thaw);

# this is the smallest bin (1 K)
use constant MIN_BIN    => 1000;
# this is the largest that any reference sequence can be (100 megabases)
use constant MAX_BIN    => 100_000_000;
use constant MAX_SEGMENT => 1_000_000_000;  # the largest a segment can get

#We have to define a limit because Berkeleydb sorts in lexicografic order,
#so all the numbers have to have the same length.	
use constant MAX_NUM_LENGTH => length(MAX_BIN);

use base 'Bio::DB::GFF::Adaptor::memory';

sub new {
  my $class = shift ;
  my ($dbdir,$preferred_groups) = rearrange([
					     [qw(DSN DB DIR DIRECTORY)],
					     'PREFERRED_GROUPS',
					    ],@_);
  $dbdir ||= $ENV{TMPDIR} ? "$ENV{TMPDIR}/test" : "/tmp/test";
  my $self = bless {},$class;
  $self->dsn($dbdir);
  $self->preferred_groups($preferred_groups) if defined $preferred_groups;
  $self->_open_databases();
  return $self;
}

sub _open_databases {
  my $self = shift;
  my $dsn  = $self->dsn;
  unless (-d $dsn) {  # try to create the directory
    mkpath($dsn) or die "Couldn't create database directory $dsn: $!";
  }
  my (%db,%iddb);
  local $DB_BTREE->{flags} = R_DUP;
  tie(%db,'DB_File',$self->_index_file,O_RDWR|O_CREAT,0666,$DB_BTREE)
    or die "Couldn't tie ",$self->_index_file,": $!";
  tie(%iddb,'DB_File',$self->_hash_file,O_RDWR|O_CREAT,0666,$DB_HASH)
    or die "Couldn't tie ",$self->_hash_fle,": $!";
  $self->{db}   = \%db;
  $self->{iddb} = \%iddb;

  if (-e $self->_fasta_file) {
    my $dna_db = Bio::DB::Fasta->new($self->_fasta_file) or $self->throw("Can't reindex sequence file: $@");
    $self->dna_db($dna_db);
  }

}

sub _close_databases {
  my $self = shift;
  delete $self->{db};
  delete $self->{iddb};
}

sub _delete_features {
  my $self        = shift;
  my @feature_ids = @_;
  my $removed = 0;
  my $last_id = $self->{db}{__index__};
  for my $id (@feature_ids) {
    next unless $id >= 0 && $id < $last_id;
    my $feat  = thaw ($self->{iddb}{$id}) or next;
    delete $self->{iddb}{$id};  # delete data
    my @keys = $self->_secondary_keys($feat);
    $self->db->del_dup($_,$id) foreach @keys;
    $removed++;
  }
  $removed;
}

sub _secondary_keys {
  my $self = shift;
  my $feat = shift;
  return (
		"__class__". $feat->{gclass},
		"__name__".(lc $feat->{gname}),
		"__bin__"."$feat->{ref}$;$feat->{bin}",
		"__type__".join(':',$feat->{method},$feat->{source}),
		map {"__attr__".$_->[0]."__".$_->[1]} @{$feat->{attributes}}
	  );
}

sub _delete {
  my $self        = shift;
  my $delete_spec = shift;
  return $self->SUPER::_delete($delete_spec) if @{$delete_spec->{segments}} or @{$delete_spec->{types}};
  $self->throw("This operation would delete all feature data and -force not specified")
    unless $delete_spec->{force};
  my $deleted = $self->{db}{__index__};
  %{$self->{iddb}} = ();
  %{$self->{db}}   = ();
  $deleted;
}

sub do_initialize {
  my $self  = shift;
  my $erase = shift;
  if ($erase) {
    $self->_close_databases;
    unlink $self->_index_file;
    unlink $self->_hash_file;
    unlink $self->_fasta_file;
    unlink $self->_fasta_file.'.index';
    $self->_open_databases;
    $self->_next_id(0);
  }
  1;
}

sub load_sequence {
  my $self = shift;
  my ($io_handle,$id) = @_;
  my $file = $self->_fasta_file;
  my $loaded = 0;

  open (F,">>$file") or $self->throw("Couldn't open $file for writing: $!");

  if (defined $id) {
    print F ">$id\n";
    $loaded++;
  }

  while (<$io_handle>) {
    $loaded++ if /^>/;
    print F $_;
  }
  close F;
  my $dna_db = Bio::DB::Fasta->new($file) or $self->throw("Can't reindex sequence file: $@");
  $self->dna_db($dna_db);
  return $loaded;
}

# with duplicates enabled, we cannot simply do $db->{__index__}++;
sub _next_id {
  my $self = shift;
  my $db = $self->{db};
  if (@_) {
    delete $db->{__index__};
    return $db->{__index__} = shift;
  } else {
    my $index = ${db}->{__index__};
    delete $db->{__index__};
    $db->{__index__} = $index + 1;
    return $index;
  }
}

sub _index_file {
  my $self = shift;
  return $self->dsn . "/features.btree";
}

sub _hash_file {
  my $self = shift;
  return $self->dsn . "/features.hash";
}

sub _fasta_file {
  my $self = shift;
  return $self->dsn . "/sequence.fa";
}

sub db {
  my $db   = shift()->{db} or return;
  return tied(%$db);
}

sub iddb {
  my $iddb   = shift()->{iddb} or return;
  return tied(%$iddb);
}

sub dsn {
  my $self = shift;
  my $d    = $self->{dsn};
  $self->{dsn} = shift if @_;
  $d;
}

sub load_gff_line {

  my ($self, $feat) = @_;

  #Continue if this is the '##sequence-region' header line.
  # return if lc($feat->{source}) eq "reference";

  $feat->{strand} = '' if $feat->{strand} && $feat->{strand} eq '.';
  $feat->{phase} = ''  if $feat->{phase}  && $feat->{phase}  eq '.';

  my $start = $feat->{start};
  my $stop = $feat->{stop};
  my $type = join(':',$feat->{method},$feat->{source});

  my $bin =  bin($feat->{start},$feat->{stop},MIN_BIN);
  $feat->{bin} = $bin;

  my $id = $feat->{feature_id} = $self->_next_id;
  $bin = $self->normalizeNumber($bin);

  my $db = $self->{db};
  $db->{"__class__". $feat->{gclass}}   = $id;
  $db->{"__name__".(lc $feat->{gname})} = $id;

  $db->{"__bin__"."$feat->{ref}$;$bin"}   = $id;
  $db->{"__type__".$type}                 = $id;

  for my $attr (@{$feat->{attributes}}) {
    my ($attr_name,$attr_value) = @$attr;
    $db->{"__attr__".$attr_name."__".$attr_value} = $id;
  }

  #warn "Storing start $start, stop $stop, bin $bin, id $id";
  $self->{iddb}{$id} = freeze($feat);
}

sub finish_load {
  my $self = shift;
  $self->iddb->sync;
  $self->db->sync;
  1;
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

sub get_types {
  my $self = shift;
  my ($srcseq,$class,$start,$stop,$want_count,$typelist) = @_;
  my (%obj,%result,$key,$value);
  $key = "__type__";

  if (!$srcseq) { # optimized full type list
    my $db = $self->db;
    my $status = $db->seq($key,$value,R_CURSOR);

    while ($status == 0 && $key =~ /^__type__(.+)/) {
      my $type = $1;
      my ($method,$source) = split ':',$type;
      $obj{$type} = Bio::DB::GFF::Typename->new($method,$source);
      $result{$type}++;

      if ($want_count) {
	$status = $db->seq($key,$value,R_NEXT);
      } else { # skip to next key set
	$key .= "\0";
	$status = $db->seq($key,$value,R_CURSOR)
      }

    }
  }

  else { # range search
    for my $feature (@{$self->_get_features_by_search_options(undef,
							      {rangetype => 'overlaps',
							       refseq    => $srcseq,
							       refclass  => ($class || undef),
							       start     => ($start || undef),
							       stop      => ($stop  || undef),
							      },
							      {}
							     )}
		    ) {
      my $type = Bio::DB::GFF::Typename->new($feature->{method},$feature->{source});
      $obj{$type} = $type;
      $result{$type}++;
    }
  }

  return $want_count ? %result : values %obj;
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
    if ($name =~ /^([^\*]+)\*$/) {
      $name = $1;
      $use_glob++;
    }
	
    else {
      $name = quotemeta($name);
      $name =~ s/\\\*/.*/g;
      $name =~ s/\\\?/.?/g;
      $use_regexp++;
    }
  }

  my @features;
  if ($use_glob) {
    my $callback = sub {my $feat = shift; $feat->{gname} =~ /^$name/i};
    @features = @{$self->retrieve_features_range (-table => "name",
						  -start => $name,
						  -do_while => $callback)
		};
  }
  elsif ($use_regexp) {
    my $filter = sub {my $feat = shift; $feat->{gname} =~ /$name/i};
    @features = @{$self->filter_features(-table => "name", -filter => $filter)};
  }

  else {
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
	
    if ($matches) {
      push @matches, $matches;
      return -1 if @matches >= $limit;
    }
    return $matches;
  };

  my @features = @{$self->filter_features(-table => "attrib__Note", -filter => $filter)};

  for (my $i=0; $i<scalar @matches; $i++)  {
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

sub _get_features_by_search_options {

  #The $data argument is not used and is preserved for superclass compatibility
  my ($self, $data,$search,$options) = @_;
  my $count = 0;

  my ($rangetype,$refseq,$class,$start,$stop,$types,$sparse,$order_by_group,$attributes) = 
    (@{$search}{qw(rangetype refseq refclass start stop types)},
     @{$options}{qw(sparse sort_by_group ATTRIBUTES)}) ;
	
  $start = 0               unless defined($start);
  $stop  = MAX_BIN         unless defined($stop);
	
  my %found;

  my $bin =  bin($start,$stop,MIN_BIN);  
  $bin = $self->normalizeNumber($bin);

  my @features;
  my $filter = sub {
    my $feature = shift;
	
    my $ref           = $feature->{ref};
    my $feature_start = $feature->{start};
    my $feature_stop  = $feature->{stop};
    my $feature_id    = $feature->{feature_id};

    return if $found{$feature_id};

    if (defined $refseq) {
      next unless lc $refseq eq lc $ref;
    }

    if (defined $start or defined $stop) {
      $start = 0               unless defined($start);
      $stop  = MAX_SEGMENT     unless defined($stop);

      if ($rangetype eq 'overlaps') {
	return unless $feature_stop >= $start && $feature_start <= $stop;
      } elsif ($rangetype eq 'contains') {
	return unless $feature_start >= $start && $feature_stop <= $stop;
      } elsif ($rangetype eq 'contained_in') {
	return unless $feature_start <= $start && $feature_stop >= $stop;
      } else {
	return unless $feature_start == $start && $feature_stop == $stop;
      }
    }

    my $feature_source = $feature->{source};
    my $feature_method = $feature->{method};

    if (defined $types && @$types){
      return unless $self->_matching_typelist($feature_method,$feature_source,$types);
    }

    my $feature_attributes = $feature->{attributes};
    if (defined $attributes){
      return unless $self->_matching_attributes($feature_attributes,$attributes);
    }

    $found{$feature_id}++;
    return 1;
  };

  if (defined $refseq && !$sparse) {
    my $tier = MAX_BIN;
    while ($tier >= MIN_BIN) {
      my ($tier_start,$tier_stop) = (bin_bot($tier,$start),bin_top($tier,$stop));
      # warn "Using $tier_start $tier_stop\n";
      if ($tier_start == $tier_stop) {
	$self->retrieve_features(-table => "bin",
				 -key => "$refseq$;$tier_start",
				 -filter => $filter,
				 -result => \@features);
      } else {
	my $callback = sub {my $feat = shift; $feat->{bin} <= $tier_stop};
	$self->retrieve_features_range(-table => "bin",
				       -start => "$refseq$;$tier_start",
				       -do_while => $callback,
				       -filter => $filter,
				       -result => \@features);
      }
	
      $tier /= 10;
    }
  }

  elsif (@$types) {
    foreach (@$types) {
      my $type = join ':',@$_;
      $self->retrieve_features_range(-table    => 'type',
				     -start    => $type,
				     -do_while => sub { my $f = shift;
							$self->_matching_typelist(
										  $f->{method},
										  $f->{source},
										  [$_]) },
				     -result => \@features);
    }
  }

  elsif (defined $attributes) {
    my ($attribute_name,$attribute_value) = each %$attributes; # pick first one
    $self->retrieve_features(-table => 'attr',
			     -key   => "${attribute_name}__${attribute_value}",
			     -filter => $filter,
			     -result  => \@features);
  }

  else { # linear search
    while (my ($key,$value) = each %{$self->{iddb}}) {
      my $feature = thaw($self->{iddb}{$key});
      next unless $filter->($feature);
      push @features,$feature;
    }
  }

  return \@features;
}

sub retrieve_features {
  my $self = shift;
  my ($table, $key, $filter, $result) = rearrange(['TABLE','KEY','FILTER', 'RESULT'],@_);

  my @result;
  $result ||= \@result;

  my $frozen;
  my @ids = $self->db->get_dup("__".$table."__".$key);
  my $iddb = $self->{iddb};

  foreach my $id (@ids) {
    my $frozen = $iddb->{$id};
    my $feat = thaw($frozen);
    next if $filter && !$filter->($feat);
    push @$result, $feat;
  }
  return $result;
}

sub retrieve_features_range {
  my ($self) = shift;
  my ($table, $start, $do_while, $filter, $result) = rearrange(['TABLE','START','DO_WHILE', 'FILTER', 'RESULT'],@_);

  my @result;
  $result ||= \@result;
  my ($id, $key, $value);

  $key = "__".$table."__".$start;
  my $db   = $self->db;

  for (my $status = $db->seq($key,$value,R_CURSOR);
       $status == 0;
       $status = $db->seq($key,$value,R_NEXT)) {

    my $feature = thaw($self->{iddb}{$value});
    last unless $do_while->($feature);

    unless ($filter) {
      push @$result, $feature;
    } else {
      my $filter_result = $filter->($feature);
      push @$result, $feature if $filter_result;
      last if $filter_result == -1;
    }
  }

  return $result;
}

sub filter_features {
  my ($self) = shift;

  my ($table, $filter) = rearrange(['TABLE','FILTER'],@_);

  my @result;
  my ($key, $frozen);

  my $iddb = $self->iddb;
  for (my $status = $iddb->seq($key,$frozen,R_FIRST);
       $status == 0;
       $status = $iddb->seq($key,$frozen,R_NEXT)) {

    my $feature = thaw($frozen);
    if (my $filter_result = $filter->($feature)) {
      push @result, $feature;
      last if $filter_result == -1;
    }
  }

  return \@result;
}


sub _basic_features_by_id{
  my $self = shift;
  my ($ids) = @_;

  $ids = [$ids] unless ref $ids =~ /ARRAY/;

  my @result;
  my $iddb = $self->{iddb};
  for my $feature_id (@$ids){
    my $frozen = $iddb->{$feature_id};
    push @result, thaw($frozen);
  }

  return wantarray() ? @result : $result[0];
}

sub normalizeNumber {
  my ($self, $num) = @_;
  while ((length $num) < MAX_NUM_LENGTH)
  {
    $num = "0".$num;
  }
  return $num;
}

1;

