package Bio::DB::GFF::Adaptor::berkeleydb;

use strict;

use Bio::DB::GFF::Util::Rearrange; # for rearrange()
use Bio::DB::GFF::Util::Binning;
use Bio::DB::Fasta;
use Bio::DB::GFF::Adaptor::berkeleydb::iterator;
use Bio::DB::GFF::Adaptor::memory::feature_serializer; # qw(feature2string string2feature @hash2array_map);

use DB_File;
use File::Path 'mkpath';

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
  my ($dbdir,$preferred_groups,$autoindex) = rearrange([
							[qw(DSN DB)],
							'PREFERRED_GROUPS',
							[qw(DIR AUTOINDEX)],
					    ],@_);
  if (defined $dbdir && defined $autoindex) {
    $class->throw("If both -dsn and -dir (or -autoindex) are specified, they must point to the same directory")
      unless $dbdir eq $autoindex;
  }

  $dbdir ||= $autoindex;
  $dbdir ||= $ENV{TMPDIR} ? "$ENV{TMPDIR}/test" : "/tmp/test";

  my $self = bless {},$class;
  $self->dsn($dbdir);
  $self->preferred_groups($preferred_groups) if defined $preferred_groups;
  $self->_autoindex                          if $autoindex;
  $self->_open_databases();
  return $self;
}

sub _autoindex {
  my $self = shift;
  my $dir    = $self->dsn;
  my %ignore = map {$_=>1} ($self->_index_file,$self->_hash_file,$self->_fasta_file,$self->_temp_file,$self->_timestamp_file);

  my $maxtime   = 0;
  my $maxfatime = 0;

  opendir (D,$dir) or $self->throw("Couldn't open directory $dir for reading: $!");

  while (defined (my $node = readdir(D))) {
    next if $node =~ /^\./;
    my $path      = "$dir/$node";
    next if $ignore{$path};
    next unless -f $path;
    my $mtime = _mtime(\*_);  # not a typo
    $maxtime   = $mtime if $mtime > $maxtime;
    $maxfatime = $mtime if $mtime > $maxfatime && $node =~ /\.(?:fa|fasta|dna)(?:\.gz)?$/;
  }

  close D;

  my $timestamp_time  = _mtime($self->_timestamp_file) || 0;
  my $all_files_exist = -e $self->_index_file && -e $self->_hash_file && (-e $self->_fasta_file || !$maxfatime);

  # to avoid rebuilding FASTA files if not changed
  my $spare_fasta     = $maxfatime > 0 && $maxfatime < $timestamp_time && -e $self->_fasta_file;  

  if ($maxtime > $timestamp_time || !$all_files_exist) {
    print STDERR __PACKAGE__,": Reindexing files in $dir. This may take a while....\n";
    $self->do_initialize(1,$spare_fasta);
    $self->load_gff($dir);
    $self->load_fasta($dir) unless $spare_fasta;
    print STDERR __PACKAGE__,": Reindexing done\n";
  }

  else {
    $self->_open_databases();
  }

}

sub _open_databases {
  my $self = shift;
  my $dsn  = $self->dsn;
  unless (-d $dsn) {  # try to create the directory
    mkpath($dsn) or die "Couldn't create database directory $dsn: $!";
  }
  my %db;
  local $DB_BTREE->{flags} = R_DUP;
  tie(%db,'DB_File',$self->_index_file,O_RDWR|O_CREAT,0666,$DB_BTREE)
    or $self->throw("Couldn't tie ".$self->_index_file.": $!");
  $self->{db}   = \%db;
  $self->{data} = FeatureStore->new($self->_data_file);

  if (-e $self->_fasta_file) {
    my $dna_db = Bio::DB::Fasta->new($self->_fasta_file) or $self->throw("Can't reindex sequence file: $@");
    $self->dna_db($dna_db);
  }

}

sub _close_databases {
  my $self = shift;
  delete $self->{db};
  delete $self->{data};
}

# do nothing!
sub setup_load { 1; }

sub _delete_features {
  my $self        = shift;
  my @feature_ids = @_;
  my $removed = 0;
  my $last_id = $self->{data}->last_id;
  for my $id (@feature_ids) {
    next unless $id >= 0 && $id < $last_id;
    my $feat  = $self->{data}->get($id) or next;
    $self->{data}->remove($id);
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
		"__name__".lc(join ":",$feat->{gclass},$feat->{gname}),
		"__bin__".lc("$feat->{ref}$;$feat->{bin}"),
		"__type__".lc(join(':',$feat->{method},$feat->{source})),
		map {"__attr__".lc(join(':',$_->[0],$_->[1]))} @{$feat->{attributes}}
	  );
}

sub _delete {
  my $self        = shift;
  my $delete_spec = shift;
  return $self->SUPER::_delete($delete_spec) if @{$delete_spec->{segments}} or @{$delete_spec->{types}};
  $self->throw("This operation would delete all feature data and -force not specified")
    unless $delete_spec->{force};
  my $deleted = $self->{db}{__count__};
  $self->{data} = FeatureStore->new($self->_data_file,1);
  %{$self->{db}}   = ();
  $deleted;
}

# with duplicates enabled, we cannot simply do $db->{__index__}++;
sub _bump_feature_count {
  my $self = shift;
  my $db = $self->{db};
  if (@_) {
    delete $db->{__count__};
    return $db->{__tount__} = shift;
  } else {
    my $index = ${db}->{__count__};
    delete $db->{__count__};
    $db->{__count__} = $index + 1;
    return $index;
  }
}

sub do_initialize {
  my $self  = shift;
  my $erase = shift;
  my $spare_fasta = shift;
  if ($erase) {
    $self->_close_databases;
    unlink $self->_index_file;
    unlink $self->_data_file;
    unless ($spare_fasta) {
      unlink $self->_fasta_file;
      unlink $self->_fasta_file.'.index';
    }
    unlink $self->_timestamp_file;
    $self->_open_databases;
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
  $self->_touch_timestamp;
  return $loaded;
}

sub _mtime {
  my $file = shift;
  my @stat = stat($file);
  return $stat[9];
}

sub _index_file {
  my $self = shift;
  return $self->dsn . "/_bdb_features.btree";
}

sub _data_file {
  my $self = shift;
  return $self->dsn . "/_bdb_features.data";
}

sub _fasta_file {
  my $self = shift;
  return $self->dsn . "/_bdb_sequence.fa";
}

sub _temp_file {
  my $self = shift;
  return $self->dsn ."/_bdb_temporary_results.btree";
}

sub _timestamp_file {
  my $self = shift;
  return $self->dsn ."/_bdb_timestamp";
}

sub db {
  my $db   = shift()->{db} or return;
  return tied(%$db);
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

  my $id = $self->{data}->put($feat);
  $bin = $self->normalizeNumber($bin);

  my $db = $self->{db};
  for my $skey ($self->_secondary_keys($feat)) {
    $db->{$skey} = $id;
  }

  $self->_bump_feature_count();

}

sub finish_load {
  my $self = shift;
  $self->db->sync;
  $self->_touch_timestamp;
  1;
}

sub _touch_timestamp {
  my $self = shift;
  my $tsf = $self->_timestamp_file;
  open (F,">$tsf") or $self->throw("Couldn't open $tsf: $!");
  print F scalar(localtime);
  close F;
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
  my @features = @{$self->retrieve_features(-table => 'name', -key=>"$class:$name")};
  if (!@features) {  # nothing matched exactly, so try aliases
    @features = @{$self->retrieve_features(-table=>'attr',-key=>"Alias:$name")};
  }

  foreach my $feature (@features){
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
    for my $feature (@{$self->_get_features_by_search_options(
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
    @features = @{$self->retrieve_features_range (-table => 'name',
						  -start => "$class:$name",
						  -do_while => $callback)
		};
  }
  elsif ($use_regexp) {
    my $filter = sub {my $feat = shift; $feat->{gname} =~ /$name/i};
    @features = @{$self->filter_features(-table =>'name', -filter => $filter)};
  }

  else {
    @features = @{$self->retrieve_features(-table => 'name', -key => "$class:$name")};
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

    $callback->(@{$feature}{@hash2array_map},0);
  }
  return $count;
}

#sub get_feature_by_attribute{
sub _feature_by_attribute{
  my $self = shift;
  my ($attributes,$callback) = @_;
  $callback || $self->throw('must provide a callback argument');
  my $count = 0;
  my $feature_group_id = undef;

  #there could be more than one set of attributes......
  while (my ($key, $value) = each %$attributes) {
	
    my @features = @{$self->retrieve_features
		       (-table => "attr", -key => "$key:$value")};

    for my $feature (@features) {
      $callback->(@{$feature}{@hash2array_map},$feature_group_id);
      $count++;
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
    return unless defined $feature->{gclass} && defined $feature->{gname}; # ignore NULL objects
    return unless $feature->{attributes};
    my @attributes = @{$feature->{attributes}};
    my @values     = map {$_->[1]} @attributes;
    my $value      = "@values";

    my @hits;
    while ($value =~ /($search)/ig) {
      push @hits,$1;
    }
	
    if (@hits) {
      push @matches, scalar @hits;
      return -1 if @matches >= $limit;
    }
    return scalar @hits;
  };

  my @features = @{$self->filter_features(-table => "attrib__note", -filter => $filter)};

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
  my ($self, $search,$options) = @_;
  my $count = 0;

  my ($rangetype,$refseq,$class,$start,$stop,$types,$sparse,$order_by_group,$attributes,$temp_table) =
    (@{$search}{qw(rangetype refseq refclass start stop types)},
     @{$options}{qw(sparse sort_by_group ATTRIBUTES temp_table)}) ;
	
  $start = 0               unless defined($start);
  $stop  = MAX_BIN         unless defined($stop);
	
  my $bin =  bin($start,$stop,MIN_BIN);  
  $bin = $self->normalizeNumber($bin);

  my ($results,@features,%found,%results_table);

  if ($temp_table) {
    local $DB_BTREE->{flags} = R_DUP;
    unlink $self->_temp_file;
    tie(%results_table,'DB_File',$self->_temp_file,O_RDWR|O_CREAT,0666,$DB_BTREE) 
      or $self->throw("Couldn't tie ".$self->_temp_file.": $!");
    $results = \%results_table;
  } else {
    $results = \@features;
  }

  my $filter = sub {
    my $feature = shift;
	
    my $ref           = $feature->{ref};
    my $feature_start = $feature->{start};
    my $feature_stop  = $feature->{stop};
    my $feature_id    = $feature->{feature_id};

    return 0 if $found{$feature_id}++;

    if (defined $refseq) {
      return 0 unless lc $refseq eq lc $ref;
      $start = 0               unless defined($start);
      $stop  = MAX_SEGMENT     unless defined($stop);

      if ($rangetype eq 'overlaps') {
	return 0 unless $feature_stop >= $start && $feature_start <= $stop;
      } elsif ($rangetype eq 'contains') {
	return 0 unless $feature_start >= $start && $feature_stop <= $stop;
      } elsif ($rangetype eq 'contained_in') {
	return 0 unless $feature_start <= $start && $feature_stop >= $stop;
      } else {
	return 0 unless $feature_start == $start && $feature_stop == $stop;
      }
    }

    my $feature_source = $feature->{source};
    my $feature_method = $feature->{method};

    if (defined $types && @$types){
      return 0 unless $self->_matching_typelist($feature_method,$feature_source,$types);
    }

    my $feature_attributes = $feature->{attributes};
    if (defined $attributes){
      return 0 unless $self->_matching_attributes($feature_attributes,$attributes);
    }

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
				 -result => $results);
      } else {
	my $callback = sub {my $feat = shift; $feat->{bin} <= $tier_stop};
	$self->retrieve_features_range(-table => "bin",
				       -start => "$refseq$;$tier_start",
				       -do_while => $callback,
				       -filter => $filter,
				       -result => $results);
      }
	
      $tier /= 10;
    }
  }

  elsif (@$types) {
    foreach (@$types) {
      my $type = join ':',@$_;
      $self->retrieve_features_range(-table    => 'type',
				     -start    => $type,
				     -filter   => $filter,
				     -do_while => sub { my $f = shift;
							lc($f->{method}) eq lc($_->[0]) 
							  &&
							lc($f->{source}||$_->[1]||'') eq lc($_->[1]||'')
						      },
				     -result => $results);
    }
  }

  elsif (defined $attributes) {
    my ($attribute_name,$attribute_value) = each %$attributes; # pick first one
    $self->retrieve_features(-table => 'attr',
			     -key   => "${attribute_name}:${attribute_value}",
			     -filter => $filter,
			     -result  => $results);
  }

  else {
    $self->filter_features(-filter => $filter,-result=>$results);
  }

  return $results;
}

sub retrieve_features {
  my $self = shift;
  my ($table, $key, $filter, $result) = rearrange(['TABLE','KEY','FILTER', 'RESULT'],@_);

  my @result;
  $result ||= \@result;

  my $frozen;
  my @ids  = $self->db->get_dup("__".lc($table)."__".lc($key));
  my $data = $self->{data};

  foreach my $id (@ids) {
    my $feat = $data->get($id);
    my $filter_result = $filter ? $filter->($feat) : 1;
    next unless $filter_result;
    if (ref $result eq 'HASH') {
     $result->{"$feat->{gclass}:$feat->{gname}"} = join ($;,$self->_hash_to_array($feat));
    } else {
      push @$result, $feat;
    }
    last if $filter_result == -1;
  }
  return $result;
}

sub retrieve_features_range {
  my ($self) = shift;
  my ($table, $start, $do_while, $filter, $result) = rearrange(['TABLE','START','DO_WHILE', 'FILTER', 'RESULT'],@_);

  my @result;
  $result ||= \@result;
  my ($id, $key, $value);

  $key = lc "__".$table."__".$start;
  my $db   = $self->db;

  for (my $status = $db->seq($key,$value,R_CURSOR);
       $status == 0;
       $status = $db->seq($key,$value,R_NEXT)) {

    my $feat = $self->{data}->get($value);
    last unless $do_while->($feat);

    my $filter_result = $filter ? $filter->($feat) : 1;
    next unless $filter_result;

    if (ref $result eq 'HASH') {
     $result->{"$feat->{gclass}:$feat->{gname}"} = join($;,$self->_hash_to_array($feat));
    } else {
      push @$result,$feat;
    }
    last if $filter_result == -1;
  }

  return $result;
}


sub filter_features {
  my ($self) = shift;

  my ($filter,$result) = rearrange(['FILTER','RESULT'],@_);

  my @result;
  $result ||= \@result;

  my ($key, $frozen);
  my $data = $self->{data};
  $data->reset;
  while (my $feat = $data->next) {

    my $filter_result = $filter ? $filter->($feat) : 1;
    next unless $filter_result;

    if (ref($result) eq 'HASH') {
     $result->{"$feat->{gclass}:$feat->{gname}"} = join($;,$self->_hash_to_array($feat));
    } else {
      push @$result,$feat;
    }
    last if $filter_result == -1;
  }

  return $result;
}


sub _basic_features_by_id{
  my $self = shift;
  my ($ids) = @_;

  $ids = [$ids] unless ref $ids =~ /ARRAY/;

  my @result;
  my $data = $self->{data};
  for my $feature_id (@$ids){
    push @result, $data->get($feature_id);
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

sub get_features_iterator {
  my $self = shift;

  my ($search,$options,$callback) = @_;
  $callback || $self->throw('must provide a callback argument');
  $options->{temp_table}++;

  my $results = $self->_get_features_by_search_options($search,$options);
  return Bio::DB::GFF::Adaptor::berkeleydb::iterator->new($results,$callback,$self->_temp_file);
}

#--------------------------------------------------------------------------#

package FeatureStore;

# This is a very specialized package that stores serialized features onto a file-based
# array. The array is indexed by the physical offset to the beginning of each serialized
# feature.

use strict;
use Fcntl qw(SEEK_SET SEEK_END);
use base 'Bio::Root::Root';
use Bio::DB::GFF::Adaptor::memory::feature_serializer; # qw(feature2string string2feature @hash2array_map);

sub new {
  my $class  = shift;
  my $dbname = shift    or $class->throw("must provide a filepath argument");
  my $overwrite = shift;

  my $mode = $overwrite ? "+>" : "+>>";

  open (F,$mode,$dbname) or $class->throw("$dbname: $!");
  my $self = bless {
		    fh        => \*F,
		    next_idx  => 0,
		    last_id   => 0,
		   },$class;
  return $self;
}

sub put {
  my $self   = shift;
  my $feature = shift;
  my $fh = $self->{fh};
  seek($fh,0,SEEK_END);
  my $offset = tell($fh) || 0;

  $self->{last_id} = $offset;

  my $id = pack("L",$offset);
  $feature->{feature_id} = $id;
  my $value = feature2string($feature);
  print $fh pack("n/a*",$value);


  return $id;
}

sub last_id {
  shift->{last_id};
}

sub get {
  my $self     = shift;
  my $idx      = shift;
  my $offset   = unpack("L",$idx);
  my $fh = $self->{fh};

  my ($value,$length);
  seek($fh,$offset,SEEK_SET);
  return unless read($fh,$length,2);
  return unless read($fh,$value,unpack("n",$length));
  $self->{next_idx} = tell($fh);
  return if substr($value,0,1) eq "\0";
  return string2feature($value);
}

sub next {
  my $self = shift;
  my $fh     = $self->{fh};
  my $result;
  do {
    $result = $self->get(pack("L",$self->{next_idx}));
  } until $result || eof($fh);
  $self->{next_idx} = 0 unless $result;
  $result;
}

sub remove {
  my $self   = shift;
  my $id     = shift;
  my $offset = unpack("L",$id);
  my $fh     = $self->{fh};
  my ($value,$length);
  seek($fh,$offset,SEEK_SET);
  return unless read($fh,$length,2);
  print $fh "\0"x$length;  # null it out
  1;
}

sub _seek {
  my $self = shift;
  my $idx  = shift;
  my $offset   = unpack("L",$idx);
  seek($self->{fh},$offset,SEEK_SET);
  $self->{next_idx} = tell($self->{fh});
}

sub reset {
  my $self = shift;
  $self->_seek(pack("L",0));
}

sub _feature2string {
  my $feature = shift;
  my @a = @{$feature}{@hash2array_map};
  push @a,map {@$_} @{$feature->{attributes}} if $feature->{attributes};
  return join $;,@a;
}

sub _string2feature {
  my $string  = shift;
  my (%feature,@attributes);

  (@feature{@hash2array_map},@attributes) = split $;,$string;
  while (@attributes) {
    my ($key,$value) = splice(@attributes,0,2);
    push @{$feature{attributes}},[$key,$value];
  }
  $feature{group_id} = undef;
  \%feature;
}


1;
