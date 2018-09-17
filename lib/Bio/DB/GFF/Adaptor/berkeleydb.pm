package Bio::DB::GFF::Adaptor::berkeleydb;


=head1 NAME

Bio::DB::GFF::Adaptor::berkeleydb -- Bio::DB::GFF database adaptor for in-memory databases

=head1 SYNOPSIS

  use Bio::DB::GFF;
  my $db = Bio::DB::GFF->new(-adaptor=> 'berkeleydb',
                             -create => 1, # on initial build you need this
			     -dsn    => '/usr/local/share/gff/dmel');

  # initialize an empty database, then load GFF and FASTA files
  $db->initialize(1);
  $db->load_gff('/home/drosophila_R3.2.gff');
  $db->load_fasta('/home/drosophila_R3.2.fa');

  # do queries
  my $segment  = $db->segment(Chromosome => '1R');
  my $subseg   = $segment->subseq(5000,6000);
  my @features = $subseg->features('gene');

See L<Bio::DB::GFF> for other methods.

=head1 DESCRIPTION

This adaptor implements a berkeleydb-indexed version of Bio::DB::GFF.
It requires the DB_File and Storable modules. It can be used to store
and retrieve short to medium-length GFF files of several million
features in length.

=head1 CONSTRUCTOR

Use Bio::DB::GFF-E<gt>new() to construct new instances of this class.
Three named arguments are recommended:

 Argument    Description
 --------    -----------

 -adaptor    Set to "berkeleydb" to create an instance of this class.

 -dsn        Path to directory where the database index files will be stored (alias -db)

 -autoindex  Monitor the indicated directory path for FASTA and GFF files, and update the
               indexes automatically if they change (alias -dir)

 -write      Set to a true value in order to update the database.

 -create     Set to a true value to create the database the first time
               (implies -write)

 -tmp        Location of temporary directory for storing intermediate files
               during certain queries.

 -preferred_groups  Specify the grouping tag. See L<Bio::DB::GFF>

The -dsn argument selects the directory in which to store the database
index files. If the directory does not exist it will be created
automatically, provided that the current process has sufficient
privileges. If no -dsn argument is specified, a database named "test"
will be created in your system's temporary files directory.

The -tmp argument specifies the temporary directory to use for storing
intermediate search results. If not specified, your system's temporary
files directory will be used. On Unix systems, the TMPDIR environment
variable is honored. Note that some queries can require a lot of
space.

The -autoindex argument, if present, selects a directory to be
monitored for GFF and FASTA files (which can be compressed with the
gzip program if desired). Whenever any file in this directory is
changed, the index files will be updated. Note that the indexing can
take a long time to run: anywhere from 5 to 10 minutes for a million
features. An alias for this argument is -dir, which gives this adaptor
a similar flavor to the "memory" adaptor.

-dsn and -dir can point to the same directory. If -dir is given but
-dsn is absent the index files will be stored into the directory
containing the source files.  For autoindexing to work, you must
specify the same -dir path each time you open the database.

If you do not choose autoindexing, then you will want to load the
database using the bp_load_gff.pl command-line tool. For example:

 bp_load_gff.pl -a berkeleydb -c -d /usr/local/share/gff/dmel dna1.fa dna2.fa features.gff

=head1 METHODS

See L<Bio::DB::GFF> for inherited methods

=head1 BUGS

The various get_Stream_* methods and the features() method with the
-iterator argument only return an iterator after the query runs
completely and the module has been able to generate a temporary
results file on disk. This means that iteration is not as big a win as
it is for the relational-database adaptors.

Like the dbi::mysqlopt adaptor, this module uses a binning scheme to
speed up range-based searches. The binning scheme used here imposes a
hard-coded 1 gigabase (1000 Mbase) limit on the size of the largest
chromosome or other reference sequence.

=head1 SEE ALSO

L<Bio::DB::GFF>, L<bioperl>

=head1 AUTHORS

Vsevolod (Simon) Ilyushchenko E<gt>simonf@cshl.eduE<lt>
Lincoln Stein E<gt>lstein@cshl.eduE<lt>

Copyright (c) 2005 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

use strict;

use DB_File;
use File::Path 'mkpath';
use File::Spec;
use File::Temp 'tempfile';

use Bio::DB::GFF::Util::Rearrange; # for rearrange()
use Bio::DB::GFF::Util::Binning;
use Bio::DB::Fasta;
use Bio::DB::GFF::Adaptor::berkeleydb::iterator;
use Bio::DB::GFF::Adaptor::memory::feature_serializer; # qw(feature2string string2feature @hash2array_map);

# this is the smallest bin (1 K)
use constant MIN_BIN    => 1000;
# this is the largest that any reference sequence can be (1000 megabases)
use constant MAX_BIN     => 1_000_000_000;
use constant MAX_SEGMENT => 1_000_000_000;  # the largest a segment can get

#We have to define a limit because Berkeleydb sorts in lexicografic order,
#so all the numbers have to have the same length.	
use constant MAX_NUM_LENGTH => length(MAX_BIN);

use base 'Bio::DB::GFF::Adaptor::memory';

sub new {
  my $class = shift ;
  my ($dbdir,$preferred_groups,$autoindex,$write,$create,$tmpdir) = rearrange([
									       [qw(DSN DB)],
									       'PREFERRED_GROUPS',
									       [qw(DIR AUTOINDEX)],
									       [qw(WRITE WRITABLE)],
									       'CREATE',
									       'TMP',
									      ],@_);
  $tmpdir ||= File::Spec->tmpdir;
  $dbdir  ||= $autoindex;
  $dbdir  ||= "$tmpdir/test";
  $write  ||= $create;

  my $self = bless {},$class;
  $self->dsn($dbdir);
  $self->tmpdir($tmpdir);
  $self->preferred_groups($preferred_groups) if defined $preferred_groups;
  $self->_autoindex($autoindex)              if $autoindex;
  $self->_open_databases($write,$create);
  return $self;
}

sub _autoindex {
  my $self    = shift;
  my $autodir = shift;

  my $dir    = $self->dsn;
  my %ignore = map {$_=>1} ($self->_index_file,$self->_data_file,
			    $self->_fasta_file,$self->_temp_file,
			    $self->_notes_file,
			    $self->_timestamp_file);

  my $maxtime   = 0;
  my $maxfatime = 0;

  opendir (my $D,$autodir) or $self->throw("Couldn't open directory $autodir for reading: $!");

  while (defined (my $node = readdir($D))) {
    next if $node =~ /^\./;
    my $path      = "$dir/$node";
    next if $ignore{$path};
    next unless -f $path;
    my $mtime = _mtime(\*_);  # not a typo
    $maxtime   = $mtime if $mtime > $maxtime;
    $maxfatime = $mtime if $mtime > $maxfatime && $node =~ /\.(?:fa|fasta|dna)(?:\.gz)?$/;
  }

  close $D;

  my $timestamp_time  = _mtime($self->_timestamp_file) || 0;
  my $all_files_exist = -e $self->_index_file && -e $self->_data_file && (-e $self->_fasta_file || !$maxfatime);

  # to avoid rebuilding FASTA files if not changed
  my $spare_fasta     = $maxfatime > 0 && $maxfatime < $timestamp_time && -e $self->_fasta_file;  

  if ($maxtime > $timestamp_time || !$all_files_exist) {
    print STDERR __PACKAGE__,": Reindexing files in $dir. This may take a while....\n";
    $self->do_initialize(1,$spare_fasta);
    $self->load_gff($autodir,1);
    $self->load_fasta($autodir,1) unless $spare_fasta;
    print STDERR __PACKAGE__,": Reindexing done\n";
  }

  else {
    $self->_open_databases();
  }

}

sub _open_databases {
  my $self   = shift;
  my ($write,$create) = @_;

  my $dsn  = $self->dsn;
  unless (-d $dsn) {  # directory does not exist
    $create or $self->throw("Directory $dsn does not exist and you did not specify the -create flag");
    mkpath($dsn) or $self->throw("Couldn't create database directory $dsn: $!");
  }

  my %db;
  local $DB_BTREE->{flags} = R_DUP;
  $DB_BTREE->{compare}     = sub { lc($_[0]) cmp lc($_[1]) };
  my $flags = O_RDONLY;
  $flags   |= O_RDWR  if $write;
  $flags   |= O_CREAT if $create;

  tie(%db,'DB_File',$self->_index_file,$flags,0666,$DB_BTREE)
    or $self->throw("Couldn't tie ".$self->_index_file.": $!");
  $self->{db}   = \%db;
  $self->{data} = FeatureStore->new($self->_data_file,$write,$create);

  if (-e $self->_fasta_file) {
    my $dna_db = Bio::DB::Fasta->new($self->_fasta_file) or $self->throw("Can't reindex sequence file: $@");
    $self->dna_db($dna_db);
  }

  my $mode =  $write  ? "+>>"
            : $create ? "+>"
            : "<";

  my $notes_file = $self->_notes_file;
  open my $F, $mode, $notes_file or $self->throw("Could not open file '$notes_file': $!");
  $self->{notes} = $F;
}

sub _close_databases {
  my $self = shift;
  delete $self->{db};
  delete $self->{data};
  delete $self->{notes};
}

sub _delete_features {
  my $self        = shift;
  my @feature_ids = @_;
  my $removed = 0;
  my $last_id = $self->{data}->last_id;
  for my $id (@feature_ids) {
    next unless $id >= 0 && $id < $last_id;
    my $feat  = $self->{data}->get($id) or next;
    $self->{data}->remove($id);
    $self->_bump_class_count($feat->{gclass},-1);
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
		"__type__".join(':',$feat->{method},$feat->{source}),
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
  $self->{data} = FeatureStore->new($self->_data_file,1,1);
  %{$self->{db}}   = ();
  $deleted;
}

# with duplicates enabled, we cannot simply do $db->{__index__}++;
sub _bump_feature_count {
  my $self = shift;
  my $db = $self->{db};
  if (@_) {
    delete $db->{__count__};
    return $db->{__count__} = shift;
  } else {
    my $index = ${db}->{__count__};
    delete $db->{__count__};
    $db->{__count__} = ($index || 0) + 1;
    return $index;
  }
}

sub _bump_class_count {
  my $self = shift;
  my ($class,$count) = @_;
  $count ||= 1;
  my $db  = $self->{db};
  my $key = "__class__$class";
  my $newcount = ($db->{$key} || 0) + $count;
  delete $db->{$key};
  $db->{$key} = $newcount;
}

sub classes {
  my $self = shift;
  my $db   = $self->db;
  my ($key,$value) = ('__class__',undef);
  my %classes;
  for (my $status = $db->seq($key,$value,R_CURSOR);
       $status == 0;
       $status = $db->seq($key,$value,R_NEXT)) {
    my ($class) = $key =~ /^__class__(.+)/ or last;
    $classes{$class}++ if $value > 0;
  }
  my @classes = sort keys %classes;
  return @classes;
}

sub do_initialize {
  my $self  = shift;
  my $erase = shift;
  my $spare_fasta = shift; # used internally only!
  if ($erase) {
    $self->_close_databases;
    unlink $self->_index_file;
    unlink $self->_data_file;
    unlink $self->_notes_file;
    unless ($spare_fasta) {
      unlink $self->_fasta_file;
      unlink $self->_fasta_file.'.index';
    }
    unlink $self->_timestamp_file;
    $self->_open_databases(1,1);
  }
  1;
}

# load_sequence($fasta_filehandle,$first_sequence_id)
sub load_sequence {
  my $self = shift;
  my ($io_handle,$id) = @_;
  my $file = $self->_fasta_file;
  my $loaded = 0;

  open my $F, '>>', $file or $self->throw("Could not append file '$file': $!");

  if (defined $id) {
    print $F ">$id\n";
    $loaded++;
  }

  while (<$io_handle>) {
    $loaded++ if /^>/;
    print $F $_;
  }
  close $F;
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
  return $self->dsn . "/bdb_features.btree";
}

sub _data_file {
  my $self = shift;
  return $self->dsn . "/bdb_features.data";
}

sub _fasta_file {
  my $self = shift;
  return $self->dsn . "/bdb_sequence.fa";
}

sub _notes_file {
  my $self = shift;
  return $self->dsn . "/bdb_notes.idx";
}

sub _temp_file {
  my $self = shift;
  local $^W=0;
  my (undef,$filename) = tempfile("bdb_temp_XXXXXX",DIR=>$self->tmpdir,OPEN=>0);
  return $filename;
}

sub _timestamp_file {
  my $self = shift;
  return $self->dsn ."/bdb_timestamp";
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

sub tmpdir {
  my $self = shift;
  my $d    = $self->{tmpdir};
  $self->{tmpdir} = shift if @_;
  $d;
}

sub load_gff_line {

  my ($self, $feat) = @_;

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

  # save searchable notes to separate index
  my $fh = $self->{notes};
  my @notes = map {$_->[1]} grep {lc $_->[0] eq 'note'} @{$feat->{attributes}};
  print $fh $_,"\t",pack("u*",$id) or $self->throw("An error occurred while updating indexes: $!")
    foreach @notes;

  $self->{records_loaded}++;
  $self->_bump_feature_count();
  $self->_bump_class_count($feat->{gclass});

}

# do nothing!
sub setup_load {
  my $self = shift;
  $self->{records_loaded} = 0;
  1;
}

sub finish_load {
  my $self = shift;
  $self->db->sync && $self->throw("An error occurred while updating indexes: $!");
  $self->_touch_timestamp;
  $self->{records_loaded};
}

sub _touch_timestamp {
  my $self = shift;
  my $tsf = $self->_timestamp_file;
  open my $F, '>', $tsf or $self->throw("Could not write file '$tsf': $!");
  print $F scalar(localtime);
  close $F;
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
  my ($use_regexp, $use_glob,$using_alias_search);

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
    @features = @{$self->retrieve_features(-table=>'name',   -key => "$class:$name")};
  }

  unless (@features) {
    $using_alias_search++;
    @features = @{$self->retrieve_features(-table=>'attr',   -key=>"Alias:$name")};
  }

  foreach my $feature (@features){
    $id++;
    next unless $using_alias_search || $feature->{gclass} eq $class;

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

  $search_string =~ tr/*?//d;

  my @results;

  my @words = map {quotemeta($_)} $search_string =~ /(\w+)/g;
  my $search = join '|',@words;

  my (%found,$found);
  my $note_index = $self->{notes};
  seek($note_index,0,0);  # back to start
  while (<$note_index>) {
    next unless /$search/;
    chomp;
    my ($note,$uu) = split "\t";
    $found{unpack("u*",$uu)}++;
    last if $limit && ++$found >= $limit;
  }

  my (@features, @matches);
  for my $idx (keys %found) {
    my $feature = $self->{data}->get($idx) or next;
    my @attributes = @{$feature->{attributes}};
    my @values     = map {lc $_->[0] eq 'note' ? $_->[1] : ()} @attributes;
    my $value      = "@values";

    my $hits;
    $hits++ while $value =~ /($search)/ig;  # count the number of times we were hit
    push @matches,$hits;
    push @features,$feature;
  }

  for (my $i=0; $i<@matches; $i++)  {
    my $feature = $features[$i];
    my $matches = $matches[$i];

    my $relevance = 10 * $matches;
    my $featname = Bio::DB::GFF::Featname->new($feature->{gclass}=>$feature->{gname});
    my $type     = Bio::DB::GFF::Typename->new($feature->{method}=>$feature->{source});
    my $note;
    $note   = join ' ',map {$_->[1]} grep {$_->[0] eq 'Note'}                @{$feature->{attributes}};
    push @results,[$featname,$note,$relevance,$type];
  }

  return @results;
}

sub _get_features_by_search_options {

  #The $data argument is not used and is preserved for superclass compatibility
  my ($self, $search,$options) = @_;
  my $count = 0;

  my ($rangetype,$refseq,$class,$start,$stop,$types,$sparse,$order_by_group,$attributes,$temp_file) =
    (@{$search}{qw(rangetype refseq refclass start stop types)},
     @{$options}{qw(sparse sort_by_group ATTRIBUTES temp_file)}) ;

  $start = 0               unless defined($start);
  $stop  = MAX_BIN         unless defined($stop);

  my $bin =  bin($start,$stop,MIN_BIN);  
  $bin = $self->normalizeNumber($bin);

  my ($results,@features,%found,%results_table);

  if ($temp_file) {
    local $DB_BTREE->{flags} = R_DUP;
    # note: there is a race condition possible here, if someone reuses the
    # same name between the time we get the tmpfile name and the time we
    # ask DB_File to open it.
    tie(%results_table,'DB_File',$temp_file,O_RDWR|O_CREAT,0666,$DB_BTREE)
      or $self->throw("Couldn't tie temporary file ".$temp_file." for writing: $!");
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
  local $^W = 0;   # because _hash_to_array() will generate lots of uninit values

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
  local $^W = 0;  # because _hash_to_array will generate lots of uninit warnings

  my @result;
  $result ||= \@result;
  my ($id, $key, $value);

  $key = "__".$table."__".$start;
  my $db   = $self->db;

  for (my $status = $db->seq($key,$value,R_CURSOR);
       $status == 0;
       $status = $db->seq($key,$value,R_NEXT)) {

    my $feat = $self->{data}->get($value);
    last unless $do_while->($feat,$key);

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
  $options->{temp_file} = $self->_temp_file;

  my $results = $self->_get_features_by_search_options($search,$options);
  return Bio::DB::GFF::Adaptor::berkeleydb::iterator->new($results,$callback,$options->{temp_file});
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
  my ($write,$create) = @_;

  my $mode =  $create  ? "+>"
            : $write   ? "+>>"
            : "<";

  open my $F, $mode, $dbname or $class->throw("Could not open file '$dbname': $!");
  my $self = bless {
		    fh        => $F,
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
  print $fh pack("n/a*",$value) or $self->throw("An error occurred while updating the data file: $!");


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
  $offset ||= 0;
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
