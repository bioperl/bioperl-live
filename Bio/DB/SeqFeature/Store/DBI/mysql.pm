package Bio::DB::SeqFeature::Store::DBI::mysql;
# $Id$
use strict;

use base 'Bio::DB::SeqFeature::Store';
use Bio::DB::SeqFeature::Store::DBI::Iterator;
use DBI;
use Memoize;
use Cwd 'abs_path';
use Time::HiRes 'time';
use Bio::DB::GFF::Util::Rearrange 'rearrange';
use constant BINSIZE => 10_000;  # most features should be smaller than this
use constant DEBUG=>0;

# from the MySQL documentation...
# WARNING: if your sequence uses coordinates greater than 2 GB, you are out of luck!
use constant MAX_INT =>  2_147_483_647;
use constant MIN_INT => -2_147_483_648;

memoize('_typeid');
memoize('_locationid');
memoize('_attributeid');

###
# object initialization
#
sub init {
  my $self          = shift;
  my ($dsn,
      $is_temporary,
      $autoindex,
      $namespace,
      $binsize,
      $index_subfeatures,
      $dump_dir,
      $dbi_options) = rearrange(['DSN',
				 ['TEMP','TEMPORARY'],
				 'AUTOINDEX',
				 'NAMESPACE',
				 'BINSIZE',
				 'INDEX_SUBFEATURES',
				 ['DUMP_DIR','DUMPDIR'],
				 'OPTIONS',
				],@_);
  $dbi_options  ||= [];

  $dsn or $self->throw("Usage: ".__PACKAGE__."->init(-dsn => \$dbh || \$dsn)");

  my $dbh;
  if (ref $dsn) {
    $dbh = $dsn;
  } else {
    $dsn = "dbi:mysql:$dsn" unless $dsn =~ /^dbi:/;
    $dbh = DBI->connect($dsn,@$dbi_options);
  }
  $self->{dbh}       = $dbh;
  $self->{is_temp}   = $is_temporary;
  $self->{namespace} = $namespace;

  $self->maybe_create_meta();
  $self->default_settings;
  $self->binsize($binsize)                       if $binsize;
  $self->autoindex($autoindex)                   if defined $autoindex;
  $self->dumpdir($dump_dir)                      if $dump_dir;
  $self->init_tmp_database()                     if $self->is_temp;
}

sub _can_store_subFeatures { 1 }

sub table_definitions {
  my $self = shift;
  return {
	  feature => <<END,
(
  id         int(10) auto_increment primary key,
  indexed    tinyint default 1,
  object     MEDIUMBLOB not null
)
END
	  location => <<END,
(
  id       int(10)      not null,
  seqid    int(10)      not null,
  start    int          not null,
  end      int          not null,
  bin      int          not null,
  index(id),
  index(seqid,bin,start,end)
)
END

	  locationlist => <<END,
(
  id         int(10)       auto_increment primary key,
  seqname    varchar(50)   not null,
  index(seqname)
)
END

	  type => <<END,
(
  id       int(10)      primary key,
  typeid   int(10)      not null,
  index(typeid)
)
END

	  typelist => <<END,
(
  id       int(10) auto_increment primary key,
  tag      varchar(40)  not null,
  index(tag)
)
END
	  name => <<END,
(
  id           int(10)       not null,
  name         varchar(128)  not null,
  display_name tinyint       default 0,
  index(id),
  index(name)
)
END

	  attribute => <<END,
(
  id               int(10)       not null,
  attribute_id     int(10)   not null,
  attribute_value  text,
  index(id),
  index(attribute_id,attribute_value(10))
)
END

	  attributelist => <<END,
(
  id       int(10) auto_increment primary key,
  tag      varchar(50)  not null,
  index(tag)
)
END
	  parent2child => <<END,
(
  id               int(10)       not null,
  child            int(10)       not null,
  index(id,child)
)
END

	  meta => <<END,
(
  name      varchar(128) primary key,
  value     varchar(128) not null
)
END
	 };
}

###
# default settings -- will create and populate meta table if needed
#
sub default_settings {
  my $self = shift;
  $self->SUPER::default_settings;
  $self->binsize(BINSIZE);
  $self->autoindex(1);
  $self->dumpdir(abs_path('.'));
}


###
# retrieve database handle
#
sub dbh {
  my $self = shift;
  my $d    = $self->{dbh};
  $self->{dbh} = shift if @_;
  $d;
}

###
# get/set binsize
#
sub binsize {
  my $self = shift;
  my $d = $self->setting('binsize');
  $self->setting(binsize=>shift) if @_;
  $d;
}

###
# get/set directory for bulk load tables
#
sub dumpdir {
  my $self = shift;
  my $d = $self->{dumpdir};
  $self->{dumpdir} = abs_path(shift) if @_;
  $d;
}

###
# table namespace (multiple dbs in one mysql db)
#
sub namespace {
  my $self = shift;
  my $d    = $self->{namespace};
  $self->{namespace} = shift if @_;
  $d;
}

###
# find a path that corresponds to a dump table
#
sub dump_path {
  my $self = shift;
  my $table = $self->_qualify(shift);
  return "$self->{dumpdir}/$table.$$";
}

###
# make a filehandle (writeable) that corresponds to a dump table
#
sub dump_filehandle {
  my $self = shift;
  my $table = shift;
  eval "require IO::File" unless IO::File->can('new');
  my $path  = $self->dump_path($table);
  my $fh = $self->{filehandles}{$path} ||= IO::File->new(">$path");
  $fh;
}

###
# find the next ID for a feature (used only during bulk loading)
#
sub next_id {
  my $self = shift;
  $self->{max_id} ||= $self->max_id;
  return ++$self->{max_id};
}

###
# find the maximum ID for a feature (used only during bulk loading)
#
sub max_id {
  my $self = shift;
  my $sth  = $self->_prepare("SELECT max(id) from feature");
  $sth->execute or $self->throw($sth->errstr);
  my ($id) = $sth->fetchrow_array;
  $id;
}

###
# wipe database clean and reinstall schema
#
sub _init_database {
  my $self = shift;
  my $erase = shift;

  my $dbh    = $self->dbh;
  my $tables = $self->table_definitions;
  foreach (keys %$tables) {
    next if $_ eq 'meta';      # don't get rid of meta data!
    my $table = $self->_qualify($_);
    $dbh->do("DROP table IF EXISTS $table") if $erase;
    my $query = "CREATE TABLE IF NOT EXISTS $table $tables->{$_}";
    $dbh->do($query) or $self->throw($dbh->errstr);
  }
  $self->subfeatures_are_indexed(1);
  1;
}

sub maybe_create_meta {
  my $self = shift;
  my $table = $self->_qualify('meta');
  my $tables = $self->table_definitions;
  $self->dbh->do("CREATE TABLE IF NOT EXISTS $table $tables->{meta}");
}

sub init_tmp_database {
  my $self = shift;
  my $dbh    = $self->dbh;
  my $tables = $self->table_definitions;
  for my $t (keys %$tables) {
    my $table = $self->_qualify($t);
    my $query = "CREATE TEMPORARY TABLE $table $tables->{$t}";
    $dbh->do($query) or $self->throw($dbh->errstr);
  }
  1;
}

###
# use temporary tables
#
sub is_temp {
  shift->{is_temp};
}

sub _store {
  my $self    = shift;

  # special case for bulk updates
  return $self->_dump_store(@_) if $self->{bulk_update_in_progress};

  my $indexed = shift;
  my $count = 0;

  my $autoindex = $self->autoindex;

  my $dbh = $self->dbh;
  local $dbh->{RaiseError} = 1;
  $dbh->begin_work;
  eval {
    for my $obj (@_) {
      $self->replace($obj,$indexed);
      $self->update_indexes($obj) if $indexed && $autoindex;
      $count++;
    }
  };

  if ($@) {
    warn "Transaction aborted because $@";
    $dbh->rollback;
  }
  else {
    $dbh->commit;
  }

  # remember whether we are have ever stored a non-indexed feature
  unless ($indexed or $self->{indexed_flag}++) {
    $self->subfeatures_are_indexed(0);
  }
  $count;
}

# we memoize this in order to avoid making zillions of calls
sub autoindex {
  my $self = shift;

  # special case for bulk update -- need to build the indexes
  # at the same time we build the main feature table
  return 1 if $self->{bulk_update_in_progress};
  my $d = $self->setting('autoindex');
  $self->setting(autoindex=>shift) if @_;
  $d;
}

sub _start_bulk_update {
  my $self = shift;
  my $dbh  = $self->dbh;
  $self->{bulk_update_in_progress}++;
}

sub _finish_bulk_update {
  my $self = shift;
  my $dbh  = $self->dbh;
  my $dir = $self->{dumpdir} || '.';
  for my $table ('feature',$self->index_tables) {
    my $fh = $self->dump_filehandle($table);
    my $path = $self->dump_path($table);
    $fh->close;
    my $qualified_table = $self->_qualify($table);
    $dbh->do("LOAD DATA INFILE '$path' REPLACE INTO TABLE $qualified_table FIELDS OPTIONALLY ENCLOSED BY '\\''") 
      or $self->throw($dbh->errstr);
    unlink $path;
  }
  delete $self->{bulk_update_in_progress};
}


###
# Add a subparts to a feature. Both feature and all subparts must already be in database.
#
sub add_SeqFeature {
  my $self     = shift;

  # special purpose method for case when we are doing a bulk update
  return $self->_dump_add_SeqFeature(@_) if $self->{bulk_update_in_progress};

  my $parent   = shift;
  my @children = @_;

  my $dbh = $self->dbh;
  local $dbh->{RaiseError} = 1;

  my $child_table = $self->_parent2child_table();
  my $count = 0;

  my $sth = $self->_prepare(<<END);
REPLACE INTO $child_table (id,child) VALUES (?,?)
END

  my $parent_id = (ref $parent ? $parent->primary_id : $parent) 
    or $self->throw("$parent should have a primary_id");

  $dbh->begin_work or $self->throw($dbh->errstr);
  eval {
    for my $child (@children) {
      my $child_id = ref $child ? $child->primary_id : $child;
      defined $child_id or die "no primary ID known for $child";
      $sth->execute($parent_id,$child_id);
      $count++;
    }
  };

  if ($@) {
    warn "Transaction aborted because $@";
    $dbh->rollback;
  }
  else {
    $dbh->commit;
  }
  $sth->finish;
  $count;
}

sub get_SeqFeatures {
  my $self   = shift;
  my $parent = shift;
  my @types  = @_;

  my $parent_id = $parent->primary_id or $self->throw("$parent should have a primary_id");
  my $feature_table = $self->_feature_table;
  my $child_table   = $self->_parent2child_table();

  my @from  = ("$feature_table as f","$child_table as c");
  my @where = ('f.id=c.child','c.id=?');
  my @args  = $parent_id;

  if (@types) {
    my ($from,$where,undef,@a) = $self->_types_sql(\@types,'c.child');
    push @from,$from   if $from;
    push @where,$where if $where;
    push @args,@a;
  }

  my $from  = join ', ',@from;
  my $where = join ' AND ',@where;


  my $sth = $self->_prepare(<<END) or $self->throw($self->dbh->errstr);
SELECT f.id,f.object
  FROM $from
  WHERE $where
END

  $sth->execute(@args) or $self->throw($sth->errstr);
  return $self->_sth2objs($sth);
}

###
# add namespace to tablename
#
sub _qualify {
  my $self = shift;
  my $table_name = shift;
  my $namespace = $self->namespace;
  return $table_name unless defined $namespace;
  return "${namespace}_${table_name}";
}

###
# Fetch a Bio::SeqFeatureI from database using its primary_id
#
sub _fetch {
  my $self       = shift;
  @_ or $self->throw("usage: fetch(\$primary_id)");
  my $primary_id = shift;
  my $features = $self->_feature_table;
  my $sth = $self->_prepare(<<END);
SELECT id,object FROM $features WHERE id=?
END
  $sth->execute($primary_id) or $self->throw($sth->errstr);
  my $obj = $self->_sth2obj($sth);
  $sth->finish;
  $obj;
}

###
# Efficiently fetch a series of IDs from the database
# Can pass an array or an array ref
#
sub _fetch_many {
  my $self       = shift;
  @_ or $self->throw('usage: fetch_many($id1,$id2,$id3...)');
  my $ids = join ',',map {ref($_) ? @$_ : $_} @_ or return;
  my $features = $self->_feature_table;

  my $sth = $self->_prepare(<<END);
SELECT id,object FROM $features WHERE id IN ($ids)
END
  $sth->execute() or $self->throw($sth->errstr);
  return $self->_sth2objs($sth);
}

sub _features {
  my $self = shift;
  my ($seq_id,$start,$end,
      $name,$class,$allow_aliases,
      $types,
      $attributes,
      $range_type,
      $fromtable,
      $iterator
     ) = rearrange([['SEQID','SEQ_ID','REF'],'START',['STOP','END'],
		    'NAME','CLASS','ALIASES',
		    ['TYPES','TYPE','PRIMARY_TAG'],
		    ['ATTRIBUTES','ATTRIBUTE'],
		    'RANGE_TYPE',
		    'FROM_TABLE',
		    'ITERATOR',
		   ],@_);

  my (@from,@where,@args,@group);
  $range_type ||= 'overlaps';

  my $feature_table         = $self->_feature_table;
  @from = "$feature_table as f";

  if (defined $name) {
    # hacky backward compatibility workaround
    $name = "$class:$name" if defined $class;
    my ($from,$where,$group,@a) = $self->_name_sql($name,$allow_aliases,'f.id');
    push @from,$from   if $from;
    push @where,$where if $where;
    push @group,$group if $group;
    push @args,@a;
  }

  if (defined $seq_id) {
    my ($from,$where,$group,@a) = $self->_location_sql($seq_id,$start,$end,$range_type,'f.id');
    push @from,$from   if $from;
    push @where,$where if $where;
    push @group,$group if $group;
    push @args,@a;
  }

  if (defined($types)) {
    my ($from,$where,$group,@a) = $self->_types_sql($types,'f.id');
    push @from,$from   if $from;
    push @where,$where if $where;
    push @group,$group if $group;
    push @args,@a;
  }

  if (defined $attributes) {
    my ($from,$where,$group,@a) = $self->_attributes_sql($attributes,'f.id');
    push @from,$from    if $from;
    push @where,$where  if $where;
    push @group,$group  if $group;
    push @args,@a;
  }

  if (defined $fromtable) {
    my ($from,$where,$group,@a) = $self->_from_table_sql($fromtable,'f.id');
    push @from,$from    if $from;
    push @where,$where  if $where;
    push @group,$group  if $group;
    push @args,@a;
  }

  # if no other criteria are specified, then
  # only fetch indexed (i.e. top level objects)
  @where = 'indexed=1' unless @where;

  my $from  = join ', ',@from;
  my $where = join ' AND ',map {"($_)"} @where;
  my $group = join ', ',@group;
  $group    = "GROUP BY $group" if @group;

  my $query = <<END;
SELECT DISTINCT f.id,f.object
  FROM $from
  WHERE $where
  $group
END

  $self->_print_query($query,@args) if DEBUG;

  my $sth = $self->_prepare($query);
  $sth->execute(@args) or $self->throw($sth->errstr);
  return $iterator ? Bio::DB::SeqFeature::Store::DBI::Iterator->new($sth,$self) : $self->_sth2objs($sth);
}

sub _name_sql {
  my $self = shift;
  my ($name,$allow_aliases,$join) = @_;
  my $name_table   = $self->_name_table;

  my $from  = "$name_table as n";
  my ($match,$string) = $self->_match_sql($name);

  my $where = "n.id=$join AND n.name $match";
  $where   .= " AND n.display_name>0" unless $allow_aliases;
  return ($from,$where,'',$string);
}

sub _match_sql {
  my $self = shift;
  my $name = shift;

  my ($match,$string);
  if ($name =~ /(?:^|[^\\])[*?]/) {
    $name =~ s/(^|[^\\])([%_])/$1\\$2/g;
    $name =~ s/(^|[^\\])\*/$1%/g;
    $name =~ s/(^|[^\\])\?/$1_/g;
    $match = "LIKE ?";
    $string  = $name;
  } else {
    $match = "= ?";
    $string  = $name;
  }
  return ($match,$string);
}

sub _from_table_sql {
  my $self = shift;
  my ($from_table,$join) = @_;
  my $from  = "$from_table as ft";
  my $where = "ft.id=$join";
  return ($from,$where,'');
}

sub _attributes_sql {
  my $self = shift;
  my ($attributes,$join) = @_;

  my ($wf,@bind_args)       = $self->_make_attribute_where('a','al',$attributes);
  my ($group_by,@group_args)= $self->_make_attribute_group('a',$attributes);

  my $attribute_table       = $self->_attribute_table;
  my $attributelist_table   = $self->_attributelist_table;

  my $from = "$attribute_table as a, $attributelist_table as al";

  my $where = <<END;
  a.id=$join
  AND   a.attribute_id=al.id
  AND ($wf)
END

  my $group = $group_by;

  my @args  = (@bind_args,@group_args);
  return ($from,$where,$group,@args);
}

sub _types_sql {
  my $self  = shift;
  my ($types,$join) = @_;
  my ($primary_tag,$source_tag);

  my @types = ref $types eq 'ARRAY' ?  @$types : $types;

  my $type_table    = $self->_type_table;
  my $typelist      = $self->_typelist_table;
  my $from = "$type_table AS t,$typelist AS tl";

  my (@matches,@args);

  for my $type (@types) {

    if (ref $type && $type->isa('Bio::DB::GFF::Typename')) {
      $primary_tag = $type->method;
      $source_tag  = $type->source;
    } else {
      ($primary_tag,$source_tag) = split ':',$type;
    }

    if (defined $source_tag) {
      push @matches,"tl.tag=?";
      push @args,"$primary_tag:$source_tag";
    } else {
      push @matches,"tl.tag LIKE ?";
      push @args,"$primary_tag:%";
    }
  }
  my $matches = join ' OR ',@matches;

  # the join table is named "f"
  my $where = <<END;
   t.id=$join
   AND   tl.id=t.typeid
   AND   ($matches)
END

  return ($from,$where,'',@args);
}

sub _location_sql {
  my $self = shift;
  my ($seq_id,$start,$end,$range_type,$join) = @_;

  my $location_table = $self->_location_table;
  my $location_list  = $self->_locationlist_table;

  # the additional join on the location_list table badly impacts performance
  # so we build a copy of the table in memory
  my $seqid = $self->_locationid($seq_id) || 0; # zero is an invalid primary ID, so will return empty

# using a join (slow)
#  my $from = "$location_table as l,$location_list as ll";
  my $from  = "$location_table as l";

  $start = MIN_INT unless defined $start;
  $end   = MAX_INT unless defined $end;

  my $bin_min = int $start/BINSIZE;
  my $bin_max = int $end/BINSIZE;

  my ($range,@range_args);
  if ($range_type eq 'overlaps') {
    $range = 'bin BETWEEN ? AND ? AND start<=? AND end>=?';
    @range_args = ($bin_min,$bin_max,$end,$start);
  } elsif ($range_type eq 'contains') {
    $range = 'bin BETWEEN ? AND ? AND start>=? AND end<=?';
    @range_args = ($bin_min,$bin_max,$start,$end);
  } elsif ($range_type eq 'contained_in') {
    $range = 'bin<=? OR bin>=? AND start<=? AND end>=?';
    @range_args = ($bin_min,$bin_max,$start,$end);
  } else {
    $self->throw("range_type must be one of 'overlaps', 'contains' or 'contained_in'");
  }


# using a join (slow)
#  my $where = <<END;
#   l.id=$join
#   AND   ll.seqname=?
#   AND   ll.id=l.seqid
#   AND   $range
#END

  my $where = <<END;
   l.id=$join
   AND   l.seqid=?
   AND   $range
END

  my $group = '';

# using a join (slow)
#  my @args  = ($seq_id,@range_args);

  my @args  = ($seqid,@range_args);
  return ($from,$where,$group,@args);
}

###
# force reindexing
#
sub reindex {
  my $self = shift;
  my $from_update_table = shift;  # if present, will take ids from "update_table"

  my $dbh  = $self->dbh;
  my $count = 0;
  my $now;
  my $last_time = time();

  # tell _delete_index() not to bother removing the index rows corresponding
  # to each individual feature
  local $self->{reindexing} = 1;

  $dbh->begin_work;
  eval {
    my $update = $from_update_table;
    for my $table ($self->index_tables) {
      my $query = $from_update_table ? "DELETE $table FROM $table,$update WHERE $table.id=$update.id"
	                             : "DELETE FROM $table";
      $dbh->do($query);
      $dbh->do("ALTER TABLE $table DISABLE KEYS");
    }
    my $iterator = $self->get_seq_stream(-from_table=>$from_update_table ? $update : undef);
    while (my $f = $iterator->next_seq) {
      if (++$count %1000 == 0) {
	$now = time();
	my $elapsed = sprintf(" in %5.2fs",$now - $last_time);
	$last_time = $now;
	print STDERR "$count features indexed$elapsed...",' 'x60;
	print STDERR -t STDOUT && !$ENV{EMACS} ? "\r" : "\n";
      }
      $self->update_indexes($f);
    }
  };
  for my $table ($self->index_tables) {
    $dbh->do("ALTER TABLE $table ENABLE KEYS");
  }
  if (@_) {
    warn "Couldn't complete transaction: $@";
    $dbh->rollback;
    return;
  } else {
    $dbh->commit;
    return 1;
  }
}

sub optimize {
  my $self = shift;
  $self->dbh->do("ANALYZE TABLE $_") foreach $self->index_tables;
}

sub all_tables {
  my $self = shift;
  my @index_tables = $self->index_tables;
  my $feature_table = $self->_feature_table;
  return ($feature_table,@index_tables);
}

sub index_tables {
  my $self = shift;
  return map {$self->_qualify($_)} qw(location type name attribute parent2child)
}

sub _firstid {
  my $self = shift;
  my $features = $self->_feature_table;
  my $query = <<END;
SELECT min(id) FROM $features
END
  my $sth=$self->_prepare($query);
  $sth->execute();
  my ($first) = $sth->fetchrow_array;
  $sth->finish;
  $first;
}

sub _nextid {
  my $self = shift;
  my $lastkey = shift;
  my $features = $self->_feature_table;
  my $query = <<END;
SELECT min(id) FROM $features WHERE id>?
END
  my $sth=$self->_prepare($query);
  $sth->execute($lastkey);
  my ($next) = $sth->fetchrow_array;
  $sth->finish;
  $next;
}

sub _existsid {
  my $self = shift;
  my $key  = shift;
  my $features = $self->_feature_table;
  my $query = <<END;
SELECT count(*) FROM $features WHERE id=?
END
  my $sth=$self->_prepare($query);
  $sth->execute($key);
  my ($count) = $sth->fetchrow_array;
  $sth->finish;
  $count > 0;
}

sub _deleteid {
  my $self = shift;
  my $key  = shift;
  my $dbh = $self->dbh;
  for my $table ($self->all_tables) {
    $dbh->do("DELETE FROM $table WHERE id=$key");
  }
}

sub _clearall {
  my $self = shift;
  my $dbh = $self->dbh;
  for my $table ($self->all_tables) {
    $dbh->do("DELETE FROM $table");
  }
}

sub _featurecount {
  my $self = shift;
  my $dbh = $self->dbh;
  my $features = $self->_feature_table;
  my $query = <<END;
SELECT count(*) FROM $features
END
  my $sth=$self->_prepare($query);
  $sth->execute();
  my ($count) = $sth->fetchrow_array;
  $sth->finish;
  $count;
}

sub setting {
  my $self = shift;
  my ($variable_name,$value) = @_;
  my $meta  = $self->_meta_table;

  if (defined $value) {
    my $query = <<END;
REPLACE INTO $meta (name,value) VALUES (?,?)
END
    my $sth = $self->_prepare($query);
    $sth->execute($variable_name,$value) or $self->throw($sth->errstr);
    $sth->finish;
    $self->{settings_cache}{$variable_name} = $value;
  }
  else {
    return $self->{settings_cache}{$variable_name} if exists $self->{settings_cache}{$variable_name};
    my $query = <<END;
SELECT value FROM $meta as m WHERE m.name=?
END
    my $sth = $self->_prepare($query);
    $sth->execute($variable_name) or $self->throw($sth->errstr);
    my ($value) = $sth->fetchrow_array;
    $sth->finish;
    return $self->{settings_cache}{$variable_name} = $value;
  }
}

###
# Replace Bio::SeqFeatureI into database.
#
sub replace {
  my $self       = shift;
  my $object     = shift;
  my $index_flag = shift || undef;

  # ?? shouldn't need to do this
  # $self->_load_class($object);
  my $id = $object->primary_id;
  my $features = $self->_feature_table;

  my $sth = $self->_prepare(<<END);
REPLACE INTO $features (id,object,indexed) VALUES (?,?,?)
END
  $sth->execute($id,$self->_freeze($object),$index_flag||0) or $self->throw($sth->errstr);
  my $dbh = $self->dbh;
  $object->primary_id($dbh->{mysql_insertid});
  $self->flag_for_indexing($dbh->{mysql_insertid}) if $self->{bulk_update_in_progress};
}

###
# Insert one Bio::SeqFeatureI into database. primary_id must be undef
#
sub insert {
  my $self = shift;
  my $object = shift;
  my $index_flag = shift || 0;

  $self->_load_class($object);
  defined $object->primary_id and $self->throw("$object already has a primary id");

  my $features = $self->_feature_table;
  my $sth = $self->_prepare(<<END);
INSERT INTO $features (id,object,indexed) VALUES (?,?,?)
END
  $sth->execute(undef,$self->_freeze($object),$index_flag) or $self->throw($sth->errstr);
  my $dbh = $self->dbh;
  $object->primary_id($dbh->{mysql_insertid});
  $self->flag_for_indexing($dbh->{mysql_insertid}) if $self->{bulk_update_in_progress};
}

###
# Update a Bio::SeqFeatureI into database. primary_id must NOT be undef
#
sub _update {
  my $self   = shift;
  my $object        = shift;
  my $primary_id    = shift;

  my $features = $self->_feature_table;
  my $sth = $self->_prepare(<<END);
UPDATE $features SET object=? where id=?
END
  # If I had Mysql 4.1 or higher, I could do this:
  # INSERT INTO features (id,object) VALUES (?,?)
  # ON DUPLICATE KEY UPDATE object=?
  $sth->execute($self->_freeze($object),$primary_id)       or $self->throw($self->dbh->errstr);
  $self->flag_for_indexing($self->dbh->{mysql_insertid})   if $self->{bulk_update_in_progress};
}

###
# This subroutine flags the given primary ID for later reindexing
#
sub flag_for_indexing {
  my $self = shift;
  my $id   = shift;
  my $needs_updating = $self->_update_table;
  my $sth = $self->_prepare("REPLACE INTO $needs_updating VALUES (?)");
  $sth->execute($id) or $self->throw($self->dbh->errstr);
}

###
# Update indexes for given object
#
sub update_indexes {
  my $self = shift;
  my $obj  = shift;
  defined (my $id   = $obj->primary_id) or return;

  if ($self->{bulk_update_in_progress}) {
    $self->_dump_update_name_index($obj,$id);
    $self->_dump_update_type_index($obj,$id);
    $self->_dump_update_location_index($obj,$id);
    $self->_dump_update_attribute_index($obj,$id);
  } else {
    $self->_update_name_index($obj,$id);
    $self->_update_type_index($obj,$id);
    $self->_update_location_index($obj,$id);
    $self->_update_attribute_index($obj,$id);
  }
}

sub _update_name_index {
  my $self = shift;
  my ($obj,$id) = @_;
  my $name = $self->_name_table;
  my $primary_id = $obj->primary_id;

  $self->_delete_index($name,$id);
  my ($names,$aliases) = $self->names($obj);

  my $sth = $self->_prepare("INSERT INTO $name (id,name,display_name) VALUES (?,?,?)");

  $sth->execute($id,$_,1) or $self->throw($sth->errstr)   foreach @$names;
  $sth->execute($id,$_,0) or $self->throw($sth->errstr) foreach @$aliases;
  $sth->finish;
}

sub _update_type_index {
  my $self = shift;
  my ($obj,$id) = @_;
  my $type    = $self->_type_table;
  $self->_delete_index($type,$id);

  my $primary_tag = $obj->primary_tag;
  my $source_tag  = $obj->source_tag || '';
  return unless defined $primary_tag;

  $primary_tag    .= ":$source_tag";
  my $sth = $self->_prepare("INSERT INTO $type (id,typeid) VALUES (?,?)");
  my $typeid = $self->_typeid($primary_tag,1);
  $sth->execute($id,$typeid) or $self->throw($sth->errstr);
  $sth->finish;
}

sub _update_location_index {
  my $self = shift;
  my ($obj,$id) = @_;
  my $location    = $self->_location_table;
  $self->_delete_index($location,$id);

  my $seqname     = $obj->seq_id || '';
  my $start       = $obj->start  || '';
  my $end         = $obj->end    || '';
  my $bin_min     = int $start/BINSIZE;
  my $bin_max     = int $end/BINSIZE;
  my $sth         = $self->_prepare("INSERT INTO $location (id,seqid,start,end,bin) VALUES (?,?,?,?,?)");
  my $seqid       = $self->_locationid($seqname,1);
  for (my $bin=$bin_min; $bin<=$bin_max; $bin++) {
    $sth->execute($id,$seqid,$start,$end,$bin) or $self->throw($sth->errstr);
  }
  $sth->finish;
}

sub _update_attribute_index {
  my $self = shift;
  my ($obj,$id) = @_;
  my $attribute = $self->_attribute_table;
  $self->_delete_index($attribute,$id);

  my $sth = $self->_prepare("INSERT INTO $attribute (id,attribute_id,attribute_value) VALUES (?,?,?)");
  for my $tag ($obj->all_tags) {
    my $tagid = $self->_attributeid($tag);
    for my $value ($obj->each_tag_value($tag)) {
      $sth->execute($id,$tagid,$value) or $self->throw($sth->errstr);
    }
  }
  $sth->finish;
}


sub _genericid {
  my $self = shift;
  my ($table,$namefield,$name,$add_if_missing) = @_;
  my $qualified_table = $self->_qualify($table);
  my $sth = $self->_prepare(<<END);
SELECT id FROM $qualified_table WHERE $namefield=?
END
  $sth->execute($name) or die $sth->errstr;
  my ($id) = $sth->fetchrow_array;
  $sth->finish;
  return $id if defined $id;
  return     unless $add_if_missing;

  $sth = $self->_prepare(<<END);
INSERT INTO $qualified_table ($namefield) VALUES (?)
END
  $sth->execute($name) or die $sth->errstr;
  my $dbh = $self->dbh;
  return $dbh->{mysql_insertid};
}

sub _typeid {
  shift->_genericid('typelist','tag',shift,1);
}
sub _locationid {
  shift->_genericid('locationlist','seqname',shift,1);
}
sub _attributeid {
  shift->_genericid('attributelist','tag',shift,1);
}

sub names {
  my $self = shift;
  my $obj  = shift;

  my $primary_id = $obj->primary_id;
  my @names = $obj->display_name;
  push @names,eval{$obj->get_tag_values('Name')};
  push @names,eval{$obj->get_tag_values('ID')};
  @names = grep {defined $_ && $_ ne $primary_id} @names;

  my @aliases = grep {defined} eval{$obj->get_tag_values('Alias')};

  return (\@names,\@aliases);
}

sub _delete_index {
  my $self = shift;
  my ($table_name,$id) = @_;
  return if $self->{reindexing};
  my $sth = $self->_prepare("DELETE FROM $table_name WHERE id=?") or $self->throw($self->dbh->errstr);
  $sth->execute($id);
}

# given a statement handler that is expected to return rows of (id,object)
# unthaw each object and return a list of 'em
sub _sth2objs {
  my $self = shift;
  my $sth  = shift;
  my @result;
  while (my ($id,$o) = $sth->fetchrow_array) {
    my $obj = $self->_thaw($o,$id);
    push @result,$obj;
  }
  $sth->finish;
  return @result;
}

# given a statement handler that is expected to return rows of (id,object)
# unthaw each object and return a list of 'em
sub _sth2obj {
  my $self = shift;
  my $sth  = shift;
  my ($id,$o) = $sth->fetchrow_array;
  return unless $o;
  my $obj = $self->_thaw($o,$id);
  $obj;
}

sub _prepare {
  my $self = shift;
  my $query = shift;
  my $dbh   = $self->dbh;
  my $sth   = $dbh->prepare_cached($query) or $self->throw($dbh->errstr);
  $sth;
}


####################################################################################################
# SQL Fragment generators
####################################################################################################

sub _feature_table       {  shift->_qualify('feature')  }
sub _location_table      {  shift->_qualify('location') }
sub _locationlist_table  {  shift->_qualify('locationlist') }
sub _type_table          {  shift->_qualify('type')     }
sub _typelist_table      {  shift->_qualify('typelist') }
sub _name_table          {  shift->_qualify('name')     }
sub _attribute_table     {  shift->_qualify('attribute')}
sub _attributelist_table {  shift->_qualify('attributelist')}
sub _parent2child_table  {  shift->_qualify('parent2child')}
sub _meta_table          {  shift->_qualify('meta')}
sub _update_table        {  shift->_qualify('update_table')}

sub _make_attribute_where {
  my $self                     = shift;
  my ($attributetable,$attributenametable,$attributes) = @_;
  my @args;
  my @sql;
  my $dbh = $self->dbh;
  foreach (keys %$attributes) {
    my @match_values;
    my @values = ref($attributes->{$_}) && ref($attributes->{$_}) eq 'ARRAY' ? @{$attributes->{$_}} : $attributes->{$_};
    foreach (@values) {  # convert * into % for wildcard matches
      s/\*/%/g;
    }
    my $match  = join ' OR ',map {
      /%/ ? "$attributetable.attribute_value LIKE ?"
	  : "$attributetable.attribute_value=?"
    } @values;
    push @sql,"($attributenametable.tag=? AND ($match))";
    push @args,($_,@values);
  }
  return (join(' OR ',@sql),@args);
}

sub _make_attribute_group {
  my $self                     = shift;
  my ($table_name,$attributes) = @_;
  my $key_count = keys %$attributes or return;
  return "$table_name.id HAVING count($table_name.id)>?",$key_count-1;
}

sub _print_query {
  my $self = shift;
  my ($query,@args) = @_;
  while ($query =~ /\?/) {
    my $arg = $self->dbh->quote(shift @args);
    $query =~ s/\?/$arg/;
  }
  warn $query,"\n";
}

###
# special-purpose store for bulk loading - write to a file rather than to the db
#
sub _dump_store {
  my $self    = shift;
  my $indexed = shift;

  my $count = 0;
  my $store_fh = $self->dump_filehandle('feature');
  my $dbh      = $self->dbh;

  my $autoindex = $self->autoindex;

  for my $obj (@_) {
    my $id       = $self->next_id;
    print $store_fh join("\t",$id,$indexed,$dbh->quote($self->_freeze($obj))),"\n";
    $obj->primary_id($id);
    $self->update_indexes($obj) if $indexed && $autoindex;
    $count++;
  }

  # remember whether we are have ever stored a non-indexed feature
  unless ($indexed or $self->{indexed_flag}++) {
    $self->subfeatures_are_indexed(0);
  }
  $count;
}

sub _dump_add_SeqFeature {
  my $self     = shift;
  my $parent   = shift;
  my @children = @_;

  my $dbh = $self->dbh;
  my $fh = $self->dump_filehandle('parent2child');
  my $parent_id = (ref $parent ? $parent->primary_id : $parent) 
    or $self->throw("$parent should have a primary_id");
  my $count = 0;

  for my $child_id (@children) {
    print $fh join("\t",$parent_id,$child_id),"\n";
    $count++;
  }
  $count;
}

sub _dump_update_name_index {
  my $self = shift;
  my ($obj,$id) = @_;
  my $fh      = $self->dump_filehandle('name');
  my $dbh     = $self->dbh;
  my ($names,$aliases) = $self->names($obj);
  print $fh join("\t",$id,$dbh->quote($_),1),"\n" foreach @$names;
  print $fh join("\t",$id,$dbh->quote($_),0),"\n" foreach @$aliases;
}

sub _dump_update_type_index {
  my $self = shift;
  my ($obj,$id) = @_;

  my $primary_tag = $obj->primary_tag;
  my $source_tag  = $obj->source_tag || '';
  return unless defined $primary_tag;
  my $dbh         = $self->dbh;

  $primary_tag    .= ":$source_tag";
  my $fh  = $self->dump_filehandle('type');
  my $typeid = $self->_typeid($primary_tag,1);
  print $fh join ("\t",$id,$typeid),"\n";
}

sub _dump_update_location_index {
  my $self = shift;
  my ($obj,$id) = @_;

  my $seqname     = $obj->seq_id || '';
  my $start       = $obj->start  || '';
  my $end         = $obj->end    || '';
  my $bin_min     = int $start/BINSIZE;
  my $bin_max     = int $end/BINSIZE;
  my $fh          = $self->dump_filehandle('location');
  my $seqid       = $self->_locationid($seqname,1);
  for (my $bin=$bin_min; $bin<=$bin_max; $bin++) {
    print $fh join("\t",$id,$seqid,$start,$end,$bin),"\n";
  }
}

sub _dump_update_attribute_index {
  my $self = shift;
  my ($obj,$id) = @_;
  my $fh        = $self->dump_filehandle('attribute');
  my $dbh       = $self->dbh;
  for my $tag ($obj->all_tags) {
    my $tagid = $self->_attributeid($tag);
    for my $value ($obj->each_tag_value($tag)) {
      print $fh join("\t",$id,$tagid,$dbh->quote($value)),"\n";
    }
  }
}

1;
