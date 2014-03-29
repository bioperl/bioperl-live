package Bio::DB::SeqFeature::Store::DBI::SQLite;

#$Id$

=head1 NAME

Bio::DB::SeqFeature::Store::DBI::SQLite -- SQLite implementation of Bio::DB::SeqFeature::Store

=head1 SYNOPSIS

  use Bio::DB::SeqFeature::Store;

  # Open the sequence database
  my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::SQLite',
                                           -dsn     => '/path/to/database.db');

  # get a feature from somewhere
  my $feature = Bio::SeqFeature::Generic->new(...);

  # store it
  $db->store($feature) or die "Couldn't store!";

  # primary ID of the feature is changed to indicate its primary ID
  # in the database...
  my $id = $feature->primary_id;

  # get the feature back out
  my $f  = $db->fetch($id);

  # change the feature and update it
  $f->start(100);
  $db->update($f) or die "Couldn't update!";

  # searching...
  # ...by id
  my @features = $db->fetch_many(@list_of_ids);

  # ...by name
  @features = $db->get_features_by_name('ZK909');

  # ...by alias
  @features = $db->get_features_by_alias('sma-3');

  # ...by type
  @features = $db->get_features_by_name('gene');

  # ...by location
  @features = $db->get_features_by_location(-seq_id=>'Chr1',-start=>4000,-end=>600000);

  # ...by attribute
  @features = $db->get_features_by_attribute({description => 'protein kinase'})

  # ...by the GFF "Note" field
  @result_list = $db->search_notes('kinase');

  # ...by arbitrary combinations of selectors
  @features = $db->features(-name => $name,
                            -type => $types,
                            -seq_id => $seqid,
                            -start  => $start,
                            -end    => $end,
                            -attributes => $attributes);

  # ...using an iterator
  my $iterator = $db->get_seq_stream(-name => $name,
                                     -type => $types,
                                     -seq_id => $seqid,
                                     -start  => $start,
                                     -end    => $end,
                                     -attributes => $attributes);

  while (my $feature = $iterator->next_seq) {
    # do something with the feature
  }

  # ...limiting the search to a particular region
  my $segment  = $db->segment('Chr1',5000=>6000);
  my @features = $segment->features(-type=>['mRNA','match']);

  # getting & storing sequence information
  # Warning: this returns a string, and not a PrimarySeq object
  $db->insert_sequence('Chr1','GATCCCCCGGGATTCCAAAA...');
  my $sequence = $db->fetch_sequence('Chr1',5000=>6000);

  # what feature types are defined in the database?
  my @types    = $db->types;

  # create a new feature in the database
  my $feature = $db->new_feature(-primary_tag => 'mRNA',
                                 -seq_id      => 'chr3',
                                 -start      => 10000,
                                 -end        => 11000);

=head1 DESCRIPTION

Bio::DB::SeqFeature::Store::SQLite is the SQLite adaptor for
Bio::DB::SeqFeature::Store. You will not create it directly, but
instead use Bio::DB::SeqFeature::Store-E<gt>new() to do so.

See L<Bio::DB::SeqFeature::Store> for complete usage instructions.

=head2 Using the SQLite adaptor

To establish a connection to the database, call
Bio::DB::SeqFeature::Store-E<gt>new(-adaptor=E<gt>'DBI::SQLite',@more_args). The
additional arguments are as follows:

  Argument name       Description
  -------------       -----------

 -dsn              The path to the SQLite database file.

 -namespace        A prefix to attach to each table. This allows you
                   to have several virtual databases in the same
                   physical database.

 -temp             Boolean flag. If true, a temporary database
                   will be created and destroyed as soon as
                   the Store object goes out of scope. (synonym -temporary)

 -autoindex        Boolean flag. If true, features in the database will be
                   reindexed every time they change. This is the default.


 -tmpdir           Directory in which to place temporary files during "fast" loading.
                   Defaults to File::Spec->tmpdir(). (synonyms -dump_dir, -dumpdir, -tmp)

 -dbi_options      A hashref to pass to DBI->connect's 4th argument, the "attributes."
                   (synonyms -options, -dbi_attr)

 -write            Pass true to open database for writing or updating.

If successful, a new instance of
Bio::DB::SeqFeature::Store::DBI::SQLite will be returned.

In addition to the standard methods supported by all well-behaved
Bio::DB::SeqFeature::Store databases, several following
adaptor-specific methods are provided. These are described in the next
sections.

=cut

use strict;

use base 'Bio::DB::SeqFeature::Store::DBI::mysql';
use Bio::DB::SeqFeature::Store::DBI::Iterator;
use DBI qw(:sql_types);
use Memoize;
use Cwd qw(abs_path getcwd);
use Bio::DB::GFF::Util::Rearrange 'rearrange';
use Bio::SeqFeature::Lite;
use File::Spec;
use constant DEBUG=>0;
use constant EXPERIMENTAL_COVERAGE=>1;

# Using same limits as MySQL adaptor so I don't have to make something up.
use constant MAX_INT =>  2_147_483_647;
use constant MIN_INT => -2_147_483_648;
use constant SUMMARY_BIN_SIZE =>  1000;  # we checkpoint coverage this often, about 20 meg overhead per feature type on hg
use constant USE_SPATIAL=>0;

# The binning scheme places each feature into a bin.
# Bins are variably sized as powers of two. For example,
# there are 585 bins of size 2**17 (131072 bases)
my (@BINS,%BINS);
{
    @BINS = map {2**$_} (17, 20, 23, 26, 29);  # TO DO: experiment with different bin sizes
    my $start=0;
    for my $b (sort {$b<=>$a} @BINS) {
	$BINS{$b} = $start;
	$start   += $BINS[-1]/$b;
    }
}


# my %BINS = (
#     2**11 => 37449,
#     2**14 =>  4681,
#     2**17 =>   585,
#     2**20 =>    73,
#     2**23 =>     9,
#     2**26 =>     1,
#     2**29 =>     0
# );
# my @BINS = sort {$a<=>$b} keys %BINS;

sub calculate_bin {
    my $self = shift;
    my ($start,$end) = @_;

    my $len = $end - $start;
    for my $bin (@BINS) {
	next if $len > $bin;
	# possibly fits here
	my $binstart = int $start/$bin;
	my $binend   = int $end/$bin;
	return $binstart+$BINS{$bin} if $binstart == $binend;
    }

    die "unreasonable coordinates ",$start+1,"..$end";
}

sub search_bins {
    my $self = shift;
    my ($start,$end) = @_;
    my @results;

    for my $bin (@BINS) {
	my $binstart = int $start/$bin;
	my $binend   = int $end/$bin;
	push @results,$binstart+$BINS{$bin}..$binend+$BINS{$bin};
    }
    return @results;
}


###
# object initialization
#
sub init {
  my $self          = shift;
  my ($dsn,
      $is_temporary,
      $autoindex,
      $namespace,
      $dump_dir,
      $user,
      $pass,
      $dbi_options,
      $writeable,
      $create,
     ) = rearrange(['DSN',
		    ['TEMP','TEMPORARY'],
		    'AUTOINDEX',
		    'NAMESPACE',
		    ['DUMP_DIR','DUMPDIR','TMP','TMPDIR'],
		    'USER',
		    ['PASS','PASSWD','PASSWORD'],
		    ['OPTIONS','DBI_OPTIONS','DBI_ATTR'],
		    ['WRITE','WRITEABLE'],
		    'CREATE',
		   ],@_);
  $dbi_options  ||= {};
  $writeable    = 1 if $is_temporary or $dump_dir;

  $dsn or $self->throw("Usage: ".__PACKAGE__."->init(-dsn => \$dbh || \$dsn)");

  my $dbh;
  if (ref $dsn) {
    $dbh = $dsn;
  } else {
    $dsn = "dbi:SQLite:$dsn" unless $dsn =~ /^dbi:/;
    $dbh = DBI->connect($dsn,$user,$pass,$dbi_options) or $self->throw($DBI::errstr);
    $dbh->do("PRAGMA synchronous = OFF;"); # makes writes much faster
    $dbh->do("PRAGMA temp_store = MEMORY;"); # less disk I/O; some speedup
    $dbh->do("PRAGMA cache_size = 20000;"); # less disk I/O; some speedup
    # Keep track of database file location
    my $cwd = getcwd;
    my ($db_file) = ($dsn =~ m/(?:db(?:name)?|database)=(.+)$/);
    $self->{dbh_file} = "$cwd/$db_file";
  }
  $self->{dbh}       = $dbh;
  $self->{is_temp}   = $is_temporary;
  $self->{namespace} = $namespace;
  $self->{writeable} = $writeable;

  $self->default_settings;
  $self->autoindex($autoindex)                   if defined $autoindex;
  $self->dumpdir($dump_dir)                      if $dump_dir;
  if ($self->is_temp) {
    $self->init_tmp_database();
  } elsif ($create) {
    $self->init_database('erase');
  }
}

sub table_definitions {
  my $self = shift;
  my $defs =
      {
	  feature => <<END,
(
  id        integer primary key autoincrement,
  typeid    integer not null,
  strand    integer default 0,
  "indexed" integer default 1,
  object    blob not null
)
END

	  locationlist => <<END,
(
  id      integer primary key autoincrement,
  seqname text    not null
);
create index index_locationlist on locationlist (seqname);
END

	  typelist => <<END,
(
  id  integer primary key autoincrement,
  tag text    not null collate nocase
);
create index index_typelist on typelist (tag);
END
	  name => <<END,
(
  id           integer not null,
  name         text    not null collate nocase,
  display_name integer default 0
);
create index index_name_id on name(id);
create index index_name_name on name(name);
END

	  attribute => <<END,
(
  id              integer not null,
  attribute_id    integer not null,
  attribute_value text collate nocase
);
create index index_attribute_id    on attribute(attribute_id);
create index index_attribute_value on attribute(attribute_value);
END

	  attributelist => <<END,
(
  id  integer primary key autoincrement,
  tag text    not null
);
create index index_attributelist_id  on attributelist(id);
create index index_attributelist_tag on attributelist(tag);
END
	  parent2child => <<END,
(
  id    integer not null,
  child integer not null
);
create unique index index_parent2child_id_child on parent2child(id,child);
END

	  meta => <<END,
(
  name  text primary key,
  value text not null
)
END
	  sequence => <<END,
(
  id       integer not null,
  offset   integer not null,
  sequence blob,
  primary key(id,offset)
)
END
	 };

  unless ($self->_has_spatial_index) {
    $defs->{feature_location} = <<END;
(
  id       int(10) primary key,
  seqid    int(10),
  bin      int,
  start    int,
  end      int
);
create index index_feature_location on feature_location(seqid,bin,start,end);
END

    }

  if (EXPERIMENTAL_COVERAGE) {
    $defs->{interval_stats} = <<END;
(
   typeid            integer not null,
   seqid             integer not null,
   bin               integer not null,
   cum_count         integer not null,
   unique(typeid,seqid,bin)
);
END
  }
  return $defs;
}

sub _init_database {
    my $self  = shift;

    # must do this first before calling table_definitions
    $self->_create_spatial_index;
    $self->SUPER::_init_database(@_);
}

sub init_tmp_database {
    my $self = shift;
    my $erase = shift;
    $self->_create_spatial_index;
    $self->SUPER::init_tmp_database(@_);
}

sub _create_spatial_index{
    my $self = shift;
    my $dbh   = $self->dbh;
    local $dbh->{PrintError} = 0;
    $dbh->do("DROP TABLE IF EXISTS feature_index"); # spatial index
    if (USE_SPATIAL) {
	$dbh->do("CREATE VIRTUAL TABLE feature_index USING RTREE(id,seqid,bin,start,end)");
    }
}

sub _has_spatial_index {
    my $self = shift;
    return $self->{'_has_spatial_index'} if exists $self->{'_has_spatial_index'};
    my $dbh  = $self->dbh;
    my ($count) = $dbh->selectrow_array("select count(*) from sqlite_master where name='feature_index'");
    return $self->{'_has_spatial_index'} = $count;
}

sub _finish_bulk_update {
  my $self = shift;
  my $dbh  = $self->dbh;
  my $dir = $self->{dumpdir} || '.';

  $self->begin_work; # making this a transaction greatly improves performance

  for my $table ('feature', $self->index_tables) {
    my $fh = $self->dump_filehandle($table);
    my $path = $self->dump_path($table);
    $fh->close;

    open $fh, '<', $path or $self->throw("Could not read file '$path': $!");
    my $qualified_table = $self->_qualify($table);

    my $sth;
    if ($table =~ /feature$/) {
      $sth = $dbh->prepare("REPLACE INTO $qualified_table VALUES (?,?,?,?,?)");

      while (<$fh>) {
        chomp();
        my ($id,$typeid,$strand,$indexed,$obj) = split(/\t/);
        $sth->bind_param(1, $id);
        $sth->bind_param(2, $typeid);
        $sth->bind_param(3, $strand);
        $sth->bind_param(4, $indexed);
        $sth->bind_param(5, pack('H*',$obj), {TYPE => SQL_BLOB});
        $sth->execute();
      }
    } else {
	my $feature_index = $self->_feature_index_table;
	if ($table =~ /parent2child$/) {
	    $sth = $dbh->prepare("REPLACE INTO $qualified_table VALUES (?,?)");
	} elsif ($table =~ /$feature_index$/) {
	    $sth = $dbh->prepare(
		$self->_has_spatial_index ?"REPLACE INTO $qualified_table VALUES (?,?,?,?,?)"
		                          :"REPLACE INTO $qualified_table (id,seqid,bin,start,end) VALUES (?,?,?,?,?)"
		);
	} else { # attribute or name
	    $sth = $dbh->prepare("REPLACE INTO $qualified_table VALUES (?,?,?)");
	}

	while (<$fh>) {
	    chomp();
	    $sth->execute(split(/\t/));
	}
    }
    $fh->close();
    unlink $path;
  }
  $self->commit; # commit the transaction
  delete $self->{bulk_update_in_progress};
  delete $self->{filehandles};
}

sub index_tables {
    my $self = shift;
    my @t    = $self->SUPER::index_tables;
    return (@t,$self->_feature_index_table);
}

sub _enable_keys  { }  # nullop
sub _disable_keys { }  # nullop

sub _fetch_indexed_features_sql {
    my $self = shift;
    my $location_table       = $self->_qualify('feature_location');
    my $feature_table        = $self->_qualify('feature');
    return <<END;
    SELECT typeid,seqid,start-1,end
  FROM $location_table as l,$feature_table as f
 WHERE l.id=f.id AND f.\"indexed\"=1
  ORDER BY typeid,seqid,start
END
}

###
# get primary sequence between start and end
#
sub _fetch_sequence {
  my $self = shift;
  my ($seqid,$start,$end) = @_;

  # backward compatibility to the old days when I liked reverse complementing
  # dna by specifying $start > $end
  my $reversed;
  if (defined $start && defined $end && $start > $end) {
    $reversed++;
    ($start,$end) = ($end,$start);
  }
  $start-- if defined $start;
  $end--   if defined $end;

  my $offset1 = $self->_offset_boundary($seqid,$start || 'left');
  my $offset2 = $self->_offset_boundary($seqid,$end   || 'right');
  my $sequence_table = $self->_sequence_table;
  my $locationlist_table = $self->_locationlist_table;

  # CROSS JOIN gives a hint to the SQLite query optimizer -- mucho speedup!
  my $sth     = $self->_prepare(<<END);
SELECT sequence,offset
   FROM $locationlist_table as ll CROSS JOIN $sequence_table as s
   WHERE ll.id=s.id
     AND ll.seqname= ?
     AND offset >= ?
     AND offset <= ?
   ORDER BY offset
END

  my $seq = '';
  $sth->execute($seqid,$offset1,$offset2) or $self->throw($sth->errstr);

  while (my($frag,$offset) = $sth->fetchrow_array) {
    substr($frag,0,$start-$offset) = '' if defined $start && $start > $offset;
    $seq .= $frag;
  }
  substr($seq,$end-$start+1) = '' if defined $end && $end-$start+1 < length($seq);
  if ($reversed) {
    $seq = reverse $seq;
    $seq =~ tr/gatcGATC/ctagCTAG/;
  }
  $sth->finish;
  $seq;
}

sub _offset_boundary {
  my $self = shift;
  my ($seqid,$position) = @_;

  my $sequence_table     = $self->_sequence_table;
  my $locationlist_table = $self->_locationlist_table;

  my $sql;
  # use "CROSS JOIN" to give a hint to the SQLite query optimizer.
  $sql =  $position eq 'left'  ? "SELECT min(offset) FROM $locationlist_table as ll CROSS JOIN $sequence_table as s ON ll.id=s.id WHERE ll.seqname=?"
         :$position eq 'right' ? "SELECT max(offset) FROM $locationlist_table as ll CROSS JOIN $sequence_table as s ON ll.id=s.id WHERE ll.seqname=?"
	 :"SELECT max(offset) FROM $locationlist_table as ll CROSS JOIN $sequence_table as s ON ll.id=s.id WHERE ll.seqname=? AND offset<=?";
  my $sth = $self->_prepare($sql);
  my @args = $position =~ /^-?\d+$/ ? ($seqid,$position) : ($seqid);
  $sth->execute(@args) or $self->throw($sth->errstr);
  my $boundary = $sth->fetchall_arrayref->[0][0];
  $sth->finish;
  return $boundary;
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
  my ($seq_id,$start,$end,$strand,
      $name,$class,$allow_aliases,
      $types,
      $attributes,
      $range_type,
      $fromtable,
      $iterator,
      $sources
     ) = rearrange([['SEQID','SEQ_ID','REF'],'START',['STOP','END'],'STRAND',
		    'NAME','CLASS','ALIASES',
		    ['TYPES','TYPE','PRIMARY_TAG'],
		    ['ATTRIBUTES','ATTRIBUTE'],
		    'RANGE_TYPE',
		    'FROM_TABLE',
		    'ITERATOR',
		    ['SOURCE','SOURCES']
		   ],@_);

  my (@from,@where,@args,@group);
  $range_type ||= 'overlaps';

  my $feature_table         = $self->_feature_table;
  @from = "$feature_table as f";

  if (defined $name) {
    # hacky backward compatibility workaround
    undef $class if $class && $class eq 'Sequence';
    $name = "$class:$name" if defined $class && length $class > 0;
    # last argument is the join field
    my ($from,$where,$group,@a) = $self->_name_sql($name,$allow_aliases,'f.id');
    push @from,$from   if $from;
    push @where,$where if $where;
    push @group,$group if $group;
    push @args,@a;
  }

  if (defined $seq_id) {
    # last argument is the name of the features table
    my ($from,$where,$group,@a) = $self->_location_sql($seq_id,$start,$end,$range_type,$strand,'f');
    push @from,$from   if $from;
    push @where,$where if $where;
    push @group,$group if $group;
    push @args,@a;
  }

  if (defined($sources)) {
    my @sources = ref($sources) eq 'ARRAY' ? @{$sources} : ($sources);
    if (defined($types)) {
        my @types = ref($types) eq 'ARRAY' ? @{$types} : ($types);
        my @final_types;
        foreach my $type (@types) {
            # *** not sure what to do if user supplies both -source and -type
            #     where the type includes a source!
            if ($type =~ /:/) {
                push(@final_types, $type);
            }
            else {
                foreach my $source (@sources) {
                    push(@final_types, $type.':'.$source);
                }
            }
        }
        $types = \@final_types;
    }
    else {
        $types = [map { ':'.$_ } @sources];
    }
  }
  if (defined($types)) {
    # last argument is the name of the features table
    my ($from,$where,$group,@a) = $self->_types_sql($types,'f');
    push @from,$from   if $from;
    push @where,$where if $where;
    push @group,$group if $group;
    push @args,@a;
  }

  if (defined $attributes) {
    # last argument is the join field
    my ($from,$where,$group,@a) = $self->_attributes_sql($attributes,'f.id');
    push @from,$from    if $from;
    push @where,$where  if $where;
    push @group,$group  if $group;
    push @args,@a;
  }

  if (defined $fromtable) {
    # last argument is the join field
    my ($from,$where,$group,@a) = $self->_from_table_sql($fromtable,'f.id');
    push @from,$from    if $from;
    push @where,$where  if $where;
    push @group,$group  if $group;
    push @args,@a;
  }

  # if no other criteria are specified, then
  # only fetch indexed (i.e. top level objects)
  @where = '"indexed"=1' unless @where;

  my $from  = join ', ',@from;
  my $where = join ' AND ',map {"($_)"} @where;
  my $group = join ', ',@group;
  $group    = "GROUP BY $group" if @group;

  my $query = <<END;
SELECT f.id,f.object
  FROM $from
  WHERE $where
  $group
END

  $self->_print_query($query,@args) if DEBUG || $self->debug;

  my $sth = $self->_prepare($query);
  $sth->execute(@args) or $self->throw($sth->errstr);
  return $iterator ? Bio::DB::SeqFeature::Store::DBI::Iterator->new($sth,$self) : $self->_sth2objs($sth);
}

sub _make_attribute_group {
  my $self                     = shift;
  my ($table_name,$attributes) = @_;
  my $key_count = keys %$attributes or return;
  my $count = $key_count-1;
  return "f.id HAVING count(f.id)>$count";
}

sub _location_sql {
  my $self = shift;
  my ($seq_id,$start,$end,$range_type,$strand,$location) = @_;

  # the additional join on the location_list table badly impacts performance
  # so we build a copy of the table in memory
  my $seqid = $self->_locationid_nocreate($seq_id) || 0; # zero is an invalid primary ID, so will return empty

  my $feature_index = $self->_feature_index_table;
  my $from  = "$feature_index as fi";

  my ($bin_where,@bin_args);
  if (defined $start && defined $end && !$self->_has_spatial_index) {
      my @bins   = $self->search_bins($start,$end);
      $bin_where = ' AND bin in ('.join(',',@bins).')';
  }

  $start = MIN_INT unless defined $start;
  $end   = MAX_INT unless defined $end;

  my ($range,@range_args);
  if ($range_type eq 'overlaps') {
      $range = "fi.end>=? AND fi.start<=?".$bin_where;
      @range_args = ($start,$end,@bin_args);
  } elsif ($range_type eq 'contains') {
      $range = "fi.start>=? AND fi.end<=?".$bin_where;
      @range_args = ($start,$end,@bin_args);
  } elsif ($range_type eq 'contained_in') {
      $range = "fi.start<=? AND fi.end>=?";
      @range_args = ($start,$end);
  } else {
      $self->throw("range_type must be one of 'overlaps', 'contains' or 'contained_in'");
  }

  if (defined $strand) {
      $range .= " AND strand=?";
      push @range_args,$strand;
  }

  my $where = <<END;
  fi.seqid=?
      AND   $location.id=fi.id
      AND   $range
END
;
  my $group = '';

  my @args  = ($seqid,@range_args);
  return ($from,$where,$group,@args);
}

sub _feature_index_table          {
    my $self = shift;
    return $self->_has_spatial_index ? $self->_qualify('feature_index')
                                     : $self->_qualify('feature_location') }

# Do a case-insensitive search a la the PostgreSQL adaptor
sub _name_sql {
  my $self = shift;
  my ($name,$allow_aliases,$join) = @_;
  my $name_table   = $self->_name_table;

  my $from  = "$name_table as n";
  my ($match,$string) = $self->_match_sql($name);

  my $where = "n.id=$join AND n.name $match COLLATE NOCASE";
  $where   .= " AND n.display_name>0" unless $allow_aliases;
  return ($from,$where,'',$string);
}

sub _search_attributes {
  my $self = shift;
  my ($search_string,$attribute_names,$limit) = @_;
  my @words               = map {quotemeta($_)} split /\s+/,$search_string;

  my $name_table          = $self->_name_table;
  my $attribute_table     = $self->_attribute_table;
  my $attributelist_table = $self->_attributelist_table;
  my $type_table          = $self->_type_table;
  my $typelist_table      = $self->_typelist_table;

  my @tags    = @$attribute_names;
  my $tag_sql = join ' OR ',("al.tag=?") x @tags;

  my $perl_regexp = join '|',@words;

  my @wild_card_words = map { "%$_%" } @words;
  my $sql_regexp = join ' OR ',("a.attribute_value LIKE ?")  x @words;
  # CROSS JOIN disables SQLite's table reordering optimization
  my $sql = <<END;
SELECT name,attribute_value,tl.tag,n.id
  FROM $attributelist_table        AS al
       CROSS JOIN $attribute_table AS a  ON al.id = a.attribute_id
       CROSS JOIN $name_table      AS n  ON n.id = a.id
       CROSS JOIN $type_table      AS t  ON t.id = n.id
       CROSS JOIN $typelist_table  AS tl ON tl.id = t.typeid
  WHERE ($tag_sql)
    AND ($sql_regexp)
    AND n.display_name=1
END
  $sql .= "LIMIT $limit" if defined $limit;
  $self->_print_query($sql,@tags,@words) if DEBUG || $self->debug;
  my $sth = $self->_prepare($sql);
  $sth->execute(@tags,@wild_card_words) or $self->throw($sth->errstr);

  my @results;
  while (my($name,$value,$type,$id) = $sth->fetchrow_array) {
    my (@hits) = $value =~ /$perl_regexp/ig;
    my @words_in_row = split /\b/,$value;
    my $score  = int(@hits*100/@words/@words_in_row);
    push @results,[$name,$value,$score,$type,$id];
  }
  $sth->finish;
  @results = sort {$b->[2]<=>$a->[2]} @results;
  return @results;
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
    $match = "= ? COLLATE NOCASE";
    $string  = $name;
  }
  return ($match,$string);
}

sub _attributes_sql {
  my $self = shift;
  my ($attributes,$join) = @_;

  my ($wf,@bind_args)       = $self->_make_attribute_where('a','al',$attributes);
  my ($group_by,@group_args)= $self->_make_attribute_group('a',$attributes);

  my $attribute_table       = $self->_attribute_table;
  my $attributelist_table   = $self->_attributelist_table;

  my $from = "$attribute_table AS a INDEXED BY index_attribute_id, $attributelist_table AS al";

  my $where = <<END;
  a.id=$join
  AND   a.attribute_id=al.id
  AND ($wf)
END

  my $group = $group_by;

  my @args  = (@bind_args,@group_args);
  return ($from,$where,$group,@args);
}

# overridden because of case-sensitivity of matches
sub _types_sql {
  my $self  = shift;
  my ($types,$type_table) = @_;
  my ($primary_tag,$source_tag);

  my @types = ref $types eq 'ARRAY' ?  @$types : $types;

  my $typelist      = $self->_typelist_table;
  my $from = "$typelist AS tl";

  my (@matches,@args);

  for my $type (@types) {

    if (ref $type && $type->isa('Bio::DB::GFF::Typename')) {
      $primary_tag = $type->method;
      $source_tag  = $type->source;
    } else {
      ($primary_tag,$source_tag) = split ':',$type,2;
    }

    if (length $source_tag) {
      push @matches,"tl.tag=? COLLATE NOCASE";
      push @args,"$primary_tag:$source_tag";
    } else {
      push @matches,"tl.tag LIKE ?";
      push @args,"$primary_tag:%";
    }
  }
  my $matches = join ' OR ',@matches;

  my $where = <<END;
   tl.id=$type_table.typeid
   AND   ($matches)
END

  return ($from,$where,'',@args);
}

sub optimize {
  my $self = shift;
  $self->dbh->do("ANALYZE $_") foreach $self->index_tables;
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
REPLACE INTO $features (id,object,"indexed",strand,typeid) VALUES (?,?,?,?,?)
END

  my ($seqid,$start,$end,$strand,$bin) = $index_flag ? $self->_get_location_and_bin($object) : (undef)x6;

  my $primary_tag = $object->primary_tag;
  my $source_tag  = $object->source_tag || '';
  $primary_tag    .= ":$source_tag";
  my $typeid   = $self->_typeid($primary_tag,1);

  my $frozen = $self->no_blobs() ? 0 : $self->freeze($object);

  $sth->bind_param(1, $id);
  $sth->bind_param(2, $frozen, {TYPE => SQL_BLOB});
  $sth->bind_param(3, $index_flag||0);
  $sth->bind_param(4, $strand);
  $sth->bind_param(5, $typeid);

  $sth->execute() or $self->throw($sth->errstr);

  my $dbh = $self->dbh;
  $object->primary_id($dbh->func('last_insert_rowid')) unless defined $id;

  $self->flag_for_indexing($dbh->func('last_insert_rowid')) if $self->{bulk_update_in_progress};
}

# doesn't work with this schema, since we have to update name and attribute
# tables which need object ids, which we can only know by replacing feats in
# the feature table one by one
sub bulk_replace {
    my $self       = shift;
    my $index_flag = shift || undef;
    my @objects    = @_;

    my $features = $self->_feature_table;

    my @insert_values;
    foreach my $object (@objects) {
        my $id = $object->primary_id;
        my (undef,undef,undef,$strand) = $index_flag ? $self->_get_location_and_bin($object) : (undef)x4;
        my $primary_tag = $object->primary_tag;
        my $source_tag  = $object->source_tag || '';
        $primary_tag    .= ":$source_tag";
        my $typeid   = $self->_typeid($primary_tag,1);

        push(@insert_values, ($id,0,$index_flag||0,$strand,$typeid));
    }

    my @value_blocks;
    for (1..@objects) {
        push(@value_blocks, '(?,?,?,?,?)');
    }
    my $value_blocks = join(',', @value_blocks);
    my $sql = qq{REPLACE INTO $features (id,object,"indexed",strand,typeid) VALUES $value_blocks};

    my $sth = $self->_prepare($sql);
    $sth->execute(@insert_values) or $self->throw($sth->errstr);
}

sub _get_location_and_bin {
    my $self = shift;
    my $obj  = shift;
    my $seqid   = $self->_locationid($obj->seq_id||'');
    my $start   = $obj->start;
    my $end     = $obj->end;
    my $strand  = $obj->strand;
    return ($seqid,$start,$end,$strand,$self->calculate_bin($start,$end));
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
INSERT INTO $features (id,object,"indexed") VALUES (?,?,?)
END
  $sth->execute(undef,$self->freeze($object),$index_flag) or $self->throw($sth->errstr);
  my $dbh = $self->dbh;
  $object->primary_id($dbh->func('last_insert_rowid'));
  $self->flag_for_indexing($dbh->func('last_insert_rowid')) if $self->{bulk_update_in_progress};
}

=head2 toplevel_types

 Title   : toplevel_types
 Usage   : @type_list = $db->toplevel_types
 Function: Get the toplevel types in the database
 Returns : array of Bio::DB::GFF::Typename objects
 Args    : none
 Status  : public

This is similar to types() but only returns the types of
INDEXED (toplevel) features.

=cut

sub toplevel_types {
    my $self = shift;
    eval "require Bio::DB::GFF::Typename"
	unless Bio::DB::GFF::Typename->can('new');
    my $typelist_table      = $self->_typelist_table;
    my $feature_table       = $self->_feature_table;
    my $sql = <<END;
SELECT distinct(tag) from $typelist_table as tl,$feature_table as f
 WHERE tl.id=f.typeid
   AND f."indexed"=1
END
;
    $self->_print_query($sql) if DEBUG || $self->debug;
    my $sth = $self->_prepare($sql);
    $sth->execute() or $self->throw($sth->errstr);

    my @results;
    while (my($tag) = $sth->fetchrow_array) {
	push @results,Bio::DB::GFF::Typename->new($tag);
    }
    $sth->finish;
    return @results;
}

sub _genericid {
  my $self = shift;
  my ($table,$namefield,$name,$add_if_missing) = @_;
  my $qualified_table = $self->_qualify($table);
  my $sth = $self->_prepare(<<END);
SELECT id FROM $qualified_table WHERE $namefield=? COLLATE NOCASE
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
  return $dbh->func('last_insert_rowid');
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
    my ($seqid,$start,$end,$strand) = $indexed ? $self->_get_location_and_bin($obj) : (undef)x4;
    my $primary_tag = $obj->primary_tag;
    my $source_tag  = $obj->source_tag || '';
    $primary_tag    .= ":$source_tag";
    my $typeid   = $self->_typeid($primary_tag,1);

    # Encode BLOB in hex so we can more easily import it into SQLite
    print $store_fh
    join("\t",$id,$typeid,$strand,$indexed,
         unpack('H*', $self->freeze($obj))),"\n";
    $obj->primary_id($id);
    $self->_update_indexes($obj) if $indexed && $autoindex;
    $count++;
  }

  # remember whether we are have ever stored a non-indexed feature
  unless ($indexed or $self->{indexed_flag}++) {
    $self->subfeatures_are_indexed(0);
  }
  $count;
}

sub _dump_update_name_index {
  my $self = shift;
  my ($obj,$id) = @_;
  my $fh      = $self->dump_filehandle('name');
  my $dbh     = $self->dbh;
  my ($names,$aliases) = $self->feature_names($obj);
  # unlike DBI::mysql, don't quote, as quotes will be quoted when loaded
  print $fh join("\t",$id,$_,1),"\n" foreach @$names;
  print $fh join("\t",$id,$_,0),"\n" foreach @$aliases;
}

sub _update_name_index {
  my $self = shift;
  my ($obj,$id) = @_;
  my $name = $self->_name_table;
  my $primary_id = $obj->primary_id;

  $self->_delete_index($name,$id);
  my ($names,$aliases) = $self->feature_names($obj);

  my $sth = $self->_prepare("INSERT INTO $name (id,name,display_name) VALUES (?,?,?)");

  $sth->execute($id,$_,1) or $self->throw($sth->errstr)   foreach @$names;
  $sth->execute($id,$_,0) or $self->throw($sth->errstr) foreach @$aliases;
  $sth->finish;
}


sub _dump_update_attribute_index {
  my $self = shift;
  my ($obj,$id) = @_;
  my $fh        = $self->dump_filehandle('attribute');
  my $dbh       = $self->dbh;
  for my $tag ($obj->all_tags) {
    my $tagid = $self->_attributeid($tag);
    for my $value ($obj->each_tag_value($tag)) {
      # unlike DBI::mysql, don't quote, as quotes will be quoted when loaded
      print $fh join("\t",$id,$tagid,$value),"\n";
    }
  }
}

sub _update_indexes {
    my $self = shift;
    my $obj  = shift;
    defined (my $id   = $obj->primary_id) or return;
    $self->SUPER::_update_indexes($obj);

    if ($self->{bulk_update_in_progress}) {
	$self->_dump_update_location_index($obj,$id);
    } else {
	$self->_update_location_index($obj,$id);
    }
}

sub _update_location_index {
    my $self = shift;
    my ($obj,$id) = @_;
    my ($seqid,$start,$end,$strand,$bin) = $self->_get_location_and_bin($obj);

    my $table = $self->_feature_index_table;
    $self->_delete_index($table,$id);

    my ($sql,@args);

    if ($self->_has_spatial_index) {
	$sql    = "INSERT INTO $table (id,seqid,bin,start,end) values (?,?,?,?,?)";
	@args   = ($id,$seqid,$bin,$start,$end);
    } else {
	$sql    = "INSERT INTO $table (id,seqid,bin,start,end) values (?,?,?,?,?)";
	@args   = ($id,$seqid,$bin,$start,$end);
    }

    my $sth  = $self->_prepare($sql);
    $sth->execute(@args);
    $sth->finish;
}

sub _dump_update_location_index {
    my $self = shift;
    my ($obj,$id) = @_;
    my $table   = $self->_feature_index_table;
    my $fh      = $self->dump_filehandle($table);
    my $dbh     = $self->dbh;
    my ($seqid,$start,$end,$strand,$bin) = $self->_get_location_and_bin($obj);
    my @args = $self->_has_spatial_index ? ($id,$seqid,$bin,$start,$end)
	                                 : ($id,$seqid,$bin,$start,$end);
    print $fh join("\t",@args),"\n";
}

sub DESTROY {
    my $self = shift;
    # Remove filehandles, so temporal files can be properly deleted
    if (%DBI::installed_drh) {
	DBI->disconnect_all;
	%DBI::installed_drh = ();
    }
    undef $self->{dbh};
}

1;

=head1 AUTHOR

Nathan Weeks - Nathan.Weeks@ars.usda.gov

Copyright (c) 2009 Nathan Weeks

Modified 2010 to support cumulative statistics by Lincoln Stein
<lincoln.stein@gmail.com>.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself. See the Bioperl license for
more details.

=cut
