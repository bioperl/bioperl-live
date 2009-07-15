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
use Cwd 'abs_path';
use Bio::DB::GFF::Util::Rearrange 'rearrange';
use Bio::SeqFeature::Lite;
use File::Spec;
use constant DEBUG=>0;

# Using same limits as MySQL adaptor so I don't have to make something up.
use constant MAX_INT =>  2_147_483_647;
use constant MIN_INT => -2_147_483_648;
use constant MAX_BIN =>  1_000_000_000;  # size of largest feature = 1 Gb
use constant MIN_BIN =>  1000;           # smallest bin we'll make - on a 100 Mb chromosome, there'll be 100,000 of these

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
  return {
	  feature => <<END,
(
  id        integer primary key autoincrement,
  typeid    integer not null,
  seqid     integer,
  start     integer,
  end       integer,
  strand    integer default 0,
  tier      integer,
  bin       integer,
  "indexed" integer default 1,
  object    blob not null
);
create index index_feature_seqid_tier_bin_typeid on feature (seqid,tier,bin,typeid);
create index index_feature_typeid on feature(typeid);
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
  tag text    not null
);
create index index_typelist on typelist (tag);
END
	  name => <<END,
(
  id           integer not null,
  name         text    not null,
  display_name integer default 0
);
create index index_name_id on name(id);
create index index_name_name on name(name);
END

	  attribute => <<END,
(
  id              integer not null,
  attribute_id    integer not null,
  attribute_value text
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
create index index_parent2child_id_child on parent2child(id,child);
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
}

sub _finish_bulk_update {
  my $self = shift;
  my $dbh  = $self->dbh;
  my $dir = $self->{dumpdir} || '.';

  $dbh->begin_work; # making this a transaction greatly improves performance
  
  for my $table ('feature', $self->index_tables) {
    my $fh = $self->dump_filehandle($table);
    my $path = $self->dump_path($table);
    $fh->close;
    open($fh, $path);
    my $qualified_table = $self->_qualify($table);

    my $sth;
    if ($table eq 'feature') {
      $sth = $dbh->prepare("REPLACE INTO $qualified_table VALUES (?,?,?,?,?,?,?,?,?,?)");

      while (<$fh>) {
        chomp();
        my ($id,$typeid,$seqid,$start,$end,$strand,$tier,$bin,$indexed,$obj) = 
            split(/\t/);
        $sth->bind_param(1, $id);
        $sth->bind_param(2, $typeid);
        $sth->bind_param(3, $seqid);
        $sth->bind_param(4, $start);
        $sth->bind_param(5, $end);
        $sth->bind_param(6, $strand);
        $sth->bind_param(7, $tier);
        $sth->bind_param(8, $bin);
        $sth->bind_param(9, $indexed);
        $sth->bind_param(10, pack('H*',$obj), {TYPE => SQL_BLOB});
        $sth->execute();
      }
    } else {
      if ($table eq 'parent2child') {
        $sth = $dbh->prepare("REPLACE INTO $qualified_table VALUES (?,?)");
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
  $dbh->commit; # commit the transaction
  delete $self->{bulk_update_in_progress};
  delete $self->{filehandles};
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
SELECT f.id,f.object,f.typeid,f.seqid,f.start,f.end,f.strand
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

# Do a case-insensitive search a la the PostgreSQL adaptor
sub _name_sql {
  my $self = shift;
  my ($name,$allow_aliases,$join) = @_;
  my $name_table   = $self->_name_table;

  my $from  = "$name_table as n";
  my ($match,$string) = $self->_match_sql($name);

  my $where = "n.id=$join AND lower(n.name) $match";
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
       CROSS JOIN $attribute_table AS a ON al.id = a.attribute_id
       CROSS JOIN $name_table      AS n ON n.id = a.id
       CROSS JOIN $type_table      AS t ON t.id = n.id
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
    $match = "= lower(?)";
    $string  = lc($name);
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

    if (defined $source_tag) {
      push @matches,"lower(tl.tag)=lower(?)";
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
REPLACE INTO $features (id,object,"indexed",seqid,start,end,strand,tier,bin,typeid) VALUES (?,?,?,?,?,?,?,?,?,?)
END

  my @location = $index_flag ? $self->_get_location_and_bin($object) : (undef)x6;

  my $primary_tag = $object->primary_tag;
  my $source_tag  = $object->source_tag || '';
  $primary_tag    .= ":$source_tag";
  my $typeid   = $self->_typeid($primary_tag,1);

  my $frozen = $self->no_blobs() ? 0 : $self->freeze($object);

  $sth->bind_param(1, $id);
  $sth->bind_param(2, $frozen, {TYPE => SQL_BLOB});
  $sth->bind_param(3, $index_flag||0);
  $sth->bind_param(4, $location[0]);
  $sth->bind_param(5, $location[1]);
  $sth->bind_param(6, $location[2]);
  $sth->bind_param(7, $location[3]);
  $sth->bind_param(8, $location[4]);
  $sth->bind_param(9, $location[5]);
  $sth->bind_param(10,$typeid);

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
        my @location = $index_flag ? $self->_get_location_and_bin($object) : (undef)x6;
        my $primary_tag = $object->primary_tag;
        my $source_tag  = $object->source_tag || '';
        $primary_tag    .= ":$source_tag";
        my $typeid   = $self->_typeid($primary_tag,1);
        
        push(@insert_values, ($id,0,$index_flag||0,@location,$typeid));
    }
    
    my @value_blocks;
    for (1..@objects) {
        push(@value_blocks, '(?,?,?,?,?,?,?,?,?,?)');
    }
    my $value_blocks = join(',', @value_blocks);
    my $sql = qq{REPLACE INTO $features (id,object,"indexed",seqid,start,end,strand,tier,bin,typeid) VALUES $value_blocks};
    
    my $sth = $self->_prepare($sql);
    $sth->execute(@insert_values) or $self->throw($sth->errstr);
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

=head2 types

 Title   : types
 Usage   : @type_list = $db->types
 Function: Get all the types in the database
 Returns : array of Bio::DB::GFF::Typename objects
 Args    : none
 Status  : public

=cut

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
    my ($seqid,$start,$end,$strand,$tier,$bin) = $indexed ? $self->_get_location_and_bin($obj) : (undef)x6;
    my $primary_tag = $obj->primary_tag;
    my $source_tag  = $obj->source_tag || '';
    $primary_tag    .= ":$source_tag";
    my $typeid   = $self->_typeid($primary_tag,1);

    # Encode BLOB in hex so we can more easily import it into SQLite
    print $store_fh
    join("\t",$id,$typeid,$seqid,$start,$end,$strand,$tier,$bin,$indexed,
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


1;

=head1 AUTHOR

Nathan Weeks <Nathan.Weeks@ars.usda.gov>

Copyright (c) 2009 Nathan Weeks

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

