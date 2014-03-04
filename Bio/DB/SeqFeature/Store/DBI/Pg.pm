
package Bio::DB::SeqFeature::Store::DBI::Pg;
use DBD::Pg qw(:pg_types);
use MIME::Base64;
# $Id: Pg.pm 14656 2008-04-14 15:05:37Z lstein $

=head1 NAME

Bio::DB::SeqFeature::Store::DBI::Pg -- PostgreSQL implementation of Bio::DB::SeqFeature::Store

=head1 SYNOPSIS

  use Bio::DB::SeqFeature::Store;

  # Open the sequence database
  my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::Pg',
                                          -dsn     => 'dbi:Pg:test');

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

Bio::DB::SeqFeature::Store::Pg is the Pg adaptor for
Bio::DB::SeqFeature::Store. You will not create it directly, but
instead use Bio::DB::SeqFeature::Store-E<gt>new() to do so.

See L<Bio::DB::SeqFeature::Store> for complete usage instructions.

=head2 Using the Pg adaptor

Before you can use the adaptor, you must use the Pgadmin tool to
create a database and establish a user account with write
permission. In order to use "fast" loading, the user account must have
"file" privileges.

To establish a connection to the database, call
Bio::DB::SeqFeature::Store-E<gt>new(-adaptor=E<gt>'DBI::Pg',@more_args). The
additional arguments are as follows:

  Argument name       Description
  -------------       -----------

 -dsn              The database name. You can abbreviate
                   "dbi:Pg:foo" as "foo" if you wish.

 -user             Username for authentication.

 -pass             Password for authentication.

 -namespace        Creates a SCHEMA for the tables. This allows you
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
Bio::DB::SeqFeature::Store::DBI::Pg will be returned.

In addition to the standard methods supported by all well-behaved
Bio::DB::SeqFeature::Store databases, several following
adaptor-specific methods are provided. These are described in the next
sections.

=cut

use strict;

use base 'Bio::DB::SeqFeature::Store::DBI::mysql';
use Bio::DB::SeqFeature::Store::DBI::Iterator;
use DBI;
use Memoize;
use Cwd 'abs_path';
use Bio::DB::GFF::Util::Rearrange 'rearrange';
use File::Copy;
use File::Spec;
use constant DEBUG=>0;

use constant MAX_INT =>  2_147_483_647;
use constant MIN_INT => -2_147_483_648;
use constant MAX_BIN =>  1_000_000_000;  # size of largest feature = 1 Gb
use constant MIN_BIN =>  1000;           # smallest bin we'll make - on a 100 Mb chromosome, there'll be 100,000 of these

###
# object initialization
#
# NOTE: most of this code can be refactored and inherited from DBI or DBI::mysql adapter
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
      $schema,
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
                    'SCHEMA'
		   ],@_);


  $dbi_options  ||= {pg_server_prepare => 0};
  $writeable    = 1 if $is_temporary or $dump_dir;

  $dsn or $self->throw("Usage: ".__PACKAGE__."->init(-dsn => \$dbh || \$dsn)");

  my $dbh;
  if (ref $dsn) {
    $dbh = $dsn;
  } else {
    $dsn = "dbi:Pg:$dsn" unless $dsn =~ /^dbi:/;
    $dbh = DBI->connect($dsn,$user,$pass,$dbi_options) or $self->throw($DBI::errstr);
  }
  $dbh->do('set client_min_messages=warning') if $dbh;
  $self->{'original_arguments'} = {
      'dsn' => $dsn,
      'user' => $user,
      'pass' => $pass,
      'dbh_options' => $dbi_options,
  };
  $self->{dbh}       = $dbh;
  $self->{dbh}->{InactiveDestroy} = 1;
  $self->{is_temp}   = $is_temporary;
  $self->{writeable} = $writeable;
  $self->{namespace} = $namespace || $schema || 'public';
  $self->schema($self->{namespace});

  $self->default_settings;
  $self->autoindex($autoindex)                   if defined $autoindex;
  $self->dumpdir($dump_dir)                      if $dump_dir;
  if ($self->is_temp) {
    # warn "creating a temp database isn't supported";
    #$self->init_tmp_database();
    $self->init_database('erase');
  } elsif ($create) {
    $self->init_database('erase');
  }
}

sub table_definitions {
  my $self = shift;
  return {
	  feature => <<END,
(
  id       serial primary key,
  typeid   int      not null,
  seqid    int,
  start    int,
  "end"      int,
  strand   int      default 0,
  tier     int,
  bin      int,
  indexed  int default 1,
  object     bytea not null
);
  CREATE INDEX feature_stuff ON feature(seqid,tier,bin,typeid);
  CREATE INDEX feature_typeid ON feature(typeid);
END

	  locationlist => <<END,
(
  id         serial primary key,
  seqname    text   not null
); CREATE INDEX locationlist_seqname ON locationlist(seqname);
END

	  typelist => <<END,
(
  id       serial primary key,
  tag      text  not null
); CREATE INDEX typelist_tab ON typelist(tag);
END
	  name => <<END,
(
  id           int       not null,
  name         text  not null,
  display_name int       default 0
);
  CREATE INDEX name_id ON name( id );
  CREATE INDEX name_name ON name( name );
  CREATE INDEX name_lcname ON name( lower(name) );
  CREATE INDEX name_lcname_varchar_patt_ops ON name USING BTREE (lower(name) varchar_pattern_ops);
END

	  attribute => <<END,
(
  id               int       not null,
  attribute_id     int   not null,
  attribute_value  text
);
  CREATE INDEX attribute_id ON attribute(id);
  CREATE INDEX attribute_id_val ON attribute(attribute_id,SUBSTR(attribute_value, 1, 10));
END

	  attributelist => <<END,
(
  id       serial primary key,
  tag      text  not null
);
  CREATE INDEX attributelist_tag ON attributelist(tag);
END
	  parent2child => <<END,
(
  id               int       not null,
  child            int       not null
);
  CREATE UNIQUE INDEX parent2child_id_child ON parent2child(id,child);
END

	  meta => <<END,
(
  name      text primary key,
  value     text not null
)
END
	  sequence => <<END,
(
  id       int not null,
  "offset"   int not null,
  sequence text,
  primary key(id,"offset")
)
END

    interval_stats => <<END,
(
   typeid            int not null,
   seqid             int not null,
   bin               int not null,
   cum_count         int not null
 );
 CREATE UNIQUE INDEX interval_stats_id ON interval_stats(typeid,seqid,bin);
END
	 };
}

sub schema {
  my ($self, $schema) = @_;
  $self->{'schema'} = $schema if defined($schema);
  if ($schema) {
    $self->_check_for_namespace();
    $self->dbh->do("SET search_path TO " . $self->{'schema'} );
  } else {
    $self->dbh->do("SET search_path TO public");
  }
  return $self->{'schema'};
}
###
# wipe database clean and reinstall schema
#
sub _init_database {
  my $self = shift;
  my $erase = shift;

  my $dbh    = $self->dbh;
  my $namespace = $self->namespace;
  my $tables = $self->table_definitions;
    my $temporary = $self->is_temp ? 'TEMPORARY' : '';
  foreach (keys %$tables) {
    next if $_ eq 'meta';      # don't get rid of meta data!
    my $table = $self->_qualify($_);
    $dbh->do("DROP TABLE IF EXISTS $table") if $erase;
    my @table_exists = $dbh->selectrow_array("SELECT * FROM pg_tables WHERE tablename = '$table' AND schemaname = '$self->namespace'");
		if (!scalar(@table_exists)) {
			my $query = "CREATE $temporary TABLE $table $tables->{$_}";
			$dbh->do($query) or $self->throw($dbh->errstr);
		}
  }
  $self->subfeatures_are_indexed(1) if $erase;
  1;
}

sub maybe_create_meta {
  my $self = shift;
  return unless $self->writeable;
  my $namespace    = $self->namespace;
  my $table        = $self->_qualify('meta');
  my $tables       = $self->table_definitions;
  my $temporary    = $self->is_temp ? 'TEMPORARY' : '';
  my @table_exists = $self->dbh->selectrow_array("SELECT * FROM pg_tables WHERE tablename = 'meta' AND schemaname = '$namespace'");
  $self->dbh->do("CREATE $temporary TABLE $table $tables->{meta}")
      unless @table_exists;
}

###
# check if the namespace/schema exists, if not create it
#

sub _check_for_namespace {
  my $self = shift;
  my $namespace = $self->namespace;
  return if $namespace eq 'public';
  my $dbh  = $self->dbh;
  my @schema_exists = $dbh->selectrow_array("SELECT * FROM pg_namespace WHERE nspname = '$namespace'");
  if (!scalar(@schema_exists)) {
      my $query = "CREATE SCHEMA $namespace";
	  $dbh->do($query) or $self->throw($dbh->errstr);

	  # if temp parameter is set and schema created for this process then enable removal in remove_namespace()
	  if ($self->is_temp) {
	      $self->{delete_schema} = 1;
	  }
  }
}

###
# Overiding inherited mysql _qualify (We do not need to qualify for PostgreSQL as we have set the search_path above)
#
sub _qualify {
  my $self = shift;
  my $table_name = shift;
  return $table_name;
}

###
# when is_temp is set and the schema did not exist beforehand then we are able to remove it
#
sub remove_namespace {
  my $self = shift;
  if ($self->{delete_schema}) {
      my $namespace = $self->namespace;
      $self->dbh->do("DROP SCHEMA $namespace") or $self->throw($self->dbh->errstr);
  }
}

####Overiding the inherited mysql function _prepare

sub _prepare {
   my $self = shift;
   my $query = shift;
   my $dbh   = $self->dbh;
   my $schema = $self->{namespace};

   if ($schema) {
     $self->_check_for_namespace();
     $dbh->do("SET search_path TO " . $self->{'schema'} );
   } else {
     $dbh->do("SET search_path TO public");
   }
   my $sth   = $dbh->prepare_cached($query, {}, 3) or
   $self->throw($dbh->errstr);
   $sth;
}

sub _finish_bulk_update {
  my $self = shift;
  my $dbh  = $self->dbh;
  my $dir  = $self->{dumpdir} || '.';
  for my $table ('feature',$self->index_tables) {
    my $fh   = $self->dump_filehandle($table);
    my $path = $self->dump_path($table);
    $fh->close;
    my $qualified_table = $self->_qualify($table);
    copy($path, "$path.bak");
    # Get stuff from file into STDIN so we don't have to be superuser
    open my $FH, '<', $path or $self->throw("Could not read file '$path': $!");
    print STDERR "Loading file $path\n";
    $dbh->do("COPY $qualified_table FROM STDIN CSV QUOTE '''' DELIMITER '\t'") or $self->throw($dbh->errstr);
    while (my $line = <$FH>) {
      $dbh->pg_putline($line);
    }
    $dbh->pg_endcopy() or $self->throw($dbh->errstr);
    close $FH;
    #unlink $path;
  }
  delete $self->{bulk_update_in_progress};
  delete $self->{filehandles};
}


###
# Add a subparts to a feature. Both feature and all subparts must already be in database.
#
sub _add_SeqFeature {
  my $self     = shift;

  # special purpose method for case when we are doing a bulk update
  return $self->_dump_add_SeqFeature(@_) if $self->{bulk_update_in_progress};

  my $parent   = shift;
  my @children = @_;

  my $dbh = $self->dbh;
  local $dbh->{RaiseError} = 1;

  my $child_table = $self->_parent2child_table();
  my $count = 0;

  my $querydel = "DELETE FROM $child_table WHERE id = ? AND child = ?";
  my $query = "INSERT INTO $child_table (id,child) VALUES (?,?)";
  my $sthdel = $self->_prepare($querydel);
  my $sth = $self->_prepare($query);

  my $parent_id = (ref $parent ? $parent->primary_id : $parent)
    or $self->throw("$parent should have a primary_id");

  $self->begin_work or $self->throw($dbh->errstr);
  eval {
    for my $child (@children) {
      my $child_id = ref $child ? $child->primary_id : $child;
      defined $child_id or die "no primary ID known for $child";
      $sthdel->execute($parent_id, $child_id);
      $sth->execute($parent_id,$child_id);
      $count++;
    }
  };

  if ($@) {
    warn "Transaction aborted because $@";
    $self->rollback;
  }
  else {
    $self->commit;
  }
  $sth->finish;
  $count;
}

# because this is a reserved word in postgresql
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

  my $sth     = $self->_prepare(<<END);
SELECT sequence,"offset"
   FROM $sequence_table as s,$locationlist_table as ll
   WHERE s.id=ll.id
     AND ll.seqname= ?
     AND "offset" >= ?
     AND "offset" <= ?
   ORDER BY "offset"
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
  $sql =  $position eq 'left'  ? "SELECT min(\"offset\") FROM $sequence_table as s,$locationlist_table as ll WHERE s.id=ll.id AND ll.seqname=?"
         :$position eq 'right' ? "SELECT max(\"offset\") FROM $sequence_table as s,$locationlist_table as ll WHERE s.id=ll.id AND ll.seqname=?"
	 :"SELECT max(\"offset\") FROM $sequence_table as s,$locationlist_table as ll WHERE s.id=ll.id AND ll.seqname=? AND \"offset\"<=?";

  my $sth = $self->_prepare($sql);
  my @args = $position =~ /^-?\d+$/ ? ($seqid,$position) : ($seqid);
  $sth->execute(@args) or $self->throw($sth->errstr);
  my $boundary = $sth->fetchall_arrayref->[0][0];
  $sth->finish;
  return $boundary;
}

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
  my $sql_regexp = join ' OR ',("a.attribute_value SIMILAR TO ?")  x @words;
  my $sql = <<END;
SELECT name,attribute_value,tl.tag,n.id
  FROM $name_table as n,$attribute_table as a,$attributelist_table as al,$type_table as t,$typelist_table as tl
  WHERE n.id=a.id
    AND al.id=a.attribute_id
    AND n.id=t.id
    AND t.typeid=tl.id
    AND n.display_name=1
    AND ($tag_sql)
    AND ($sql_regexp)
END
  $sql .= "LIMIT $limit" if defined $limit;
  $self->_print_query($sql,@tags,@wild_card_words) if DEBUG || $self->debug;
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

# overridden here because the mysql adapter uses
# a non-standard query hint
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

sub _match_sql {
  my $self = shift;
  my $name = shift;

  my ($match,$string);
  if ($name =~ /(?:^|[^\\])[*?]/) {
    $name =~ s/(^|[^\\])([%_])/$1\\$2/g;
    $name =~ s/(^|[^\\])\*/$1%/g;
    $name =~ s/(^|[^\\])\?/$1_/g;
    $match = "LIKE ?";
    $string  = lc($name);
  } else {
    $match = "= lower(?)";
    $string  = lc($name);
  }
  return ($match,$string);
}

# overridden because of differences between LIKE behavior in mysql and postgres
# as well as case-sensitivity of matches
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

    if ($source_tag) {
      push @matches,"lower(tl.tag)=lower(?)";
      push @args,"$primary_tag:$source_tag";
    } else {
      push @matches,"tl.tag ILIKE ?";
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

# overridden because mysql adapter uses the non-standard REPLACE syntax
sub setting {
  my $self = shift;
  my ($variable_name,$value) = @_;
  my $meta  = $self->_meta_table;

  if (defined $value && $self->writeable) {
    my $querydel = "DELETE FROM $meta WHERE name = ?";
    my $query    = "INSERT INTO $meta (name,value) VALUES (?,?)";
    my $sthdel   = $self->_prepare($querydel);
    my $sth      = $self->_prepare($query);
    $sthdel->execute($variable_name);
    $sth->execute($variable_name,$value) or $self->throw($sth->errstr);
    $sth->finish;
    $self->{settings_cache}{$variable_name} = $value;
  }
  else {
    return $self->{settings_cache}{$variable_name} if exists $self->{settings_cache}{$variable_name};
    my $query = "SELECT value FROM $meta as m WHERE m.name=?";
    my $sth = $self->_prepare($query);
#    $sth->execute($variable_name) or $self->throw($sth->errstr);
    unless ($sth->execute($variable_name)) {
      my $errstr = $sth->errstr;
      $sth = $self->_prepare("SHOW search_path");
      $sth->execute();
      $errstr .= "With search_path " . $sth->fetchrow_arrayref->[0] . "\n";
      $self->throw($errstr);
    }

    my ($value) = $sth->fetchrow_array;
    $sth->finish;
    return $self->{settings_cache}{$variable_name} = $value;
  }
}

# overridden because of use of REPLACE in mysql adapter
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

  my $query = "INSERT INTO $features (id,object,indexed,seqid,start,\"end\",strand,tier,bin,typeid) VALUES (?,?,?,?,?,?,?,?,?,?)";
  my $query_noid = "INSERT INTO $features (object,indexed,seqid,start,\"end\",strand,tier,bin,typeid) VALUES (?,?,?,?,?,?,?,?,?)";
  my $querydel = "DELETE FROM $features WHERE id = ?";

  my $sthdel = $self->_prepare($querydel);
  my $sth = $self->_prepare($query);
  my $sth_noid = $self->_prepare($query_noid);

  my @location = $index_flag ? $self->_get_location_and_bin($object) : (undef)x6;

  my $primary_tag = $object->primary_tag;
  my $source_tag  = $object->source_tag || '';
  $primary_tag    .= ":$source_tag";
  my $typeid   = $self->_typeid($primary_tag,1);

  if ($id) {
    $sthdel->execute($id);
    $sth->execute($id,encode_base64($self->freeze($object), ''),$index_flag||0,@location,$typeid) or $self->throw($sth->errstr);
  } else {
    $sth_noid->execute(encode_base64($self->freeze($object), ''),$index_flag||0,@location,$typeid) or $self->throw($sth->errstr);
  }

  my $dbh = $self->dbh;

  $object->primary_id($dbh->last_insert_id(undef, undef, undef, undef, {sequence=>$features."_id_seq"})) unless defined $id;

  $self->flag_for_indexing($dbh->last_insert_id(undef, undef, undef, undef, {sequence=>$features."_id_seq"})) if $self->{bulk_update_in_progress};
}

=head2 types

 Title   : types
 Usage   : @type_list = $db->types
 Function: Get all the types in the database
 Returns : array of Bio::DB::GFF::Typename objects
 Args    : none
 Status  : public

=cut

# overridden because "offset" is reserved in postgres
###
# Insert a bit of DNA or protein into the database
#
sub _insert_sequence {
  my $self = shift;
  my ($seqid,$seq,$offset) = @_;
  my $id = $self->_locationid($seqid);
  my $seqtable = $self->_sequence_table;
  my $sthdel = $self->_prepare("DELETE FROM $seqtable WHERE id = ? AND \"offset\" = ?");
  my $sth = $self->_prepare(<<END);
INSERT INTO $seqtable (id,"offset",sequence) VALUES (?,?,?)
END
  $sthdel->execute($id,$offset);
  $sth->execute($id,$offset,$seq) or $self->throw($sth->errstr);
}

# overridden because of mysql adapter's use of REPLACE
###
# This subroutine flags the given primary ID for later reindexing
#
sub flag_for_indexing {
  my $self = shift;
  my $id   = shift;
  my $needs_updating = $self->_update_table;

  my $querydel = "DELETE FROM $needs_updating WHERE id = ?";
  my $query = "INSERT INTO $needs_updating VALUES (?)";
  my $sthdel = $self->_prepare($querydel);
  my $sth = $self->_prepare($query);

  $sthdel->execute($id);
  $sth->execute($id) or $self->throw($self->dbh->errstr);
}

# overridden because of the different ways that mysql and postgres
# handle id sequences
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
  return $dbh->last_insert_id(undef, undef, undef, undef, {sequence=>$qualified_table."_id_seq"});
}

# overridden because of differences in binding between mysql and postgres adapters
# given a statement handler that is expected to return rows of (id,object)
# unthaw each object and return a list of 'em
sub _sth2objs {
  my $self = shift;
  my $sth  = shift;
  my @result;
  my ($id, $o);
	$sth->bind_col(1, \$id);
	$sth->bind_col(2, \$o, { pg_type => PG_BYTEA});
  #while (my ($id,$o) = $sth->fetchrow_array) {
  while ($sth->fetch) {
    my $obj = $self->thaw(decode_base64($o) ,$id);
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
  my $obj = $self->thaw(decode_base64($o) ,$id);
  $obj;
}

####################################################################################################
# SQL Fragment generators
####################################################################################################

# overridden because of base64 encoding needed here
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

		my $frozen_object = encode_base64($self->freeze($obj), '');
                # TODO: Fix this, why does frozen object start with quote but not end with one
    print $store_fh join("\t",$id,$typeid,$seqid,$start,$end,$strand,$tier,$bin,$indexed,$frozen_object),"\n";
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

sub _enable_keys  { }  # nullop
sub _disable_keys { }  # nullop

sub _add_interval_stats_table {
    my $self = shift;
    my $tables          = $self->table_definitions;
    my $interval_stats  = $self->_interval_stats_table;
    ##check to see if it exists yet; if it does, just return because
    ##there is a drop from in the next step
    my $dbh = $self->dbh;
    my @table_exists = $dbh->selectrow_array("SELECT * FROM pg_tables WHERE tablename
 = '$interval_stats' AND schemaname = '".$self->namespace."'");
    if (!scalar(@table_exists)) {
        my $query = "CREATE TABLE $interval_stats $tables->{interval_stats}";
        $dbh->do($query) or $self->throw($dbh->errstr);
    }
}

sub _fetch_indexed_features_sql {
    my $self     = shift;
    my $features = $self->_feature_table;
    return <<END;
SELECT typeid,seqid,start-1,"end"
  FROM $features as f
 WHERE f.indexed=1
  ORDER BY typeid,seqid,start
END
}



1;
