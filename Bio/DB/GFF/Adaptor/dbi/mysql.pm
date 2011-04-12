package Bio::DB::GFF::Adaptor::dbi::mysql;

=head1 NAME

Bio::DB::GFF::Adaptor::dbi::mysql -- Database adaptor for a specific mysql schema

=head1 SYNOPSIS

See L<Bio::DB::GFF>

=cut

# a simple mysql adaptor
use strict;
use Bio::DB::GFF::Util::Rearrange; # for rearrange()
use Bio::DB::GFF::Util::Binning;
use base qw(Bio::DB::GFF::Adaptor::dbi);

use constant MAX_SEGMENT => 100_000_000;  # the largest a segment can get

use constant GETSEQCOORDS =><<END;
SELECT fref,
       IF(ISNULL(gclass),'Sequence',gclass),
       min(fstart),
       max(fstop),
       fstrand,
       gname
  FROM fdata,fgroup
  WHERE fgroup.gname=?
    AND fgroup.gclass=?
    AND fgroup.gid=fdata.gid
    GROUP BY fref,fstrand,gname
END
;

use constant GETALIASCOORDS =><<END;
SELECT fref,
       IF(ISNULL(gclass),'Sequence',gclass),
       min(fstart),
       max(fstop),
       fstrand,
       gname
  FROM fdata,fgroup,fattribute,fattribute_to_feature
  WHERE fattribute_to_feature.fattribute_value=?
    AND fgroup.gclass=?
    AND fgroup.gid=fdata.gid
    AND fattribute.fattribute_name='Alias'
    AND fattribute_to_feature.fattribute_id=fattribute.fattribute_id
    AND fattribute_to_feature.fid=fdata.fid
    GROUP BY fref,fstrand,gname
END
;

use constant GETALIASLIKE =><<END;
SELECT fref,
       IF(ISNULL(gclass),'Sequence',gclass),
       min(fstart),
       max(fstop),
       fstrand,
       gname
  FROM fdata,fgroup,fattribute,fattribute_to_feature
  WHERE fattribute_to_feature.fattribute_value LIKE ?
    AND fgroup.gclass=?
    AND fgroup.gid=fdata.gid
    AND fattribute.fattribute_name='Alias'
    AND fattribute_to_feature.fattribute_id=fattribute.fattribute_id
    AND fattribute_to_feature.fid=fdata.fid
    GROUP BY fref,fstrand,gname
END
;

use constant GETFORCEDSEQCOORDS =><<END;
SELECT fref,
       IF(ISNULL(gclass),'Sequence',gclass),
       min(fstart),
       max(fstop),
       fstrand
  FROM fdata,fgroup
  WHERE fgroup.gname=?
    AND fgroup.gclass=?
    AND fdata.fref=?
    AND fgroup.gid=fdata.gid
    GROUP BY fref,fstrand
END
;

use constant FULLTEXTSEARCH => <<END;
SELECT distinct gclass,gname,fattribute_value,MATCH(fattribute_value) AGAINST (?) as score,fmethod,fsource
  FROM fgroup,fattribute_to_feature,fdata,ftype
  WHERE fgroup.gid=fdata.gid
    AND fdata.fid=fattribute_to_feature.fid
    AND fdata.ftypeid=ftype.ftypeid
    AND MATCH(fattribute_value) AGAINST (?)
END
;

=head1 DESCRIPTION

This adaptor implements a specific mysql database schema that is
compatible with Bio::DB::GFF.  It inherits from
Bio::DB::GFF::Adaptor::dbi, which itself inherits from Bio::DB::GFF.

The schema uses several tables:

=over 4

=item fdata

This is the feature data table.  Its columns are:
-
    fid	           feature ID (integer)
    fref           reference sequence name (string)
    fstart         start position relative to reference (integer)
    fstop          stop postion relative to reference (integer)
    ftypeid        feature type ID (integer)
    fscore         feature score (float); may be null
    fstrand        strand; one of "+" or "-"; may be null
    fphase         phase; one of 0, 1 or 2; may be null
    gid            group ID (integer)
    ftarget_start  for similarity features, the target start position (integer)
    ftarget_stop   for similarity features, the target stop position (integer)

Note that it would be desirable to normalize the reference sequence
name, since there are usually many features that share the same
reference feature.  However, in the current schema, query performance
suffers dramatically when this additional join is added.

=item fgroup

This is the group table. There is one row for each group.  Columns:

    gid	      the group ID (integer)
    gclass    the class of the group (string)
    gname     the name of the group (string)

The group table serves multiple purposes.  As you might expect, it is
used to cluster features that logically belong together, such as the
multiple exons of the same transcript.  It is also used to assign a
name and class to a singleton feature.  Finally, the group table is
used to identify the target of a similarity hit.  This is consistent
with the way in which the group field is used in the GFF version 2
format.

The fgroup.gid field joins with the fdata.gid field. 

Examples:

  mysql> select * from fgroup where gname='sjj_2L52.1';
  +-------+-------------+------------+
  | gid   | gclass      | gname      |
  +-------+-------------+------------+
  | 69736 | PCR_product | sjj_2L52.1 |
  +-------+-------------+------------+
  1 row in set (0.70 sec)

  mysql> select fref,fstart,fstop from fdata,fgroup 
            where gclass='PCR_product' and gname = 'sjj_2L52.1' 
                  and fdata.gid=fgroup.gid;
  +---------------+--------+-------+
  | fref          | fstart | fstop |
  +---------------+--------+-------+
  | CHROMOSOME_II |   1586 |  2355 |
  +---------------+--------+-------+
  1 row in set (0.03 sec)

=item ftype

This table contains the feature types, one per row.  Columns are:

    ftypeid      the feature type ID (integer)
    fmethod      the feature type method name (string)
    fsource      the feature type source name (string)

The ftype.ftypeid field joins with the fdata.ftypeid field.  Example:

  mysql> select fref,fstart,fstop,fmethod,fsource from fdata,fgroup,ftype 
         where gclass='PCR_product' 
               and gname = 'sjj_2L52.1'
               and fdata.gid=fgroup.gid
               and fdata.ftypeid=ftype.ftypeid;
  +---------------+--------+-------+-------------+-----------+
  | fref          | fstart | fstop | fmethod     | fsource   |
  +---------------+--------+-------+-------------+-----------+
  | CHROMOSOME_II |   1586 |  2355 | PCR_product | GenePairs |
  +---------------+--------+-------+-------------+-----------+
  1 row in set (0.08 sec)

=item fdna

This table holds the raw DNA of the reference sequences.  It has three
columns:

    fref          reference sequence name (string)
    foffset       offset of this sequence
    fdna          the DNA sequence (longblob)

To overcome problems loading large blobs, DNA is automatically
fragmented into multiple segments when loading, and the position of
each segment is stored in foffset.  The fragment size is controlled by
the -clump_size argument during initialization.

=item fattribute_to_feature

This table holds "attributes", which are tag/value pairs stuffed into
the GFF line.  The first tag/value pair is treated as the group, and
anything else is treated as an attribute (weird, huh?).

 CHR_I assembly_tag Finished     2032 2036 . + . Note "Right: cTel33B"
 CHR_I assembly_tag Polymorphism 668  668  . + . Note "A->C in cTel33B"

The columns of this table are:

    fid                 feature ID (integer)
    fattribute_id       ID of the attribute (integer)
    fattribute_value    text of the attribute (text)

The fdata.fid column joins with fattribute_to_feature.fid.

=item fattribute

This table holds the normalized names of the attributes.  Fields are:

  fattribute_id      ID of the attribute (integer)
  fattribute_name    Name of the attribute (varchar)

=back

=head2 Data Loading Methods

In addition to implementing the abstract SQL-generating methods of
Bio::DB::GFF::Adaptor::dbi, this module also implements the data
loading functionality of Bio::DB::GFF.

=cut


=head2 new

 Title   : new
 Usage   : $db = Bio::DB::GFF->new(@args)
 Function: create a new adaptor
 Returns : a Bio::DB::GFF object
 Args    : see below
 Status  : Public

The new constructor is identical to the "dbi" adaptor's new() method,
except that the prefix "dbi:mysql" is added to the database DSN identifier
automatically if it is not there already.

  Argument       Description
  --------       -----------

  -dsn           the DBI data source, e.g. 'dbi:mysql:ens0040' or "ens0040"

  -user          username for authentication

  -pass          the password for authentication

=cut

#'

sub new {
  my $class = shift;
  my ($dsn,$other) = rearrange([
				[qw(FEATUREDB DB DSN)],
			       ],@_);
  $dsn = "dbi:mysql:$dsn" if !ref($dsn) && $dsn !~ /^(?:dbi|DBI):/;
  my $self = $class->SUPER::new(-dsn=>$dsn,%$other);
  $self;
}

=head2 get_dna

 Title   : get_dna
 Usage   : $string = $db->get_dna($name,$start,$stop,$class)
 Function: get DNA string
 Returns : a string
 Args    : name, class, start and stop of desired segment
 Status  : Public

This method performs the low-level fetch of a DNA substring given its
name, class and the desired range.  This should probably be moved to
the parent class.

=cut

sub getseqcoords_query {
   my $self = shift;
   return GETSEQCOORDS ;
}

sub getaliascoords_query{
  my $self = shift;
  return GETALIASCOORDS ;
}


sub getforcedseqcoords_query{
  my $self = shift;
  return GETFORCEDSEQCOORDS ;
}


sub getaliaslike_query{
  my $self = shift;
  return GETALIASLIKE ;
}


# override parent
sub get_abscoords_bkup {
  my $self = shift;
  my ($name,$class,$refseq)  = @_;

  my $result = $self->SUPER::get_abscoords(@_);
  return $result if $result;

  my $sth;
  if ($name =~ s/\*/%/g) {
    $sth = $self->dbh->do_query(GETALIASLIKE,$name,$class);
  } else {
    $sth = $self->dbh->do_query(GETALIASCOORDS,$name,$class);
  }
  my @result;
  while (my @row = $sth->fetchrow_array) { push @result,\@row }
  $sth->finish;

  if (@result == 0) {
    $self->error("$name not found in database");
    return;
  } else {
    return \@result;
  }

}



sub make_features_select_part {
  my $self = shift;
  my $options = shift || {};
  my $s;
  if (my $b = $options->{bin_width}) {

    $s = <<END;
fref,
  1+$b*floor(fstart/$b)   as fstart,
  $b*(1+floor(fstart/$b)) as fstop,
  IF(ISNULL(fsource),fmethod,concat(fmethod,':',fsource)),'bin',
  count(*) as fscore,
  '.','.','bin',
  IF(ISNULL(fsource),concat(fref,':',fmethod),concat(fref,':',fmethod,':',fsource)),
  NULL,NULL,NULL,NULL
END
;
  } else {
    $s = <<END;
fref,fstart,fstop,fsource,fmethod,fscore,fstrand,fphase,gclass,gname,ftarget_start,ftarget_stop,fdata.fid,fdata.gid
END
;
}
  $s .= ",count(fdata.fid)" if $options->{attributes} && keys %{$options->{attributes}}>1;
  $s;
}


# IMPORTANT NOTE:
# WHETHER OR NOT THIS WORKS IS CRITICALLY DEPENDENT ON THE RELATIVE MAGNITUDE OF THE
sub make_features_from_part {
  my $self = shift;
  my $sparse_types  = shift;
  my $options       = shift || {};
  my $sparse_groups = $options->{sparse_groups};
  my $index =  $sparse_groups ? ' USE INDEX(gid)'
             : $sparse_types  ? ' USE INDEX(ftypeid)'
             : '';
  return $options->{attributes} ? "fdata${index},ftype,fgroup,fattribute,fattribute_to_feature\n"
                                : "fdata${index},ftype,fgroup\n";
}

=head2 search_notes

 Title   : search_notes
 Usage   : @search_results = $db->search_notes("full text search string",$limit)
 Function: Search the notes for a text string, using mysql full-text search
 Returns : array of results
 Args    : full text search string, and an optional row limit
 Status  : public

This is a mysql-specific method.  Given a search string, it performs a
full-text search of the notes table and returns an array of results.
Each row of the returned array is a arrayref containing the following fields:

  column 1     A Bio::DB::GFF::Featname object, suitable for passing to segment()
  column 2     The text of the note
  column 3     A relevance score.

=cut

sub search_notes {
  my $self = shift;
  my ($search_string,$limit) = @_;

  $search_string =~ tr/*?//d;

  my $query = FULLTEXTSEARCH;
  $query .= " limit $limit" if defined $limit;
  my $sth = $self->dbh->do_query($query,$search_string,$search_string);
  my @results;
  while (my ($class,$name,$note,$relevance,$method,$source) = $sth->fetchrow_array) {
     next unless $class && $name;    # sorry, ignore NULL objects
     $relevance = sprintf("%.2f",$relevance);  # trim long floats
     my $featname = Bio::DB::GFF::Featname->new($class=>$name);
     my $type     = Bio::DB::GFF::Typename->new($method,$source);
     push @results,[$featname,$note,$relevance,$type];
  }

  #added result filtering so that this method returns the expected results
  #this section of code used to be in GBrowse's do_keyword_search method

  my $match_sub = 'sub {';
  foreach (split /\s+/,$search_string) {
    $match_sub .= "return unless \$_[0] =~ /\Q$_\E/i; ";
  }
  $match_sub .= "};";
  my $match = eval $match_sub;

  my @matches = grep { $match->($_->[1]) } @results;

  return @matches;
}



################################ loading and initialization ##################################

=head2 schema

 Title   : schema
 Usage   : $schema = $db->schema
 Function: return the CREATE script for the schema
 Returns : a list of CREATE statemetns
 Args    : none
 Status  : protected

This method returns a list containing the various CREATE statements
needed to initialize the database tables.

=cut

sub schema {
    my $self = shift;
    my $dbh  = $self->dbh;
    my ($version) = $dbh->selectrow_array('select version()');
    my ($major, $minor) = split /\./, $version;
    $version = "$major.$minor";
    my $engine    = $version >= 4.1 ? 'ENGINE' : 'TYPE';
    my %schema = (
		fdata =>{ 
table=> qq{
 create table fdata (
    fid	                int not null  auto_increment,
    fref                varchar(100) not null,
    fstart              int not null,
    fstop               int not null,
    fbin                double precision,
    ftypeid             int not null,
    fscore              float,
    fstrand             enum('+','-'),
    fphase              enum('0','1','2'),
    gid                 int not null,
    ftarget_start       int,
    ftarget_stop        int,
    primary key(fid),
    unique index(fref,fbin,fstart,fstop,ftypeid,gid),
    index(ftypeid),
    index(gid)
		   ) $engine=MyISAM
}  # fdata table
}, # fdata

		fgroup =>{ 
table=> qq{
create table fgroup (
    gid	    int not null  auto_increment,
    gclass  varchar(100),
    gname   varchar(100),
    primary key(gid),
    unique(gclass,gname)
)  $engine=MyISAM
}
},

          ftype => {
table=> qq{
create table ftype (
    ftypeid      int not null   auto_increment,
    fmethod       varchar(100) not null,
    fsource       varchar(100),
    primary key(ftypeid),
    index(fmethod),
    index(fsource),
    unique ftype (fmethod,fsource)
)  $engine=MyISAM
}  #ftype table
}, #ftype

         fdna => {
table=> qq{
create table fdna (
		fref    varchar(100) not null,
	        foffset int(10) unsigned not null,
	        fdna    longblob,
		primary key(fref,foffset)
)   $engine=MyISAM
} # fdna table
},#fdna

        fmeta => {
table=> qq{
create table fmeta (
		fname   varchar(255) not null,
	        fvalue  varchar(255) not null,
		primary key(fname)
)  $engine=MyISAM
} # fmeta table
},#fmeta

       fattribute => {
table=> qq{
create table fattribute (
	fattribute_id     int(10)         unsigned not null auto_increment,
        fattribute_name   varchar(255)    not null,
	primary key(fattribute_id)
)  $engine=MyISAM
} #fattribute table
},#fattribute

       fattribute_to_feature => {
table=> qq{
create table fattribute_to_feature (
        fid              int(10) not null,
        fattribute_id    int(10) not null,
	fattribute_value text,
        key(fid,fattribute_id),
	key(fattribute_value(48)),
        fulltext(fattribute_value)
)  $engine=MyISAM
} # fattribute_to_feature table
},# fattribute_to_feature

       finterval_stats => {
table=> qq{
create table finterval_stats (
   ftypeid            integer not null,
   fref               varchar(100) not null,
   fbin               integer not null,
   fcum_count         integer not null,
   primary key(ftypeid,fref,fbin)
)  $engine=MyISAM
} # finterval_stats table
},# finterval_stats

);
  return \%schema;
}



=head2 make_classes_query

 Title   : make_classes_query
 Usage   : ($query,@args) = $db->make_classes_query
 Function: return query fragment for generating list of reference classes
 Returns : a query and args
 Args    : none
 Status  : public

=cut

sub make_classes_query {
  my $self = shift;
  return 'SELECT DISTINCT gclass FROM fgroup WHERE NOT ISNULL(gclass)';
}


=head2 make_meta_set_query

 Title   : make_meta_set_query
 Usage   : $sql = $db->make_meta_set_query
 Function: return SQL fragment for setting a meta parameter
 Returns : SQL fragment
 Args    : none
 Status  : public

By default this does nothing; meta parameters are not stored or
retrieved.

=cut

sub make_meta_set_query {
   return 'REPLACE INTO fmeta VALUES (?,?)';
}

=head2 setup_load

 Title   : setup_load
 Usage   : $db->setup_load
 Function: called before load_gff_line()
 Returns : void
 Args    : none
 Status  : protected

This method performs schema-specific initialization prior to loading a
set of GFF records.  It prepares a set of DBI statement handlers to be 
used in loading the data.

=cut

sub setup_load {
  my $self      = shift;

  my $dbh = $self->features_db;

  if ($self->lock_on_load) {
    my @tables = map { "$_ WRITE"} $self->tables;
    my $tables = join ', ',@tables;
    $dbh->do("LOCK TABLES $tables");
  }
#  for my $table (qw(fdata)) {
#    $dbh->do("alter table $table disable keys");
#  }

  my $lookup_type = $dbh->prepare_delayed('SELECT ftypeid FROM ftype WHERE fmethod=? AND fsource=?');
  my $insert_type = $dbh->prepare_delayed('INSERT INTO ftype (fmethod,fsource) VALUES (?,?)');

  my $lookup_group = $dbh->prepare_delayed('SELECT gid FROM fgroup WHERE gname=? AND gclass=?');
  my $insert_group = $dbh->prepare_delayed('INSERT INTO fgroup (gname,gclass) VALUES (?,?)');

  my $lookup_attribute = $dbh->prepare_delayed('SELECT fattribute_id FROM fattribute WHERE fattribute_name=?');
  my $insert_attribute = $dbh->prepare_delayed('INSERT INTO fattribute (fattribute_name) VALUES (?)');
  my $insert_attribute_value = $dbh->prepare_delayed('INSERT INTO fattribute_to_feature (fid,fattribute_id,fattribute_value) VALUES (?,?,?)');

  my $insert_data  = $dbh->prepare_delayed(<<END);
INSERT INTO fdata (fref,fstart,fstop,fbin,ftypeid,fscore,
		   fstrand,fphase,gid,ftarget_start,ftarget_stop)
       VALUES(?,?,?,?,?,?,?,?,?,?,?)
END
;


  $self->{load_stuff}{sth}{lookup_ftype}     = $lookup_type;
  $self->{load_stuff}{sth}{insert_ftype}     = $insert_type;
  $self->{load_stuff}{sth}{lookup_fgroup}    = $lookup_group;
  $self->{load_stuff}{sth}{insert_fgroup}    = $insert_group;
  $self->{load_stuff}{sth}{insert_fdata}     = $insert_data;
  $self->{load_stuff}{sth}{lookup_fattribute} = $lookup_attribute;
  $self->{load_stuff}{sth}{insert_fattribute} = $insert_attribute;
  $self->{load_stuff}{sth}{insert_fattribute_value} = $insert_attribute_value;
  $self->{load_stuff}{types}  = {};
  $self->{load_stuff}{groups} = {};
  $self->{load_stuff}{counter} = 0;
}

=head2 load_gff_line

 Title   : load_gff_line
 Usage   : $db->load_gff_line($fields)
 Function: called to load one parsed line of GFF
 Returns : true if successfully inserted
 Args    : hashref containing GFF fields
 Status  : protected

This method is called once per line of the GFF and passed a series of
parsed data items that are stored into the hashref $fields.  The keys are:

 ref          reference sequence
 source       annotation source
 method       annotation method
 start        annotation start
 stop         annotation stop
 score        annotation score (may be undef)
 strand       annotation strand (may be undef)
 phase        annotation phase (may be undef)
 group_class  class of annotation's group (may be undef)
 group_name   ID of annotation's group (may be undef)
 target_start start of target of a similarity hit
 target_stop  stop of target of a similarity hit
 attributes   array reference of attributes, each of which is a [tag=>value] array ref

=cut

sub load_gff_line {
  my $self = shift;
  my $gff = shift;

  my $s    = $self->{load_stuff};
  my $dbh  = $self->features_db;
  local $dbh->{PrintError} = 0;

  defined(my $typeid  = $self->get_table_id('ftype', $gff->{method} => $gff->{source})) or return;
  defined(my $groupid = $self->get_table_id('fgroup',$gff->{gname}  => $gff->{gclass})) or return;

  if ($gff->{stop}-$gff->{start}+1 > $self->max_bin) {
    warn "$gff->{gclass}:$gff->{gname} is ",$gff->{stop}-$gff->{start}+1,
      " bp long, but the maximum indexable feature is set to ",$self->max_bin," bp.\n";
    warn "Please set the maxbin value to a length at least as large as the largest feature you wish to store.\n";
    warn "\n* You will need to reinitialize the database from scratch.\n";
    warn "* With the Perl API you do this using the -max_bin argument to \$db->initialize().\n";
    warn "* With the command-line tools you do with this with --maxfeature option.\n";
  }

  my $bin =  bin($gff->{start},$gff->{stop},$self->min_bin);
  my $result = $s->{sth}{insert_fdata}->execute($gff->{ref},
					       $gff->{start},$gff->{stop},$bin,
					       $typeid,
					       $gff->{score},$gff->{strand},$gff->{phase},
					       $groupid,
					       $gff->{tstart},$gff->{tstop});

  warn $dbh->errstr,"\n" && return unless $result;

  my $fid = $dbh->{mysql_insertid}
    || $self->get_feature_id($gff->{ref},$gff->{start},$gff->{stop},$typeid,$groupid);


  # insert attributes
  foreach (@{$gff->{attributes}}) {
    defined(my $attribute_id = $self->get_table_id('fattribute',$_->[0])) or return;
    $s->{sth}{insert_fattribute_value}->execute($fid,$attribute_id,$_->[1]);
  }

  if ( (++$s->{counter} % 1000) == 0) {
    print STDERR "$s->{counter} records loaded...";
    print STDERR -t STDOUT && !$ENV{EMACS} ? "\r" : "\n";
  }

  $fid;
}

sub finish_load {
  my $self = shift;
  my $dbh = $self->features_db;
  local $dbh->{PrintError} = 0;
#  for my $table (qw(fdata)) {
#    $dbh->do("alter table $table enable keys");
#  }
  $self->SUPER::finish_load;
}


sub insert_sequence {
  my $self = shift;
  my($id,$offset,$seq) = @_;
  my $sth = $self->{_insert_sequence}
    ||= $self->dbh->prepare_delayed('replace into fdna values (?,?,?)');
  $sth->execute($id,$offset,$seq) or $self->throw($sth->errstr);
}


=head2 get_table_id

 Title   : get_table_id
 Usage   : $integer = $db->get_table_id($table,@ids)
 Function: get the ID of a group or type
 Returns : an integer ID or undef
 Args    : none
 Status  : private

This internal method is called by load_gff_line to look up the integer
ID of an existing feature type or group.  The arguments are the name
of the table, and two string identifiers.  For feature types, the
identifiers are the method and source.  For groups, the identifiers
are group name and class.

This method requires that a statement handler named I<lookup_$table>,
have been created previously by setup_load().  It is here to overcome
deficiencies in mysql's INSERT syntax.

=cut

#'
# get the object ID from a named table
sub get_table_id {
  my $self   = shift;
  my $table  = shift;
  my @ids    = @_;

  # irritating warning for null id
  my $id_key;
  {
    local $^W=0;
    $id_key = join ':',@ids;
  }

  my $s   = $self->{load_stuff};
  my $sth = $s->{sth};
  my $dbh = $self->features_db;

  unless (defined($s->{$table}{$id_key})) {

    #########################################
    # retrieval of the last inserted id is now located at the adaptor and not in caching_handle
    #######################################
    if ( (my $result = $sth->{"lookup_$table"}->execute(@ids)) > 0) {
      $s->{$table}{$id_key} = ($sth->{"lookup_$table"}->fetchrow_array)[0];
    } else {
      $sth->{"insert_$table"}->execute(@ids)
	&& ($s->{$table}{$id_key} = $self->insertid($sth->{"insert_$table"}));
	#&& ($s->{$table}{$id_key} = $sth->{"insert_$table"}{sth}{mysql_insertid});
	#&& ($s->{$table}{$id_key} = $sth->{"insert_$table"}->insertid);
    }
  }

  my $id = $s->{$table}{$id_key};
  unless (defined $id) {
    warn "No $table id for $id_key ",$dbh->errstr," Record skipped.\n";
    return;
  }
  $id;
}

sub insertid {
  my $self = shift;
  my $s = shift ;
  $s->{mysql_insertid};
}


=head2 get_feature_id

 Title   : get_feature_id
 Usage   : $integer = $db->get_feature_id($ref,$start,$stop,$typeid,$groupid)
 Function: get the ID of a feature
 Returns : an integer ID or undef
 Args    : none
 Status  : private

This internal method is called by load_gff_line to look up the integer
ID of an existing feature.  It is ony needed when replacing a feature
with new information.

=cut

# this method is called when needed to look up a feature's ID
sub get_feature_id {
  my $self = shift;
  my ($ref,$start,$stop,$typeid,$groupid) = @_;
  my $s = $self->{load_stuff};
  unless ($s->{get_feature_id}) {
    my $dbh = $self->features_db;
    $s->{get_feature_id} =
      $dbh->prepare_delayed('SELECT fid FROM fdata WHERE fref=? AND fstart=? AND fstop=? AND ftypeid=? AND gid=?');
  }
  my $sth = $s->{get_feature_id} or return;
  $sth->execute($ref,$start,$stop,$typeid,$groupid) or return;
  my ($fid) = $sth->fetchrow_array;
  return $fid;
}

sub _add_interval_stats_table {
    my $self              = shift;
    my $schema            = $self->schema;
    my $create_table_stmt = $schema->{'finterval_stats'}{'table'};
    $create_table_stmt    =~ s/create table/create table if not exists/i;
    my $dbh               = $self->features_db;
    $dbh->do($create_table_stmt) ||  warn $dbh->errstr;
}

sub disable_keys {
    my $self = shift;
    my $table = shift;
    my $dbh   = $self->dbh;
    $dbh->do("alter table $table disable keys");
}
sub enable_keys {
    my $self = shift;
    my $table = shift;
    my $dbh   = $self->dbh;
    $dbh->do("alter table $table enable keys");
}


1;



__END__

=head1 BUGS

none ;-)

=head1 SEE ALSO

L<Bio::DB::GFF>, L<bioperl>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2002 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

