package Bio::DB::GFF::Adaptor::dbi::mysqlcmap;

=head1 NAME

Bio::DB::GFF::Adaptor::dbi::mysqlcmap -- Database adaptor for an integraded
CMap/GBrowse mysql schema

=head1 SYNOPSIS

See L<Bio::DB::GFF>

=cut

# a simple mysql adaptor
use strict;
use Data::Dumper;
use Bio::DB::GFF::Adaptor::dbi;
use Bio::DB::GFF::Util::Rearrange; # for rearrange()
use Bio::DB::GFF::Util::Binning;
use base qw(Bio::DB::GFF::Adaptor::dbi::mysql);
require Bio::DB::GFF::Adaptor::dbi::mysql;

use constant GETSEQCOORDS =><<END;
SELECT fref,
       IF(ISNULL(gclass),'Sequence',gclass),
       min(fstart),
       max(fstop),
       fstrand,
       feature_name as gname
  FROM fdata,cmap_feature
  WHERE cmap_feature.feature_name=?
    AND cmap_feature.gclass=?
    AND cmap_feature.feature_id=fdata.feature_id
    GROUP BY fref,fstrand,feature_name
END
;

use constant GETALIASCOORDS =><<END;
SELECT fref,
       IF(ISNULL(gclass),'Sequence',gclass),
       min(fstart),
       max(fstop),
       fstrand,
       feature_name as gname
  FROM fdata,cmap_feature,fattribute,fattribute_to_feature
  WHERE fattribute_to_feature.fattribute_value=?
    AND cmap_feature.gclass=?
    AND cmap_feature.feature_id=fdata.feature_id
    AND fattribute.fattribute_name='Alias'
    AND fattribute_to_feature.fattribute_id=fattribute.fattribute_id
    AND fattribute_to_feature.fid=fdata.fid
    GROUP BY fref,fstrand,feature_name
END
;

use constant GETALIASLIKE =><<END;
SELECT fref,
       IF(ISNULL(gclass),'Sequence',gclass),
       min(fstart),
       max(fstop),
       fstrand,
       feature_name as gname
  FROM fdata,cmap_feature,fattribute,fattribute_to_feature
  WHERE fattribute_to_feature.fattribute_value LIKE ?
    AND cmap_feature.gclass=?
    AND cmap_feature.feature_id=fdata.feature_id
    AND fattribute.fattribute_name='Alias'
    AND fattribute_to_feature.fattribute_id=fattribute.fattribute_id
    AND fattribute_to_feature.fid=fdata.fid
    GROUP BY fref,fstrand,feature_name
END
;

use constant GETFORCEDSEQCOORDS =><<END;
SELECT fref,
       IF(ISNULL(gclass),'Sequence',gclass),
       min(fstart),
       max(fstop),
       fstrand
  FROM fdata,cmap_feature
  WHERE cmap_feature.feature_name=?
    AND cmap_feature.gclass=?
    AND fdata.fref=?
    AND cmap_feature.feature_id=fdata.feature_id
    GROUP BY fref,fstrand
END
;

use constant FULLTEXTSEARCH => <<END;
SELECT distinct gclass,feature_name,fattribute_value,MATCH(fattribute_value) AGAINST (?) as score
  FROM cmap_feature,fattribute_to_feature,fdata
  WHERE cmap_feature.feature_id=fdata.feature_id
    AND fdata.fid=fattribute_to_feature.fid
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
    feature_id     group ID used to be 'gid' (integer)
    ftarget_start  for similarity features, the target start position (integer)
    ftarget_stop   for similarity features, the target stop position (integer)

Note that it would be desirable to normalize the reference sequence
name, since there are usually many features that share the same
reference feature.  However, in the current schema, query performance
suffers dramatically when this additional join is added.

=item cmap_feature (replaces fgroup)

This is the group table. There is one row for each group.  This is the 
shared table between CMap and GBrowse.  There are many CMap related 
columns but only a few that GBrowse uses.  

GBrowse Columns:

    feature_id     the group ID (integer)
    gclass         the class of the group (string)
    feature_name   the name of the group (string)

The group table serves multiple purposes.  As you might expect, it is
used to cluster features that logically belong together, such as the
multiple exons of the same transcript.  It is also used to assign a
name and class to a singleton feature.  Finally, the group table is
used to identify the target of a similarity hit.  This is consistent
with the way in which the group field is used in the GFF version 2
format.

The cmap_feature.feature_id field joins with the fdata.feature_id field. 

Examples:

  mysql> select * from cmap_feature where feature_name='sjj_2L52.1';
  +--------------+-------------+--------------+
  | feature_id   | gclass      | feature_name |
  +--------------+-------------+--------------+
  | 69736        | PCR_product | sjj_2L52.1   |
  +--------------+-------------+--------------+
  1 row in set (0.70 sec)

  mysql> select fref,fstart,fstop from fdata,cmap_feature 
            where gclass='PCR_product' and feature_name = 'sjj_2L52.1' 
                  and fdata.feature_id=cmap_feature.feature_id;
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

  mysql> select fref,fstart,fstop,fmethod,fsource from fdata,cmap_feature,ftype 
         where gclass='PCR_product' 
               and feature_name = 'sjj_2L52.1'
               and fdata.feature_id=cmap_feature.feature_id
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

#Defined in mysql.pm

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
fref,fstart,fstop,fsource,fmethod,fscore,fstrand,fphase,gclass,feature_name as gname,ftarget_start,ftarget_stop,fdata.fid,fdata.feature_id
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
  my $index =  $sparse_groups ? ' USE INDEX(feature_id)'
             : $sparse_types  ? ' USE INDEX(ftypeid)'
             : '';
  return $options->{attributes} ? "fdata${index},ftype,cmap_feature,fattribute,fattribute_to_feature\n"
                                : "fdata${index},ftype,cmap_feature\n";
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
  my %schema = (
		fdata =>{ 
table=> q{
#create table fdata (
#    fid	         int not null  auto_increment,
#    fref         varchar(100)    not null,
#    fstart       int unsigned   not null,
#    fstop        int unsigned   not null,
#    ftypeid      int not null,
#    fscore        float,
#    fstrand       enum('+','-'),
#    fphase        enum('0','1','2'),
#    feature_id          int not null,
#    ftarget_start int unsigned,
#    ftarget_stop  int unsigned,
#    primary key(fid),
#    unique index(fref,fstart,fstop,ftypeid,feature_id),
#    index(ftypeid),
#    index(feature_id)
#) type=MyISAM


 create table fdata (
    fid	                int not null  auto_increment,
    fref                varchar(100) not null,
    fstart              int unsigned   not null,
    fstop               int unsigned   not null,
    fbin                double(20,6)  not null,
    ftypeid             int not null,
    fscore              float,
    fstrand             enum('+','-'),
    fphase              enum('0','1','2'),
    feature_id                 int not null,
    ftarget_start       int unsigned,
    ftarget_stop        int unsigned,
    primary key(fid),
    unique index(fref,fbin,fstart,fstop,ftypeid,feature_id),
    index(ftypeid),
    index(feature_id)
		   ) type=MyISAM
}  # fdata table
}, # fdata

          ftype => {
table=> q{
create table ftype (
    ftypeid      int not null   auto_increment,
    fmethod       varchar(100) not null,
    fsource       varchar(100),
    primary key(ftypeid),
    index(fmethod),
    index(fsource),
    unique ftype (fmethod,fsource)
)type=MyISAM
}  #ftype table
}, #ftype

         fdna => {
table=> q{
create table fdna (
		fref    varchar(100) not null,
	        foffset int(10) unsigned not null,
	        fdna    longblob,
		primary key(fref,foffset)
)type=MyISAM
} # fdna table
},#fdna

        fmeta => {
table=> q{
create table fmeta (
		fname   varchar(255) not null,
	        fvalue  varchar(255) not null,
		primary key(fname)
)type=MyISAM
} # fmeta table
},#fmeta

       fattribute => {
table=> q{
create table fattribute (
	fattribute_id     int(10)         unsigned not null auto_increment,
        fattribute_name   varchar(255)    not null,
	primary key(fattribute_id)
)type=MyISAM
} #fattribute table
},#fattribute

       fattribute_to_feature => {
table=> q{
create table fattribute_to_feature (
        fid              int(10) not null,
        fattribute_id    int(10) not null,
	fattribute_value text,
        key(fid,fattribute_id),
	key(fattribute_value(48)),
        fulltext(fattribute_value)
)type=MyISAM
} # fattribute_to_feature table
    }, # fattribute_to_feature


cmap_attribute => {
table=>q{
create table cmap_attribute (
  attribute_id int(11) NOT NULL default '0',
  table_name varchar(30) NOT NULL default '',
  object_id int(11) NOT NULL default '0',
  display_order int(11) NOT NULL default '1',
  is_public tinyint(4) NOT NULL default '1',
  attribute_name varchar(200) NOT NULL default '',
  attribute_value text NOT NULL,
  PRIMARY KEY  (attribute_id),
  KEY table_name (table_name,object_id,display_order,attribute_name)
) TYPE=MyISAM;
} # table
},

cmap_correspondence_evidence => {
table=>q{
create table cmap_correspondence_evidence (
  correspondence_evidence_id int(11) NOT NULL default '0',
  accession_id varchar(20) NOT NULL default '',
  feature_correspondence_id int(11) NOT NULL default '0',
  evidence_type_accession varchar(20) NOT NULL default '0',
  score double(8,2) default NULL,
  rank int(11) NOT NULL default '0',
  PRIMARY KEY  (correspondence_evidence_id),
  UNIQUE KEY accession_id (accession_id),
  KEY feature_correspondence_id (feature_correspondence_id)
) TYPE=MyISAM;
} # table
},


cmap_correspondence_lookup => {
table=>q{
create table cmap_correspondence_lookup (
  feature_id1 int(11) default NULL,
  feature_id2 int(11) default NULL,
  feature_correspondence_id int(11) default NULL,
  start_position1 double(11,2) default NULL,
  start_position2 double(11,2) default NULL,
  stop_position1 double(11,2) default NULL,
  stop_position2 double(11,2) default NULL,
  map_id1 int(11) default NULL,
  map_id2 int(11) default NULL,
  feature_type_accession1 varchar(20) default NULL,
  feature_type_accession2 varchar(20) default NULL,
  KEY feature_id1 (feature_id1),
  KEY corr_id (feature_correspondence_id),
  KEY cl_map_id1 (map_id1),
  KEY cl_map_id2 (map_id2),
  KEY cl_map_id1_map_id2 (map_id1,map_id2),
  KEY cl_map_id2_map_id1 (map_id2,map_id1)
) TYPE=MyISAM;
} # table
},


cmap_correspondence_matrix => {
table=>q{
create table cmap_correspondence_matrix (
  reference_map_aid varchar(20) NOT NULL default '0',
  reference_map_name varchar(32) NOT NULL default '',
  reference_map_set_aid varchar(20) NOT NULL default '0',
  reference_species_aid varchar(20) NOT NULL default '0',
  link_map_aid varchar(20) default NULL,
  link_map_name varchar(32) default NULL,
  link_map_set_aid varchar(20) NOT NULL default '0',
  link_species_aid varchar(20) NOT NULL default '0',
  no_correspondences int(11) NOT NULL default '0'
) TYPE=MyISAM;
} # table
},


cmap_feature => {
table=>q{
create table cmap_feature (
  feature_id int(11) NOT NULL default '0',
  accession_id varchar(20) NOT NULL default '',
  map_id int(11) default NULL,
  feature_type_accession varchar(20) NOT NULL default '0',
  feature_name varchar(32) NOT NULL default '',
  is_landmark tinyint(4) NOT NULL default '0',
  start_position double(11,2) NOT NULL default '0.00',
  stop_position double(11,2) default NULL,
  default_rank int(11) NOT NULL default '1',
  direction tinyint(4) NOT NULL default '1',
  gclass varchar(100) default NULL,
  PRIMARY KEY  (feature_id),
  UNIQUE KEY gclass (gclass,feature_name),
  UNIQUE KEY accession_id (accession_id),
  KEY feature_name (feature_name),
  KEY feature_id_map_id (feature_id,map_id),
  KEY feature_id_map_id_start (feature_id,map_id,start_position),
  KEY map_id (map_id),
  KEY map_id_feature_id (map_id,feature_id)
) TYPE=MyISAM;
} # table
},


cmap_feature_alias => {
table=>q{
create table cmap_feature_alias (
  feature_alias_id int(11) NOT NULL default '0',
  feature_id int(11) NOT NULL default '0',
  alias varchar(255) default NULL,
  PRIMARY KEY  (feature_alias_id),
  UNIQUE KEY feature_id_2 (feature_id,alias),
  KEY feature_id (feature_id),
  KEY alias (alias)
) TYPE=MyISAM;
} # table
},


cmap_feature_correspondence => {
table=>q{
create table cmap_feature_correspondence (
  feature_correspondence_id int(11) NOT NULL default '0',
  accession_id varchar(20) NOT NULL default '',
  feature_id1 int(11) NOT NULL default '0',
  feature_id2 int(11) NOT NULL default '0',
  is_enabled tinyint(4) NOT NULL default '1',
  PRIMARY KEY  (feature_correspondence_id),
  UNIQUE KEY accession_id (accession_id),
  KEY feature_id1 (feature_id1),
  KEY cmap_feature_corresp_idx (is_enabled,feature_correspondence_id)
) TYPE=MyISAM;
} # table
},


cmap_map => {
table=>q{
create table cmap_map (
  map_id int(11) NOT NULL default '0',
  accession_id varchar(20) NOT NULL default '',
  map_set_id int(11) NOT NULL default '0',
  map_name varchar(32) NOT NULL default '',
  display_order int(11) NOT NULL default '1',
  start_position double(11,2) default NULL,
  stop_position double(11,2) default NULL,
  PRIMARY KEY  (map_id),
  UNIQUE KEY accession_id (accession_id),
  UNIQUE KEY map_id (map_id,map_set_id,map_name,accession_id),
  KEY map_set_id_index (map_set_id)
) TYPE=MyISAM;
} # table
},


cmap_map_set => {
table=>q{
create table cmap_map_set (
  map_set_id int(11) NOT NULL default '0',
  accession_id varchar(20) NOT NULL default '',
  map_set_name varchar(64) NOT NULL default '',
  short_name varchar(30) NOT NULL default '',
  map_type_accession varchar(20) NOT NULL default '0',
  species_id int(11) NOT NULL default '0',
  published_on date default NULL,
  can_be_reference_map tinyint(4) NOT NULL default '1',
  display_order int(11) NOT NULL default '1',
  is_enabled tinyint(4) NOT NULL default '1',
  shape varchar(12) default NULL,
  color varchar(20) default NULL,
  width int(11) default NULL,
  map_units varchar(12) NOT NULL default '',
  is_relational_map tinyint(11) NOT NULL default '0',
  PRIMARY KEY  (map_set_id),
  UNIQUE KEY accession_id (accession_id),
  UNIQUE KEY map_set_id (map_set_id,species_id,short_name,accession_id),
  KEY cmap_map_set_idx (can_be_reference_map,is_enabled,species_id,display_order,published_on,short_name)
) TYPE=MyISAM;
} # table
},


cmap_next_number => {
table=>q{
create table cmap_next_number (
  table_name varchar(40) NOT NULL default '',
  next_number int(11) NOT NULL default '0',
  PRIMARY KEY  (table_name)
) TYPE=MyISAM;
}, # table
insert=>{next_num=>q[ insert into cmap_next_number (table_name,next_number) VALUES ('cmap_feature',82);]}
},


cmap_species => {
table=>q{
create table cmap_species (
  species_id int(11) NOT NULL default '0',
  accession_id varchar(20) NOT NULL default '',
  common_name varchar(64) NOT NULL default '',
  full_name varchar(64) NOT NULL default '',
  display_order int(11) NOT NULL default '1',
  PRIMARY KEY  (species_id),
  KEY acc_id_species_id (accession_id,species_id)
) TYPE=MyISAM;
} # table
},


cmap_xref => {
table=>q{
create table cmap_xref (
  xref_id int(11) NOT NULL default '0',
  table_name varchar(30) NOT NULL default '',
  object_id int(11) default NULL,
  display_order int(11) NOT NULL default '1',
  xref_name varchar(200) NOT NULL default '',
  xref_url text NOT NULL,
  PRIMARY KEY  (xref_id),
  KEY table_name (table_name,object_id,display_order)
) TYPE=MyISAM;
} # table
},


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
  return 'SELECT DISTINCT gclass FROM cmap_feature WHERE NOT ISNULL(gclass)';
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

#xx1
  my $lookup_type = $dbh->prepare_delayed('SELECT ftypeid FROM ftype WHERE fmethod=? AND fsource=?');
  my $insert_type = $dbh->prepare_delayed('INSERT INTO ftype (fmethod,fsource) VALUES (?,?)');

  my $lookup_group = $dbh->prepare_delayed('SELECT feature_id FROM cmap_feature WHERE feature_name=? AND gclass=?');
  my $insert_group = $dbh->prepare_delayed(' INSERT into cmap_feature (feature_id, accession_id,feature_name, gclass ) VALUES (?,feature_id,?,?)');
  my $aux_insert_group = $dbh->prepare_delayed(' update cmap_next_number set next_number = next_number +1 where table_name=\'cmap_feature\'');
  my $next_id_group = $dbh->prepare_delayed('select next_number from cmap_next_number where table_name=\'cmap_feature\'');

  my $lookup_attribute = $dbh->prepare_delayed('SELECT fattribute_id FROM fattribute WHERE fattribute_name=?');
  my $insert_attribute = $dbh->prepare_delayed('INSERT INTO fattribute (fattribute_name) VALUES (?)');
  my $insert_attribute_value = $dbh->prepare_delayed('INSERT INTO fattribute_to_feature (fid,fattribute_id,fattribute_value) VALUES (?,?,?)');

  my $insert_data  = $dbh->prepare_delayed(<<END);
INSERT INTO fdata (fref,fstart,fstop,fbin,ftypeid,fscore,
		   fstrand,fphase,feature_id,ftarget_start,ftarget_stop)
       VALUES(?,?,?,?,?,?,?,?,?,?,?)
END
;


  $self->{load_stuff}{sth}{lookup_ftype}     = $lookup_type;
  $self->{load_stuff}{sth}{insert_ftype}     = $insert_type;
  #$self->{load_stuff}{sth}{lookup_fgroup}    = $lookup_group;
  #$self->{load_stuff}{sth}{insert_fgroup}    = $insert_group;
  $self->{load_stuff}{sth}{lookup_cmap_feature}     = $lookup_group;
  $self->{load_stuff}{sth}{insert_cmap_feature}     = $insert_group;
  $self->{load_stuff}{sth}{aux_insert_cmap_feature} = $aux_insert_group;
  $self->{load_stuff}{sth}{next_id_cmap_feature}   = $next_id_group;
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
  defined(my $groupid = $self->get_table_id('cmap_feature',$gff->{gname}  => $gff->{gclass})) or return;

  if ($gff->{stop}-$gff->{start}+1 > $self->max_bin) {
    warn "$gff->{gclass}:$gff->{gname} is longer than ",$self->maxbin,".\n";
    warn "Please set the maxbin value to a larger length than the largest feature you wish to store.\n";
    warn "With the command-line tools you do with this with --maxfeature option.\n";
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
      $dbh->prepare_delayed('SELECT fid FROM fdata WHERE fref=? AND fstart=? AND fstop=? AND ftypeid=? AND feature_id=?');
  }
  my $sth = $s->{get_feature_id} or return;
  $sth->execute($ref,$start,$stop,$typeid,$groupid) or return;
  my ($fid) = $sth->fetchrow_array;
  return $fid;
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
      if (defined($sth->{"next_id_$table"})){

        $sth->{"insert_$table"}->execute(3,'string1','string2');
        # Can't use auto incrementing
        $sth->{"next_id_$table"}->execute();
        $s->{$table}{$id_key} = ($sth->{"next_id_$table"}->fetchrow_array)[0];
        if ($s->{$table}{$id_key}){
            $sth->{"insert_$table"}->execute($s->{$table}{$id_key},@ids);
            $sth->{"aux_insert_$table"}->execute() if $sth->{"aux_insert_$table"};
        }
      }
      else{
          $sth->{"insert_$table"}->execute(@ids);
          $s->{$table}{$id_key} = $self->insertid($sth->{"insert_$table"}) unless $s->{$table}{$id_key};
          $sth->{"aux_insert_$table"}->execute() if $sth->{"aux_insert_$table"};
      }
    }
  }

  my $id = $s->{$table}{$id_key};
  unless (defined $id) {
    warn "No $table id for $id_key ",$dbh->errstr," Record skipped.\n";
    return;
  }
  $id;
}



#-----------------------------------

=head2 make_features_by_name_where_part

 Title   : make_features_by_name_where_part
 Usage   : $db->make_features_by_name_where_part
 Function: create the SQL fragment needed to select a feature by its group name & class
 Returns : a SQL fragment and bind arguments
 Args    : see below
 Status  : Protected

=cut

sub make_features_by_name_where_part {
  my $self = shift;
  my ($class,$name) = @_;
  if ($name =~ /\*/) {
    $name =~ tr/*/%/;
    return ("cmap_feature.gclass=? AND cmap_feature.feature_name LIKE ?",$class,$name);
  } else {
    return ("cmap_feature.gclass=? AND cmap_feature.feature_name=?",$class,$name);
  }
}

=head2 make_features_join_part

 Title   : make_features_join_part
 Usage   : $string = $db->make_features_join_part()
 Function: make join part of the features query
 Returns : a string
 Args    : none
 Status  : protected

This method creates the part of the features query that immediately
follows the WHERE keyword.

=cut

sub make_features_join_part {
  my $self = shift;
  my $options = shift || {};
  return !$options->{attributes} ? <<END1 : <<END2;
  cmap_feature.feature_id = fdata.feature_id
  AND ftype.ftypeid = fdata.ftypeid
END1
  cmap_feature.feature_id = fdata.feature_id
  AND ftype.ftypeid = fdata.ftypeid
  AND fattribute.fattribute_id=fattribute_to_feature.fattribute_id
  AND fdata.fid=fattribute_to_feature.fid
END2
}

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

  my @words  = $search_string =~ /(\w+)/g;
  my $regex  = join '|',@words;
  my @searches = map {"fattribute_value LIKE '%${_}%'"} @words;
  my $search   = join(' OR ',@searches);

  my $query = <<END;
SELECT distinct gclass,feature_name as gname,fattribute_value
  FROM cmap_feature,fattribute_to_feature,fdata
  WHERE cmap_feature.feature_id=fdata.feature_id
     AND fdata.fid=fattribute_to_feature.fid
END
;
  $query .= " AND ($search) " if ($search);

  my $sth = $self->dbh->do_query($query);
  my @results;
  while (my ($class,$name,$note) = $sth->fetchrow_array) {
     next unless $class && $name;    # sorry, ignore NULL objects
     my @matches = $note =~ /($regex)/g;
     my $relevance = 10*@matches;
     my $featname = Bio::DB::GFF::Featname->new($class=>$name);
     push @results,[$featname,$note,$relevance];
     last if $limit && @results >= $limit;
  }
  @results;
}

# sub search_notes {
#   my $self = shift;
#   my ($search_string,$limit) = @_;
#   my $query = FULLTEXTSEARCH;
#   $query .= " limit $limit" if defined $limit;
#   my $sth = $self->dbh->do_query($query,$search_string,$search_string);
#   my @results;
#   while (my ($class,$name,$note,$relevance) = $sth->fetchrow_array) {
#      next unless $class && $name;    # sorry, ignore NULL objects
#      $relevance = sprintf("%.2f",$relevance);  # trim long floats
#      my $featname = Bio::DB::GFF::Featname->new($class=>$name);
#      push @results,[$featname,$note,$relevance];
#   }
#   @results;
# }

=head2 make_features_order_by_part

 Title   : make_features_order_by_part
 Usage   : ($query,@args) = $db->make_features_order_by_part()
 Function: make the ORDER BY part of the features() query
 Returns : a SQL fragment and bind arguments, if any
 Args    : none
 Status  : protected

This method creates the part of the features query that immediately
follows the ORDER BY part of the query issued by features() and
related methods.

=cut

sub make_features_order_by_part {
  my $self = shift;
  my $options = shift || {};
  return "cmap_feature.feature_name";
}

=head2 create_cmap_viewer_link

 Title   : create_cmap_viewer_link
 Usage   : $link_str = $db->create_cmap_viewer_link(data_source=>$ds,group_id=>$gid)
 Function: 
 Returns : 
 Args    : 
 Status  : 


=cut

sub create_cmap_viewer_link {
  my $self = shift;
  my %args = @_;
  my $data_source = $args{'data_source'};
  my $gid         = $args{'group_id'};
  my $link_str    = undef;

  my $db = $self->features_db;
  my $sql_str = qq[
    select f.feature_name, 
        f.feature_type_accession feature_type_aid,
        m.accession_id as map_aid,
        ms.accession_id as map_set_aid 
    from cmap_feature f, 
        cmap_map m, 
        cmap_map_set ms 
    where f.map_id=m.map_id 
        and ms.map_set_id=m.map_set_id 
        and f.feature_id=$gid
    ];

  my $result_ref = $db->selectrow_hashref($sql_str,{ Columns => {} });
  
  if ( $result_ref ) {
    $link_str='/cgi-bin/cmap/viewer?ref_map_set_aid='
      . $result_ref->{'map_set_aid'}
      . '&ref_map_aids='
      . $result_ref->{'map_aid'}
      . '&data_source='
      . $data_source
      . '&highlight='
      .$result_ref->{'feature_name'}
      . '&feature_type_'
      .$result_ref->{'feature_type_aid'}
      . '=2';
  }

  return $link_str;
}

1;


__END__

=head1 BUGS

none ;-)

=head1 SEE ALSO

L<Bio::DB::GFF>, L<bioperl>

=head1 AUTHOR

Ben Faga E<lt>faga@cshl.orgE<gt>.

Modified from mysql.pm by:

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2002 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

