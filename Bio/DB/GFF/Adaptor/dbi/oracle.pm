package Bio::DB::GFF::Adaptor::dbi::oracle;

=head1 NAME

Bio::DB::GFF::Adaptor::dbi::oracle -- Database adaptor for a specific oracle schema

=head1 SYNOPSIS

See L<Bio::DB::GFF>

=cut

# a simple oracle adaptor
use strict;
#use Bio::DB::GFF::Adaptor::dbi::mysql;
#use Bio::DB::GFF::Adaptor::dbi::mysqlopt;
use Bio::DB::GFF::Util::Binning;
use Bio::DB::GFF::Util::Rearrange; # for rearrange()
use base qw(Bio::DB::GFF::Adaptor::dbi);

use constant MAX_SEGMENT => 100_000_000;  # the largest a segment can get
use constant DEFAULT_CHUNK => 2000;

use constant GETSEQCOORDS =><<END;
SELECT fref,
       NVL(gclass,'Sequence'),
       min(fstart),
       max(fstop),
       fstrand,
       gname
  FROM fdata,fgroup
  WHERE fgroup.gname=?
    AND fgroup.gclass=?
    AND fgroup.gid=fdata.gid
    GROUP BY fref,fstrand,gclass,gname
END
;

use constant GETALIASCOORDS =><<END;
SELECT fref,
       NVL(gclass,'Sequence'),
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
    GROUP BY fref,fstrand,gclass,gname
END
;

use constant GETALIASLIKE =><<END;
SELECT fref,
       NVL(gclass,'Sequence'),
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
       NVL(gclass,'Sequence'),
       min(fstart),
       max(fstop),
       fstrand
  FROM fdata,fgroup
  WHERE fgroup.gname=?
    AND fgroup.gclass=?
    AND fdata.fref=?
    AND fgroup.gid=fdata.gid
    GROUP BY fref,fstrand,gclass
END
;

########################
# moved from mysqlopt.pm
########################

# this is the largest that any reference sequence can be (100 megabases)
use constant MAX_BIN    => 100_000_000;

# this is the smallest bin (1 K)
use constant MIN_BIN    => 1000;

# size of range over which it is faster to force mysql to use the range for indexing
use constant STRAIGHT_JOIN_LIMIT => 200_000;

##############################################################################

=head1 DESCRIPTION

This adaptor implements a specific oracle database schema that is
compatible with Bio::DB::GFF.  It inherits from
Bio::DB::GFF::Adaptor::dbi, which itself inherits from Bio::DB::GFF.

The schema uses several tables:

=over 4

=item fdata

This is the feature data table.  Its columns are:

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

  sql> select * from fgroup where gname='sjj_2L52.1';
  +-------+-------------+------------+
  | gid   | gclass      | gname      |
  +-------+-------------+------------+
  | 69736 | PCR_product | sjj_2L52.1 |
  +-------+-------------+------------+
  1 row in set (0.70 sec)

  sql> select fref,fstart,fstop from fdata,fgroup 
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

  sql> select fref,fstart,fstop,fmethod,fsource from fdata,fgroup,ftype 
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
except that the prefix "dbi:oracle" is added to the database DSN identifier
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
  $dsn = "dbi:Oracle:$dsn" if !ref($dsn) && $dsn !~ /^(dbi|DBI):/;
  my $self = $class->SUPER::new(-dsn=>$dsn,%$other);
  $self;
}

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
create table fdata (
  fid INTEGER  NOT NULL,
  fref VARCHAR(100) DEFAULT '' NOT NULL,
  fstart INTEGER DEFAULT '0' NOT NULL,
  fstop INTEGER DEFAULT '0' NOT NULL,
  fbin NUMBER DEFAULT '0.000000' NOT NULL,
  ftypeid INTEGER DEFAULT '0' NOT NULL,
  fscore NUMBER  ,
  fstrand VARCHAR2(3)   CHECK (fstrand IN ('+','-')),
  fphase VARCHAR2(3)   CHECK (fphase IN ('0','1','2')),
  gid INTEGER DEFAULT '0' NOT NULL,
  ftarget_start INTEGER  ,
  ftarget_stop INTEGER  ,
  CONSTRAINT fdata_pk PRIMARY KEY (fid)
)
}, # fdata table

index=>{
		fdata_fref_idx => q{
CREATE UNIQUE INDEX fdata_fref_idx ON fdata (fref,fbin,fstart,fstop,ftypeid,gid)
},
	
		fdata_ftypeid_idx => q{
CREATE INDEX fdata_ftypeid_idx ON fdata (ftypeid)
},

		fdata_gid_idx => q{
CREATE  INDEX fdata_gid_idx ON fdata (gid)
}
	 }, # fdata indexes

sequence=> {
		fdata_fid_sq => q{
CREATE SEQUENCE fdata_fid_sq START WITH 1
}
	    }, # fdata sequences

trigger=> {
		fdata_fid_ai => q{
CREATE OR REPLACE TRIGGER fdata_fid_ai
BEFORE INSERT ON fdata
FOR EACH ROW WHEN (new.fid IS NULL OR new.fid = 0)
BEGIN
   SELECT fdata_fid_sq.nextval INTO :new.fid FROM dual;
END;
}
	   }# fdata triggers
			
}, # fdata



		fgroup => { 
table => q{
CREATE TABLE fgroup (
  gid INTEGER  NOT NULL,
  gclass VARCHAR(100)  ,
  gname VARCHAR(100)  ,
  CONSTRAINT fgroup_pk PRIMARY KEY (gid)
)
}, # fgroup table

index => {
		fgroup_gclass_idx => q{
CREATE UNIQUE INDEX fgroup_gclass_idx ON fgroup (gclass,gname)
}
	   }, # fgroup indexes

sequence => {

		fgroup_gid_sq => q{
CREATE SEQUENCE fgroup_gid_sq START WITH 1
}
	     }, # fgroup sequences


trigger => {
		fgroup_gid_ai => q{
CREATE OR REPLACE TRIGGER fgroup_gid_ai
BEFORE INSERT ON fgroup
FOR EACH ROW WHEN (new.gid IS NULL OR new.gid = 0)
BEGIN
   SELECT fgroup_gid_sq.nextval INTO :new.gid FROM dual;
END;
}
	    } # fgroup triggers

}, # fgroup

		ftype => { 
table => q{
CREATE TABLE ftype (
  ftypeid INTEGER  NOT NULL,
  fmethod VARCHAR(100) DEFAULT '' NOT NULL,
  fsource VARCHAR(100),
  CONSTRAINT ftype_pk PRIMARY KEY (ftypeid)
)
}, # ftype table

index => {
		ftype_fmethod_idx => q{
CREATE  INDEX ftype_fmethod_idx ON ftype (fmethod)
},

		ftype_fsource_idx => q{
CREATE  INDEX ftype_fsource_idx ON ftype (fsource)
},
	
		ftype_ftype_idx => q{
CREATE UNIQUE INDEX ftype_ftype_idx ON ftype (fmethod,fsource)
}
	   }, # ftype indexes

sequence => {
		ftype_ftypeid_sq => q{
CREATE SEQUENCE ftype_ftypeid_sq START WITH 1
}
	     }, #ftype sequences

trigger => {
		ftype_ftypeid_ai => q{
CREATE OR REPLACE TRIGGER ftype_ftypeid_ai
BEFORE INSERT ON ftype
FOR EACH ROW WHEN (new.ftypeid IS NULL OR new.ftypeid = 0)
BEGIN
   SELECT ftype_ftypeid_sq.nextval INTO :new.ftypeid FROM dual;
END;
}
	    } #ftype triggers
}, # ftype


         fdna => {
table => q{
CREATE TABLE fdna (
  fref VARCHAR(100) DEFAULT '' NOT NULL,
  foffset INTEGER DEFAULT '0' NOT NULL,
  fdna LONG         /* LONGBLOB */  ,
  CONSTRAINT fdna_pk PRIMARY KEY (fref,foffset)
)
} #fdna table
		 }, #fdna 

        fmeta => {
table => q{
CREATE TABLE fmeta (
  fname VARCHAR(255) DEFAULT '' NOT NULL,
  fvalue VARCHAR(255) DEFAULT '' NOT NULL,
  CONSTRAINT fmeta_pk PRIMARY KEY (fname)
)
} # fmeta table
		 }, # fmeta


       fattribute => {
table => q{
CREATE TABLE fattribute (
  fattribute_id INTEGER  NOT NULL,
  fattribute_name VARCHAR(255) DEFAULT '' NOT NULL,
  CONSTRAINT fattribute_pk PRIMARY KEY (fattribute_id)
)
}, # fattribute table

sequence=> {
       fattribute_fattribute_id_sq => q{
CREATE SEQUENCE fattribute_fattribute_id_sq START WITH 1
}
	    }, # fattribute sequences

trigger => {
       fattribute_fattribute_id_ai => q{
CREATE OR REPLACE TRIGGER fattribute_fattribute_id_ai
BEFORE INSERT ON fattribute
FOR EACH ROW WHEN (new.fattribute_id IS NULL OR new.fattribute_id = 0)
BEGIN
   SELECT fattribute_fattribute_id_sq.nextval INTO :new.fattribute_id FROM dual;
END;
}
	    } # fattribute triggers
}, # fattribute

       fattribute_to_feature => {
table => q{
CREATE TABLE fattribute_to_feature (
  fid INTEGER DEFAULT '0' NOT NULL,
  fattribute_id INTEGER DEFAULT '0' NOT NULL,
  fattribute_value VARCHAR2(255) /* TEXT */  
)
}, # fattribute_to_feature table

index => {
       fattribute_to_feature_fid => q{
CREATE  INDEX fattribute_to_feature_fid ON fattribute_to_feature (fid,fattribute_id)
}
	   } # fattribute_to_feature indexes
}, # fattribute_to_feature  

       finterval_stats => {
table=> q{
CREATE TABLE "finterval_stats" (
   "ftypeid"          integer DEFAULT '0' NOT NULL,
   "fref"             VARCHAR(100) DEFAULT '' NOT NULL,
   "fbin"             integer DEFAULT '0' NOT NULL,
   "fcum_count"       integer DEFAULT '0' NOT NULL,
   CONSTRAINT finterval_stats_pk PRIMARY KEY (ftypeid,fref,fbin)
)
} # finterval_stats table
},# finterval_stats

);
  return \%schema;
}


=head2 do_initialize

 Title   : do_initialize
 Usage   : $success = $db->do_initialize($drop_all)
 Function: initialize the database
 Returns : a boolean indicating the success of the operation
 Args    : a boolean indicating whether to delete existing data
 Status  : protected

This method will load the schema into the database.  If $drop_all is
true, then any existing data in the tables known to the schema will be
deleted.

Internally, this method calls schema() to get the schema data.

=cut

# Create the schema from scratch.
# You will need create privileges for this.
#sub do_initialize {
#  my $self = shift;
#  my $erase = shift;
#  $self->drop_all if $erase;

#  my $dbh = $self->features_db;
#  my $schema = $self->schema;
 
#  foreach my $table_name(keys %$schema) {
#    my $create_table_stmt = $$schema{$table_name}{table} ;
#    $dbh->do($create_table_stmt) ||  warn $dbh->errstr;    
#  }
#  1;
#}



=head2 drop_all

 Title   : drop_all
 Usage   : $db->drop_all
 Function: empty the database
 Returns : void
 Args    : none
 Status  : protected

This method drops the tables known to this module.  Internally it
calls the abstract tables() method.

=cut

# Drop all the GFF tables -- dangerous!
#sub drop_all {
#  my $self = shift;
#  my $dbh = $self->features_db;
#  local $dbh->{PrintError} = 0;
#  foreach ($self->tables) {
#    $dbh->do("drop table $_");
#  }
#}






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
  my $schema = $self->schema; 

  my $dbh = $self->features_db;

  if ($self->lock_on_load) {
    my @tables = map { "$_ WRITE"} $self->tables;
    my $tables = join ', ',@tables;
    $dbh->do("LOCK TABLES $tables");
  }

  my $lookup_type = $dbh->prepare_delayed('SELECT ftypeid FROM ftype WHERE fmethod=? AND fsource=?');
  my $insert_type = $dbh->prepare_delayed('INSERT INTO ftype (fmethod,fsource) VALUES (?,?)');
  my $sequence_type = (keys %{$schema->{ftype}{sequence}})[0];
  my $insertid_type = $dbh->prepare_delayed("SELECT $sequence_type.CURRVAL FROM dual");

  my $lookup_group = $dbh->prepare_delayed('SELECT gid FROM fgroup WHERE gname=? AND gclass=?');
  my $insert_group = $dbh->prepare_delayed('INSERT INTO fgroup (gname,gclass) VALUES (?,?)');
  my $sequence_group = (keys %{$schema->{fgroup}{sequence}})[0];
  my $insertid_group = $dbh->prepare_delayed("SELECT $sequence_group.CURRVAL FROM dual");

  my $lookup_attribute = $dbh->prepare_delayed('SELECT fattribute_id FROM fattribute WHERE fattribute_name=?');
  my $insert_attribute = $dbh->prepare_delayed('INSERT INTO fattribute (fattribute_name) VALUES (?)');
  my $sequence_attribute = (keys %{$schema->{fattribute}{sequence}})[0];
  my $insertid_attribute = $dbh->prepare_delayed("SELECT $sequence_attribute.CURRVAL FROM dual");

  my $insert_attribute_value = $dbh->prepare_delayed('INSERT INTO fattribute_to_feature (fid,fattribute_id,fattribute_value) VALUES (?,?,?)');

  my $insert_data  = $dbh->prepare_delayed(<<END);
INSERT INTO fdata (fref,fstart,fstop,fbin,ftypeid,fscore,
		   fstrand,fphase,gid,ftarget_start,ftarget_stop)
       VALUES(?,?,?,?,?,?,?,?,?,?,?)
END
;
  my $delete_existing_data = $dbh->prepare_delayed('DELETE FROM fdata WHERE fref=? AND fstart=? AND fstop=? AND fbin=? AND ftypeid=? AND GID=?');
  my $sequence_data =  (keys %{$schema->{fdata}{sequence}})[0];
  my $insertid_data = $dbh->prepare_delayed("SELECT $sequence_data.CURRVAL FROM dual");



  $self->{load_stuff}{sth}{lookup_ftype}     = $lookup_type;
  $self->{load_stuff}{sth}{insert_ftype}     = $insert_type;
  $self->{load_stuff}{sth}{insertid_ftype}   = $insertid_type;
  $self->{load_stuff}{sth}{lookup_fgroup}    = $lookup_group;
  $self->{load_stuff}{sth}{insert_fgroup}    = $insert_group;
  $self->{load_stuff}{sth}{insertid_fgroup}  = $insertid_group;
  $self->{load_stuff}{sth}{insert_fdata}     = $insert_data;
  $self->{load_stuff}{sth}{insertid_fdata}   = $insertid_data;
  $self->{load_stuff}{sth}{delete_existing_fdata} = $delete_existing_data;
  $self->{load_stuff}{sth}{lookup_fattribute} = $lookup_attribute;
  $self->{load_stuff}{sth}{insert_fattribute} = $insert_attribute;
  $self->{load_stuff}{sth}{insertid_fattribute} = $insertid_attribute;
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
  
  if (defined $gff->{phase}){
     chomp($gff->{phase}); 
     undef($gff->{phase}) if $gff->{phase} eq '.';
   }

  if (defined $gff->{strand} && $gff->{strand} eq '.'){undef($gff->{strand})}; 
  if (defined $gff->{score}  && $gff->{score} eq '.'){undef($gff->{score})};

  my $s    = $self->{load_stuff};
  my $dbh  = $self->features_db;
  local $dbh->{PrintError} = 0;

  defined(my $typeid  = $self->get_table_id('ftype', $gff->{method} => $gff->{source})) or return;
  defined(my $groupid = $self->get_table_id('fgroup',$gff->{gname}  => $gff->{gclass})) or return;

  my $bin =  bin($gff->{start},$gff->{stop},$self->min_bin);
  my $result = $s->{sth}{insert_fdata}->execute($gff->{ref},
					       $gff->{start},$gff->{stop},$bin,
					       $typeid,
					       $gff->{score},$gff->{strand},$gff->{phase},
					       $groupid,
					       $gff->{tstart},$gff->{tstop});
  if (defined ($dbh->errstr)){
    print  $dbh->errstr,"\n" ,%$gff,"\n";
    if ($dbh->errstr =~ /ORA-02290: check constraint/){
      print "PHASE=$gff->{phase}"."===","\n";
    }

    if ($dbh->errstr =~ /ORA-00001: unique constraint/){
      $result = $s->{sth}{delete_existing_fdata}->execute($gff->{ref},
    							   $gff->{start},$gff->{stop},$bin,
    							   $typeid,
    							   $groupid);
    
      print "delete row result=$result\n";
      $result = $s->{sth}{insert_fdata}->execute($gff->{ref},
    					       $gff->{start},$gff->{stop},$bin,
    					       $typeid,
    					       $gff->{score},$gff->{strand},$gff->{phase},
    					       $groupid,
    					       $gff->{tstart},$gff->{tstop}); 
    
      print "insert row result=$result\n";
    }
  }
  warn $dbh->errstr,"\n" and print "ref=",$gff->{ref}," start=",$gff->{start}," stop=",$gff->{stop}," bin=",$bin," typeid=",$typeid," groupid=",$groupid,"\n" 
    and return unless $result;
  
  my $fid = $self->insertid($s->{sth},'fdata')
    || $self->get_feature_id($gff->{ref},$gff->{start},$gff->{stop},$typeid,$groupid);


  # insert attributes

  #  print STDERR map {"$fid attribute:". $_->[0]."=".$_->[1]."\n"} @{$gff->{attributes}};

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
    $sth->{"lookup_$table"}->execute(@ids);
    my @result = $sth->{"lookup_$table"}->fetchrow_array;
    if (@result > 0) {
      $s->{$table}{$id_key} = $result[0];
    } else {
      $sth->{"insert_$table"}->execute(@ids)
	&& ($s->{$table}{$id_key} = $self->insertid($sth,$table));
	#&& ($s->{$table}{$id_key} = $self->insertid($sth->{"insertid_$table"}));
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
  my $sth = shift ;
  my $table = shift;

  my $insert_id;
  if ($sth->{"insertid_$table"}->execute()){
     $insert_id = ($sth->{"insertid_$table"}->fetchrow_array)[0];
  }
  else{
    warn "No CURRVAL for SEQUENCE of table $table ",$sth->errstr,"\n";
    return;
  }
  return $insert_id;
}


#sub insertid {
#  my $self = shift;
#  my $insertid_sth = shift ;
#  my $insert_id;
#  if ($insertid_sth->execute){
#     $insert_id = ($insertid_sth->fetchrow_array)[0];
#  }
#  else{
#    warn "No CURRVAL for SEQUENCE  ",$insertid_sth->errstr,"\n";
#    return;
#  }
#  return $insert_id;
#}

sub insert_sequence {
  my $self = shift;
  my($id,$offset,$seq) = @_;
  my $sth = $self->{_insert_sequence}
    ||= $self->dbh->prepare_delayed('insert into fdna values (?,?,?)');
  $sth->execute($id,$offset,$seq) or $self->throw($sth->errstr);
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
  column 4     A Bio::DB::GFF::Typename object

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
SELECT distinct gclass,gname,fattribute_value,fmethod,fsource
  FROM fgroup,fattribute_to_feature,fdata,ftype
  WHERE fgroup.gid=fdata.gid
     AND fdata.fid=fattribute_to_feature.fid
     AND fdata.ftypeid=ftype.ftypeid
     AND ($search)
END
;

  my $sth = $self->dbh->do_query($query);
  my @results;
  while (my ($class,$name,$note,$method,$source) = $sth->fetchrow_array) {
     next unless $class && $name;    # sorry, ignore NULL objects
     my @matches = $note =~ /($regex)/g;
     my $relevance = 10*@matches;
     my $featname = Bio::DB::GFF::Featname->new($class=>$name);
     my $type     = Bio::DB::GFF::Typename->new($method,$source);
     push @results,[$featname,$note,$relevance,$type];
     last if $limit && @results >= $limit;
  }
  @results;
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
   return 'INSERT INTO fmeta VALUES (?,?)';
}

sub make_classes_query {
  my $self = shift;
  return 'SELECT DISTINCT gclass FROM fgroup WHERE NOT gclass IS NULL';
}


sub chunk_size {
  my $self = shift;
  $self->meta('chunk_size') || DEFAULT_CHUNK;
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


sub make_features_select_part {
  my $self = shift;
  my $options = shift || {};
  my $s;
  if (my $b = $options->{bin_width}) {

    $s = <<END;
fref,
  1+$b*floor(fstart/$b)   as fstart,
  $b*(1+floor(fstart/$b)) as fstop,
  NVL2(fsource,fmethod||':'||fsource,fmethod),'bin',
  count(*) as fscore,
  '.','.','bin',
  NVL2(fsource , fref||':'||fmethod||':'||fsource , fref||':'||fmethod),
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

sub make_features_from_part_bkup {
  my $self = shift;
  my $sparse = shift;
  my $options = shift || {};
  #my $index = $sparse ? ' USE INDEX(ftypeid)': '';
  my $index =  '';
  return $options->{attributes} ? "fdata${index},ftype,fgroup,fattribute,fattribute_to_feature\n"
                                : "fdata${index},ftype,fgroup\n";
}


####################################
# moved from mysqlopt.pm
###################################
# meta values
sub default_meta_values {
  my $self = shift;
  my @values = $self->SUPER::default_meta_values;
  return (
	  @values,
	  max_bin => MAX_BIN,
	  min_bin => MIN_BIN,
	  straight_join_limit => STRAIGHT_JOIN_LIMIT,
	 );
}

sub min_bin {
  my $self = shift;
  return $self->meta('min_bin') || MIN_BIN;
}
sub max_bin {
  my $self = shift;
  return $self->meta('max_bin') || MAX_BIN;
}
sub straight_join_limit {
  my $self = shift;
  return $self->meta('straight_join_limit') || STRAIGHT_JOIN_LIMIT;
}

1;
