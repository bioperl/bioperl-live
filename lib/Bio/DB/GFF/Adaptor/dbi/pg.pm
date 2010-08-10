package Bio::DB::GFF::Adaptor::dbi::pg;

=head1 NAME

Bio::DB::GFF::Adaptor::dbi::pg -- Database adaptor for a specific postgres schema

=head1 NOTES 

SQL commands that need to be executed before this adaptor will work:

  CREATE DATABASE <dbname>;

Also, select permission needs to be granted for each table in the
database to the owner of the httpd process (usually 'nobody', but
for some RedHat systems it is 'apache') if this adaptor is to be used
with the Generic Genome Browser (gbrowse):

  CREATE USER nobody;
  GRANT SELECT ON TABLE fmeta                 TO nobody;
  GRANT SELECT ON TABLE fgroup                TO nobody;
  GRANT SELECT ON TABLE fdata                 TO nobody;
  GRANT SELECT ON TABLE fattribute_to_feature TO nobody;
  GRANT SELECT ON TABLE fdna                  TO nobody;
  GRANT SELECT ON TABLE fattribute            TO nobody;
  GRANT SELECT ON TABLE ftype                 TO nobody;

=head2 Optimizing the database

PostgreSQL generally requires some tuning before you get very good
performance for large databases.  For general information on tuning
a PostgreSQL server, see http://www.varlena.com/GeneralBits/Tidbits/perf.html
Of particular importance is executing VACUUM FULL ANALYZE whenever
you change the database.

Additionally, for a GFF database, there are a few items you can tune.
For each automatic class in your GBrowse conf file, there will be one
or two searches done when searching for a feature.  If there are lots 
of features, these search can take several seconds.  To speed these searches,
do two things:

=over

=item 1

Set 'enable_seqscan = false' in your postgresql.conf file (and restart
your server).

=item 2

Create 'partial' indexes for each automatic class, doing this for the
example class 'Allele':

  CREATE INDEX partial_allele_gclass ON 
    fgroup (lower('gname')) WHERE gclass='Allele';

And be sure to run VACUUM FULL ANALYZE after creating the indexes.

=back

=cut

# a simple postgres adaptor
use strict;
use Bio::DB::GFF::Util::Binning; 
use Bio::DB::GFF::Util::Rearrange; # for rearrange()
use base qw(Bio::DB::GFF::Adaptor::dbi);

use constant MAX_SEGMENT => 100_000_000;  # the largest a segment can get
use constant DEFAULT_CHUNK => 2000;

use constant GETSEQCOORDS =><<END;
SELECT fref,
       COALESCE(gclass,'Sequence'),
       min(fstart),
       max(fstop),
       fstrand,
       gname
  FROM fdata,fgroup
  WHERE lower(fgroup.gname) = lower(?)
    AND fgroup.gclass=?
    AND fgroup.gid=fdata.gid
    GROUP BY fref,fstrand,gclass,gname
END
;

use constant GETALIASCOORDS =><<END;
SELECT fref,
       COALESCE(gclass,'Sequence'),
       min(fstart),
       max(fstop),
       fstrand,
       gname
  FROM fdata,fgroup,fattribute,fattribute_to_feature
  WHERE lower(fattribute_to_feature.fattribute_value)=lower(?)
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
       COALESCE(gclass,'Sequence'),
       min(fstart),
       max(fstop),
       fstrand,
       gname
  FROM fdata,fgroup,fattribute,fattribute_to_feature
  WHERE lower(fattribute_to_feature.fattribute_value) LIKE lower(?)
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
       COALESCE(gclass,'Sequence'),
       min(fstart),
       max(fstop),
       fstrand
  FROM fdata,fgroup
  WHERE lower(fgroup.gname) = lower(?)
    AND fgroup.gclass=?
    AND lower(fdata.fref) = lower(?)
    AND fgroup.gid=fdata.gid
    GROUP BY fref,fstrand,gclass
END
;

use constant FULLTEXTWILDCARD => <<END;
SELECT distinct gclass,gname,fattribute_value
    FROM fgroup,fattribute_to_feature,fdata
     WHERE fgroup.gid=fdata.gid
       AND fdata.fid=fattribute_to_feature.fid
       AND lower(fattribute_to_feature.fattribute_value) LIKE lower(?)
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

This adaptor implements a specific postgres database schema that is
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
except that the prefix "dbi:pg" is added to the database DSN identifier
automatically if it is not there already.

  Argument       Description
  --------       -----------

  -dsn           the DBI data source, e.g. 'dbi:Pg:dbname=:ens0040' or "ens0040"

  -user          username for authentication

  -pass          the password for authentication

=cut

#'

sub new {
  my $class = shift;
  my ($dsn,$other) = rearrange([
				[qw(FEATUREDB DB DSN)],
			       ],@_);
  $dsn = "dbi:Pg:dbname=$dsn" if !ref($dsn) && $dsn !~ /^(dbi|DBI):/;
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
CREATE TABLE "fdata" (
  "fid" serial NOT NULL,
  "fref" character varying(100) DEFAULT '' NOT NULL,
  "fstart" integer DEFAULT '0' NOT NULL,
  "fstop" integer DEFAULT '0' NOT NULL,
  "fbin" double precision DEFAULT '0.000000' NOT NULL,
  "ftypeid" integer DEFAULT '0' NOT NULL,
  "fscore" double precision DEFAULT NULL,
  "fstrand" character varying(3) DEFAULT NULL,
  "fphase" character varying(3) DEFAULT NULL,
  "gid" integer DEFAULT '0' NOT NULL,
  "ftarget_start" integer DEFAULT NULL,
  "ftarget_stop" integer DEFAULT NULL,
  CONSTRAINT chk_fdata_fstrand CHECK (fstrand IN ('+','-')),
  CONSTRAINT chk_fdata_fphase CHECK (fphase IN ('0','1','2')),
  CONSTRAINT pk_fdata PRIMARY KEY (fid)
)
}, # fdata table

#CONSTRAINT fref_fdata UNIQUE (fref, fbin, fstart, fstop, ftypeid, gid)
# fdata_fref_idx => q{ CREATE UNIQUE INDEX fdata_fref_idx ON fdata (fref,fbin,fstart,fstop,ftypeid,gid)}, 

index=>{
                fdata_fref_idx => q{
CREATE INDEX fdata_fref_idx ON fdata (fref,fbin,fstart,fstop,ftypeid,gid)
},

		fdata_ftypeid_idx => q{
CREATE INDEX fdata_ftypeid_idx ON fdata (ftypeid)
},

		fdata_gid_idx => q{
CREATE INDEX fdata_gid_idx ON fdata (gid)
}
	 }, # fdata indexes

}, # fdata



		fgroup => { 
table => q{
CREATE TABLE "fgroup" (
  "gid" serial NOT NULL,
  "gclass" character varying(100) DEFAULT NULL,
  "gname" character varying(100) DEFAULT NULL,
  CONSTRAINT pk_fgroup PRIMARY KEY (gid)
)
}, # fgroup table

index => {
		fgroup_gclass_idx => q{
CREATE UNIQUE INDEX fgroup_gclass_idx ON fgroup (gclass,gname)
},
                fgroup_gname_idx => q{
CREATE INDEX fgroup_gname_idx ON fgroup(gname)
},
                fgroup_lower_gname_idx => q{
CREATE INDEX fgroup_lower_gname_idx ON fgroup (lower(gname))
},
	   }, # fgroup indexes

}, # fgroup

		ftype => { 
table => q{
CREATE TABLE "ftype" (
  "ftypeid" serial NOT NULL,
  "fmethod" character varying(100) DEFAULT '' NOT NULL,
  "fsource" character varying(100) DEFAULT NULL,
  CONSTRAINT pk_ftype PRIMARY KEY (ftypeid),
  CONSTRAINT ftype_ftype UNIQUE (fmethod, fsource)
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

}, # ftype


         fdna => {
table => q{
CREATE TABLE "fdna" (
  "fref" character varying(100) DEFAULT '' NOT NULL,
  "foffset" integer DEFAULT '0' NOT NULL,
  "fdna" bytea,
  CONSTRAINT pk_fdna PRIMARY KEY (fref, foffset)
)
} #fdna table
		 }, #fdna 

        fmeta => {
table => q{
CREATE TABLE "fmeta" (
  "fname" character varying(255) DEFAULT '' NOT NULL,
  "fvalue" character varying(255) DEFAULT '' NOT NULL,
  CONSTRAINT pk_fmeta PRIMARY KEY (fname)
)
} # fmeta table
		 }, # fmeta


       fattribute => {
table => q{
CREATE TABLE "fattribute" (
  "fattribute_id" serial NOT NULL,
  "fattribute_name" character varying(255) DEFAULT '' NOT NULL,
  CONSTRAINT pk_fattribute PRIMARY KEY (fattribute_id)
)
}, # fattribute table

}, # fattribute

       fattribute_to_feature => {
table => q{
CREATE TABLE "fattribute_to_feature" (
  "fid" integer DEFAULT '0' NOT NULL,
  "fattribute_id" integer DEFAULT '0' NOT NULL,
  "fattribute_value" text
)
}, # fattribute_to_feature table

index => {
       fattribute_to_feature_fid => q{
CREATE  INDEX fattribute_to_feature_fid ON fattribute_to_feature (fid,fattribute_id)
},
       fattribute_txt_idx => q{
CREATE INDEX fattribute_txt_idx ON fattribute_to_feature (fattribute_value)
},
       fattribute_lower_idx => q{
CREATE INDEX fattribute_lower_idx ON fattribute_to_feature (lower(fattribute_value))
},
	   } # fattribute_to_feature indexes
}, # fattribute_to_feature  

       finterval_stats => {
table=> q{
CREATE TABLE "finterval_stats" (
   "ftypeid"          integer DEFAULT '0' NOT NULL,
   "fref"             character varying(100) DEFAULT '' NOT NULL,
   "fbin"             integer DEFAULT '0' NOT NULL,
   "fcum_count"       integer DEFAULT '0' NOT NULL,
   CONSTRAINT pk_finterval_stats PRIMARY KEY (ftypeid,fref,fbin)
)
} # finterval_stats table
},# finterval_stats




);
  return \%schema;
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
  my $schema = $self->schema; 

  my $dbh = $self->features_db;

  if ($self->lock_on_load) {
    my @tables = map { "$_ WRITE"} $self->tables;
    my $tables = join ', ',@tables;
    $dbh->do("LOCK TABLES $tables");
  }

  my $lookup_type = $dbh->prepare_delayed('SELECT ftypeid FROM ftype WHERE fmethod=? AND fsource=?');
  my $insert_type = $dbh->prepare_delayed('INSERT INTO ftype (fmethod,fsource) VALUES (?,?)');
  my $insertid_type = $dbh->prepare_delayed("SELECT currval('ftype_ftypeid_seq')");

  my $lookup_group = $dbh->prepare_delayed('SELECT gid FROM fgroup WHERE lower(gname)=lower(?) AND gclass=?');
  my $insert_group = $dbh->prepare_delayed('INSERT INTO fgroup (gname,gclass) VALUES (?,?)');
  my $insertid_group = $dbh->prepare_delayed("SELECT currval('fgroup_gid_seq')");

  my $lookup_attribute = $dbh->prepare_delayed('SELECT fattribute_id FROM fattribute WHERE fattribute_name=?');
  my $insert_attribute = $dbh->prepare_delayed('INSERT INTO fattribute (fattribute_name) VALUES (?)');
  my $insertid_attribute = $dbh->prepare_delayed("SELECT currval('fattribute_fattribute_id_seq')");

  my $insert_attribute_value = $dbh->prepare_delayed('INSERT INTO fattribute_to_feature (fid,fattribute_id,fattribute_value) VALUES (?,?,?)');

  my $insert_data  = $dbh->prepare_delayed(<<END);
INSERT INTO fdata (fref,fstart,fstop,fbin,ftypeid,fscore,
		   fstrand,fphase,gid,ftarget_start,ftarget_stop)
       VALUES(?,?,?,?,?,?,?,?,?,?,?)
END
;
  my $delete_existing_data = $dbh->prepare_delayed('DELETE FROM fdata WHERE fref=? AND fstart=? AND fstop=? AND fbin=? AND ftypeid=? AND GID=?');
  my $insertid_data = $dbh->prepare_delayed("SELECT currval('fdata_fid_seq')");

  $self->{load_stuff}{sth}{lookup_ftype}     = $lookup_type;
  $self->{load_stuff}{sth}{insert_ftype}     = $insert_type;
  $self->{load_stuff}{sth}{insertid_ftype}   = $insertid_type;
  $self->{load_stuff}{sth}{lookup_fgroup}    = $lookup_group;
  $self->{load_stuff}{sth}{insert_fgroup}    = $insert_group;
  $self->{load_stuff}{sth}{insertid_fgroup}    = $insertid_group;
  $self->{load_stuff}{sth}{insertid_fdata}   = $insertid_data;
  $self->{load_stuff}{sth}{insert_fdata}     = $insert_data;
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

  warn $dbh->errstr,"\n" and print "ref=",$gff->{ref}," start=",$gff->{start}," stop=",$gff->{stop}," bin=",$bin," typeid=",$typeid," groupid=",$groupid,"\n" 
    and return unless $result;
 
  my $fid = $self->insertid($s->{sth},'fdata') 
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

=head2 range_query

 Title   : range_query
 Usage   : $db->range_query($range_type,$refseq,$refclass,$start,$stop,$types,$order_by_group,$attributes,$binsize)
 Function: create statement handle for range/overlap queries
 Returns : a DBI statement handle
 Args    : see below
 Status  : Protected

This method constructs the statement handle for this module's central
query: given a range and/or a list of feature types, fetch their GFF
records.  It overrides a method in dbi.pm so that the overlaps query
can write SQL optimized for Postgres.  Specifically, instead of writing
the bin related section as a set of ORs, each bin piece is place in 
a separate select and then they are UNIONed together.  This subroutine
requires several replacements for other subroutines in dbi.pm.  In this 
module, they are named the same as those in dbi.pm but prefixed with 
"pg_".

The positional arguments are as follows:

  Argument               Description

  $isrange               A flag indicating that this is a range.
                         query.  Otherwise an overlap query is
                         assumed.

  $refseq                The reference sequence name (undef if no range).

  $refclass              The reference sequence class (undef if no range).

  $start                 The start of the range (undef if none).

  $stop                  The stop of the range (undef if none).

  $types                 Array ref containing zero or feature types in the
                         format [method,source].

  $order_by_group        A flag indicating that statement handler should group
                         the features by group id (handy for iterative fetches)

  $attributes            A hash containing select attributes.

  $binsize               A bin size for generating tables of feature density.

=cut

sub range_query {
  my $self = shift;
  my($rangetype,$refseq,$class,$start,$stop,$types,$sparse,$order_by_group,$attributes,$bin) = @_;

  my $dbh = $self->features_db;

  #  my @bin_parts = split /\n\s+OR/, $self->bin_query($start,$stop);
  #  warn "bin_part: @bin_parts\n";

  my %a             = (refseq=>$refseq,class=>$class,start=>$start,stop=>$stop,types=>$types,attributes=>$attributes,bin_width=>$bin);
  my ($query, @args, $order_by);

  if ($rangetype ne 'overlaps') {

    my $select        = $self->make_features_select_part(\%a);
    my $from          = $self->make_features_from_part($sparse,\%a);
    my $join          = $self->make_features_join_part(\%a);
    my $where;
       ($where,@args) = $self->make_features_by_range_where_part($rangetype,\%a);
    my ($group_by,@more_args) = $self->make_features_group_by_part(\%a);
       $order_by      = $self->make_features_order_by_part(\%a) if $order_by_group;

    $query         = "SELECT $select FROM $from WHERE $join";
    $query           .= " AND $where" if $where;

    if ($group_by) {
      $query           .= " GROUP BY $group_by";
      push @args,@more_args;
    }

  } else {  # most common case: overlaps query

    my @bin_parts            = split /\s*OR/, $self->bin_query($start,$stop);
    my $select               = $self->make_features_select_part(\%a);
    my $from                 = $self->make_features_from_part($sparse,\%a);
    my $join                 = $self->make_features_join_part(\%a);
    my $where;
    ($where,@args)           = $self->pg_make_features_by_range_where_part($rangetype,\%a);
    my ($group_by,@more_args)= $self->make_features_group_by_part(\%a);
    $order_by                = $self->pg_make_features_order_by_part(\%a) if $order_by_group;

    my @temp_args;
    my @query_pieces; 
    foreach my $bin (@bin_parts) {
      my $temp_query = "SELECT $select FROM $from WHERE $join AND $where AND $bin\n"; 
      push @temp_args, @args;

      if ($group_by) {
        $temp_query    .= " GROUP BY $group_by";
        push @temp_args,@more_args;
      }

      push @query_pieces, $temp_query;
    }
    
    @args             = @temp_args;
    $query            = join("UNION\n", @query_pieces); 

  }

  $query           .= " ORDER BY $order_by" if $order_by;

  $self->dbh->do('set enable_seqscan=off');
  my $sth = $self->dbh->do_query($query,@args);
  $sth;
}

sub pg_make_features_by_range_where_part {
  my $self = shift;
  my ($rangetype,$options) = @_;

  return unless $rangetype eq 'overlaps';

  $options ||= {};
  my ($refseq,$class,$start,$stop,$types,$attributes) =
    @{$options}{qw(refseq class start stop types attributes)};

  my (@query,@args);

  if ($refseq) {
    my ($q,@a) = $self->refseq_query($refseq,$class);
    push @query,$q;
    push @args,@a;
  }

  if (defined $start or defined $stop) {
    $start = 0               unless defined($start);
    $stop  = MAX_SEGMENT     unless defined($stop);

    my ($range_query,@range_args) = $self->pg_overlap_query($start,$stop);

    push @query,$range_query;
    push @args,@range_args;
  }

  if (defined $types && @$types) {
    my ($type_query,@type_args) = $self->types_query($types);
    push @query,$type_query;
    push @args,@type_args;
  }

  if ($attributes) {
    my ($attribute_query,@attribute_args) = $self->make_features_by_attribute_where_part($attributes);
    push @query,"($attribute_query)";
    push @args,@attribute_args;
  }

  my $query = join "AND",@query;
  return wantarray ? ($query,@args) : $self->dbh->dbi_quote($query,@args);
}

sub pg_overlap_query {
  my $self = shift;
  my ($start,$stop) = @_;

  my ($iq,@iargs) = $self->overlap_query_nobin($start,$stop);
  my $query = "\n$iq\n";
  my @args = @iargs;

  return wantarray ? ($query,@args) : $self->dbh->dbi_quote($query,@args);
}

sub pg_make_features_order_by_part {
  my $self = shift;
  my $options = shift || {};
  return "gname";
}

=head2 search_notes

This PostgreSQL adaptor does not implement the search notes method
because it can be very slow (although the code for the method is
contained in this method but commented out).
There is, however, a PostgreSQL adaptor that does implement it in
a more efficient way: L<Bio::DB::GFF::Adaptor::dbi::pg_fts>,
which inherits from this adaptor and uses the optional PostgreSQL
module TSearch2 for full text indexing.  See that adaptor's
documentation for more information.

See also L<Bio::DB::GFF>

 Title   : search_notes
 Usage   : @search_results = $db->search_notes("full text search string",$limit)
 Function: Search the notes for a text string, using mysql full-text search
 Returns : array of results
 Args    : full text search string, and an optional row limit
 Status  : public

This is a replacement for the mysql-specific method.  Given a search string, it
performs a ILIKE search of the notes table and returns an array of results.
Each row of the returned array is a arrayref containing the following fields:

  column 1     A Bio::DB::GFF::Featname object, suitable for passing to segment()
  column 2     The text of the note
  column 3     A relevance score.

Note that for large databases this can be very slow and may result in
time out or 500-cgi errors.  If this is happening on a regular basis,
you should look into using L<Bio::DB::GFF::Adaptor::dbi::pg_fts> which
implements the TSearch2 full text indexing scheme.

=cut

sub search_notes{
#  my $self = shift;
#  my ($search_string,$limit) = @_;
#
#  $search_string =~ tr/*/%/s;
#  $search_string =  '%'.$search_string unless $search_string =~ /^\%/;
#  $search_string =  $search_string.'%' unless $search_string =~ /\%$/;
#  warn "search_string:$search_string";
#  my $query = FULLTEXTWILDCARD;
#  $query   .= " limit $limit" if defined $limit;
#  my $sth   = $self->dbh->do_query($query,$search_string);
#
#  my @results;
#  while (my ($class,$name,$note) = $sth->fetchrow_array) {
#
#     next unless $class && $name;    # sorry, ignore NULL objects
#     my $featname = Bio::DB::GFF::Featname->new($class=>$name);
#
#     push @results,[$featname,$note,0]; #gbrowse expects a score, but
#                                        #pg doesn't give one, thus the 0
#  }
#  warn @results;
#
#  return @results;
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
  CASE WHEN fsource IS NULL THEN fmethod
       ELSE fmethod||':'||fsource,
  'bin',
  count(*) as fscore,
  '.','.','bin',
  CASE WHEN fsource IS NULL THEN fref||':'||fmethod
       ELSE fref||':'||fmethod||':'||fsource,
  NULL,NULL,NULL,NULL
END
;
  } else {
    $s = <<END;
fref,fstart,fstop,fsource,fmethod,fscore,fstrand,fphase,gclass,fgroup.gname,ftarget_start,ftarget_stop,fdata.fid,fdata.gid
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


sub _feature_by_name {
  my $self = shift;
  my ($class,$name,$location,$callback) = @_;
  $callback || $self->throw('must provide a callback argument');

  my @bin_parts      = split /\s*OR/, $self->bin_query($location->[1],$location->[2]) if $location;
  my $select         = $self->make_features_select_part;
  my $from           = $self->make_features_from_part(undef,{sparse_groups=>1});
  my ($where,@args)  = $self->make_features_by_name_where_part($class,$name);
  my $join           = $self->make_features_join_part;
  my $range          = $self->pg_make_features_by_range_where_part('overlaps',
                                                                {refseq=>$location->[0],
                                                                 class =>'',
                                                                 start=>$location->[1],
                                                                 stop =>$location->[2]}) if $location;

  my @temp_args;
  my @query_pieces;
  my $query;
  if (@bin_parts) {
    foreach my $bin (@bin_parts) {
      my $temp_query = "SELECT $select FROM $from WHERE $join AND $where AND $range AND $bin\n";
      push @temp_args, @args;
      push @query_pieces, $temp_query;
    }

    @args  = @temp_args;
    $query = join("UNION\n", @query_pieces);

  } else {
    $query  = "SELECT $select FROM $from WHERE $where AND $join";
  }

  my $sth    = $self->dbh->do_query($query,@args);

  my $count = 0;
  while (my @row = $sth->fetchrow_array) {
    $callback->(@row);
    $count++;
  }
  $sth->finish;
  return $count;
}

sub update_sequences {
  my $self = shift;
  my $dbh  = $self->features_db;
 
  $dbh->do("SELECT setval('public.fdata_fid_seq', max(fid)+1) FROM fdata");
  $dbh->do("SELECT setval('public.fattribute_fattribute_id_seq', max(fattribute_id)+1) FROM fattribute");
  $dbh->do("SELECT setval('public.fgroup_gid_seq', max(gid)+1) FROM fgroup");
  $dbh->do("SELECT setval('public.ftype_ftypeid_seq', max(ftypeid)+1) FROM ftype");

  1;
}

=head2 make_features_by_name_where_part

 Title   : make_features_by_name_where_part
 Usage   : $db->make_features_by_name_where_part
 Function: Overrides a function in Bio::DB::GFF::Adaptor::dbi to insure
           that searches will be case insensitive. It creates the SQL
           fragment needed to select a feature by its group name & class
 Returns : a SQL fragment and bind arguments
 Args    : see below
 Status  : Protected

=cut

sub make_features_by_name_where_part {
  my $self = shift;
  my ($class,$name) = @_;

  if ($name !~ /\*/) {
    #allows utilization of an index on lower(gname)
    return ("fgroup.gclass=? AND lower(fgroup.gname) = lower(?)",$class,$name);
  }
  else {
    $name =~ tr/*/%/;
    return ("fgroup.gclass=? AND lower(fgroup.gname) LIKE lower(?)",$class,$name);
  }
}

#
# Methods from dbi.pm that need to be overridden to make
# searching for fref case insensitive
#
#
sub get_dna {
  my $self = shift;
  my ($ref,$start,$stop,$class) = @_;

  my ($offset_start,$offset_stop);

  my $has_start = defined $start;
  my $has_stop  = defined $stop;

  my $reversed;
  if ($has_start && $has_stop && $start > $stop) {
    $reversed++;
    ($start,$stop) = ($stop,$start);
  }

  # turn start and stop into 0-based offsets
  my $cs = $self->dna_chunk_size;
  $start -= 1;  $stop -= 1;
  $offset_start = int($start/$cs)*$cs;
  $offset_stop  = int($stop/$cs)*$cs;

  my $sth;
  # special case, get it all
  if (!($has_start || $has_stop)) {
    $sth = $self->dbh->do_query('select fdna,foffset from fdna where lower(fref)=lower(?) order by foffset',$ref);
  }

  elsif (!$has_stop) {
    $sth = $self->dbh->do_query('select fdna,foffset from fdna where lower(fref)=lower(?) and foffset>=? order by foffset',
                                $ref,$offset_start);
  }

  else {  # both start and stop defined
    $sth = $self->dbh->do_query('select fdna,foffset from fdna where lower(fref)=lower(?) and foffset>=? and foffset<=? order by foffset',
                                $ref,$offset_start,$offset_stop);
  }

  my $dna = '';
  while (my($frag,$offset) = $sth->fetchrow_array) {
      substr($frag,0,$start-$offset) = '' if $has_start && $start > $offset;
      $dna .= $frag;
  }
  substr($dna,$stop-$start+1) = '' if $has_stop && $stop-$start+1 < length($dna);
  if ($reversed) {
    $dna = reverse $dna;
    $dna =~ tr/gatcGATC/ctagCTAG/;
  }

  $sth->finish;
  $dna;
}


sub refseq_query {
  my $self = shift;
  my ($refseq,$refclass) = @_;
  my $query = "lower(fdata.fref)=lower(?)";
  return wantarray ? ($query,$refseq) : $self->dbh->dbi_quote($query,$refseq);
}

sub make_types_where_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count,$typelist) = @_;
  my (@query,@args);
  if (defined($srcseq)) {
    push @query,'lower(fdata.fref)=lower(?)';
    push @args,$srcseq;
    if (defined $start or defined $stop) {
      $start = 1           unless defined $start;
      $stop  = MAX_SEGMENT unless defined $stop;
      my ($q,@a) = $self->overlap_query($start,$stop);
      push @query,"($q)";
      push @args,@a;
    }
  }
  if (defined $typelist && @$typelist) {
    my ($q,@a) = $self->types_query($typelist);
    push @query,($q);
    push @args,@a;
  }
  my $query = @query ? join(' AND ',@query) : '1=1';
  return wantarray ? ($query,@args) : $self->dbh->dbi_quote($query,@args);
}

sub get_feature_id {
  my $self = shift;
  my ($ref,$start,$stop,$typeid,$groupid) = @_;
  my $s = $self->{load_stuff};
  unless ($s->{get_feature_id}) {
    my $dbh = $self->features_db;
    $s->{get_feature_id} =
      $dbh->prepare_delayed('SELECT fid FROM fdata WHERE lower(fref)=lower(?) AND fstart=? AND fstop=? AND ftypeid=? AND gid=?');
  }
  my $sth = $s->{get_feature_id} or return;
  $sth->execute($ref,$start,$stop,$typeid,$groupid) or return;
  my ($fid) = $sth->fetchrow_array;
  return $fid;
}

sub _delete {
  my $self = shift;
  my $delete_spec = shift;
  my $ranges      = $delete_spec->{segments} || [];
  my $types       = $delete_spec->{types}    || [];
  my $force       = $delete_spec->{force};
  my $range_type  = $delete_spec->{range_type};
  my $dbh         = $self->features_db;

  my $query = 'delete from fdata';
  my @where;

  my @range_part;
  for my $segment (@$ranges) {
    my $ref   = $dbh->quote($segment->abs_ref);
    my $start = $segment->abs_start;
    my $stop  = $segment->abs_stop;
    my $range =  $range_type eq 'overlaps'     ? $self->overlap_query($start,$stop)
               : $range_type eq 'contains'     ? $self->contains_query($start,$stop)
               : $range_type eq 'contained_in' ? $self->contained_in_query($start,$stop)
               : $self->throw("Invalid range type '$range_type'");
    push @range_part,"(lower(fref)=lower($ref) AND $range)";
  }
  push @where,'('. join(' OR ',@range_part).')' if @range_part;

  # get all the types
  if (@$types) {
    my $types_where = $self->types_query($types);
    my $types_query = "select ftypeid from ftype where $types_where";
    my $result      = $dbh->selectall_arrayref($types_query);
    my @typeids     = map {$_->[0]} @$result;
    my $typelist    = join ',',map{$dbh->quote($_)} @typeids;
    $typelist ||= "0"; # don't cause DBI to die with invalid SQL when
                       # unknown feature types were requested.
    push @where,"(ftypeid in ($typelist))";
  }
  $self->throw("This operation would delete all feature data and -force not specified")
    unless @where || $force;
  $query .= " where ".join(' and ',@where) if @where;
  warn "$query\n" if $self->debug;
  my $result = $dbh->do($query);
  defined $result or $self->throw($dbh->errstr);
  $result;
}

sub make_abscoord_query {
  my $self = shift;
  my ($name,$class,$refseq) = @_;
  #my $query = GETSEQCOORDS;
  my $query = $self->getseqcoords_query();
  my $getforcedseqcoords = $self->getforcedseqcoords_query() ;
  if ($name =~ /\*/) {
    $name =~ s/%/\\%/g;
    $name =~ s/_/\\_/g;
    $name =~ tr/*/%/;
    $query =~ s/gname\) = lower/gname) LIKE lower/;
  }
  defined $refseq 
    ? $self->dbh->do_query($getforcedseqcoords,$name,$class,$refseq)
    : $self->dbh->do_query($query,$name,$class);
}

sub make_aliasabscoord_query {
  my $self = shift;
  my ($name,$class) = @_;
  #my $query = GETALIASCOORDS;
  my $query = $self->getaliascoords_query();
  if ($name =~ /\*/) {
    $name =~ s/%/\\%/g;
    $name =~ s/_/\\_/g;
    $name =~ tr/*/%/;
    $query =~ s/gname\) = lower/gname) LIKE lower/;
  }
  $self->dbh->do_query($query,$name,$class);
}


1;
