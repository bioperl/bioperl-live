package Bio::DB::GFF::Adaptor::dbi::mysql;

=head1 NAME

Bio::DB::GFF::Adaptor::dbi::mysql -- Database adaptor for a specific mysql schema

=head1 SYNOPSIS

See L<Bio::DB::GFF>

=cut

# a simple mysql adaptor
use strict;
use Bio::DB::GFF::Adaptor::dbi;
use Bio::DB::GFF::Util::Rearrange; # for rearrange()
use vars qw($VERSION @ISA $IN_ITERATOR);
@ISA = qw(Bio::DB::GFF::Adaptor::dbi);
$VERSION = '0.26';

$IN_ITERATOR = 0;

use constant MAX_SEGMENT => 100_000_000;  # the largest a segment can get

use constant GETSEQCOORDS =><<END;
SELECT fref,
       IF(ISNULL(gclass),'Sequence',gclass),
       min(fstart),
       max(fstop),
       fstrand
  FROM fdata,fgroup
  WHERE fgroup.gname=?
    AND fgroup.gclass=?
    AND fgroup.gid=fdata.gid
    GROUP BY fref,fstrand
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

=head1 DESCRIPTION

This adaptor implements a specific mysql database schema that is
compatible with Bio::DB::GFF.  It inherits from
Bio::DB::GFF::Adaptor::dbi, which itself inherits from Bio::DB::GFF.

The schema uses four tables:

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

This table holds the raw DNA of the reference sequences.  It has two
columns:

    fref          reference sequence name (string)
    fdna          the DNA sequence (longblob)

Note that the GFF module will not help you load this table.  You must
load it yourself.

=item fnote

This table holds "notes", which some groups have used the GFF group
field to represent.  Notes are created by creating a group class of
"Note" and a group value containing the text of the note, as shown in
this example:

 CHR_I assembly_tag Finished     2032 2036 . + . Note "Right: cTel33B"
 CHR_I assembly_tag Polymorphism 668  668  . + . Note "A->C in cTel33B"

The columns of this table are:

    fid      feature ID (integer)
    fnote    text of the note (text)

The fdata.fid column joins with fnote.fid.

=back

=head2 The fast_queries() Method

The Perl mysql driver offers a statement handler attribute called
i<mysql_use_result>.  If this is set to true, then the server will
keep the result of the query in memory rather than downloading it to
the client in one go.  This is a big win if you are iterating through
a very large result set.

The drawback is that while you are iterating through the set, then you
cannot trigger any other SQL queries.  Doing so will result in a
"command out of synch" error from the mysql server. For this reason,
fast queries are turned off by default.

=over 4

=item $flag = $db->fast_queries([$flag])

Get or set the fast_queries flag.

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
  if (!ref($dsn) && $dsn !~ /^dbi:mysql/) {
      $dsn = "dbi:mysql:$dsn";
  }
  $class->SUPER::new(-dsn=>$dsn,%$other);
}

=head2 make_dna_query

 Title   : make_dna_query
 Usage   : $sth = $db->make_dna_query($name,$class,$start,$stop);
 Function: create query that returns raw DNA sequence from database
 Returns : a DBI statement handle
 Args    : reference sequence name and class, and a range
 Status  : protected

The statement handler should return rows containing just one field,
the extracted DNA string.

=cut

sub make_dna_query {
  my $self = shift;
  my ($name,$start,$stop,$class) = @_;
  # simply ignore class for now
  $self->do_query('SELECT substring(fdna.fdna,?,?) FROM fdna WHERE fref=?',$name,$start,$stop);
}

=head2 make_abscoord_query

 Title   : make_abscoord_query
 Usage   : $sth = $db->make_abscoord_query($name,$class);
 Function: create query that finds the reference sequence coordinates given a landmark & classa
 Returns : a DBI statement handle
 Args    : name and class of landmark
 Status  : protected

The statement handler should return rows containing five fields:

  1. reference sequence name
  2. reference sequence class
  3. start position
  4. stop position
  5. strand ("+" or "-")

This query always returns "Sequence" as the class of the reference
sequence.

=cut

# given sequence name, return (reference,start,stop,strand)
sub make_abscoord_query {
  my $self = shift;
  my ($name,$class,$refseq) = @_;
  defined $refseq ? $self->do_query(GETFORCEDSEQCOORDS,$name,$class,$refseq) 
    : $self->do_query(GETSEQCOORDS,$name,$class);
}

=head2 make_features_byname_where_part

 Title   : make_features_byname_where_part
 Usage   : $db->make_features_byname_where_part
 Function: create the SQL fragment needed to select a feature by its group name & class
 Returns : a SQL fragment and bind arguments
 Args    : see below
 Status  : Protected

=cut

sub make_features_byname_where_part {
  my $self = shift;
  my ($class,$name) = @_;
  return ("fgroup.gclass=? AND fgroup.gname=?",$class,$name);
}

=head2 make_features_select_part

 Title   : make_features_select_part
 Usage   : $string = $db->make_features_select_part()
 Function: make select part of the features query
 Returns : a string
 Args    : none
 Status  : protected

This method creates the part of the features query that immediately
follows the SELECT keyword.

=cut

sub make_features_select_part {
  my $self = shift;
  return <<END;
fref,fstart,fstop,fsource,fmethod,fscore,fstrand,fphase,gclass,gname,ftarget_start,ftarget_stop,fid
END
}

=head2 make_features_from_part

 Title   : make_features_from_part
 Usage   : $string = $db->make_features_from_part()
 Function: make from part of the features query
 Returns : a string
 Args    : none
 Status  : protected

This method creates the part of the features query that immediately
follows the FROM keyword.

=cut

sub make_features_from_part {
  my $self = shift;
  return "fdata,ftype,fgroup\n";
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
  return <<END;
  fgroup.gid = fdata.gid 
  AND ftype.ftypeid = fdata.ftypeid
END
}

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
  return "fdata.gid";
}

=head2 refseq_query

 Title   : refseq_query
 Usage   : ($query,@args) = $db->refseq_query($name,$class)
 Function: create SQL fragment that selects the desired reference sequence
 Returns : a list containing the query and bind arguments
 Args    : reference sequence name and class
 Status  : protected

This method is called by make_features_byrange_where_part() to construct the
part of the select WHERE section that selects a particular reference
sequence.  It returns a mult-element list in which the first element
is the SQL fragment and subsequent elements are bind values.

The current schema does not distinguish among different classes of
reference sequence.

=cut

# IMPORTANT NOTE: THE MYSQL SCHEMA IGNORES THE SEQUENCE CLASS
# THIS SHOULD BE FIXED
sub refseq_query {
  my $self = shift;
  my ($refseq,$refclass) = @_;
  my $query = "fdata.fref=?";
  return wantarray ? ($query,$refseq) : $self->dbi_quote($query,$refseq);
}

=head2 notes

 Title   : notes
 Usage   : @notes = $db->notes($id)
 Function: return the list of notes corresponding to a feature ID
 Returns : a list of strings
 Args    : feature ID
 Status  : protected

This method is called by Bio::DB::GFF->notes() to retrieve the notes
corresponding to the internal feature ID.

=cut

sub notes {
  my $self = shift;
  my $id = shift;
  my $sth = $self->do_query('SELECT fnote FROM fnote WHERE fid=?',$id) or return;

  my @result;
  while (my($note) = $sth->fetchrow_array) {
    push @result,$note;
  }
  return @result;
  $sth->finish;
}

=head2 overlap_query

 Title   : overlap_query
 Usage   : ($query,@args) = $db->overlap_query($start,$stop)
 Function: create SQL fragment that selects the desired features by range
 Returns : a list containing the query and bind arguments
 Args    : the start and stop of a range, inclusive
 Status  : protected

This method is called by make_features_byrange_where_part() to construct the
part of the select WHERE section that selects a set of features that
overlap a range. It returns a multi-element list in which the first
element is the SQL fragment and subsequent elements are bind values.

=cut

# find features that overlap a given range
sub overlap_query {
  my $self = shift;
  my ($start,$stop) = @_;

  my $query    = qq(fdata.fstop>=? AND fdata.fstart<=?);
  return wantarray ? ($query,$start,$stop) : $self->dbi_quote($query,$start,$stop);
}

=head2 contains_query

 Title   : contains_query
 Usage   : ($query,@args) = $db->contains_query($start,$stop)
 Function: create SQL fragment that selects the desired features by range
 Returns : a list containing the query and bind arguments
 Args    : the start and stop of a range, inclusive
 Status  : protected

This method is called by make_features_byrange_where_part() to construct the
part of the select WHERE section that selects a set of features
entirely enclosed by a range. It returns a multi-element list in which
the first element is the SQL fragment and subsequent elements are bind
values.

=cut

# find features that are completely contained within a range
sub contains_query {
  my $self = shift;
  my ($start,$stop) = @_;
  my $query    = qq(fdata.fstart>=? AND fdata.fstop<=?);
  return wantarray ? ($query,$start,$stop) : $self->dbi_quote($query,$start,$stop);
}

=head2 contained_in_query

 Title   : contained_in_query
 Usage   : ($query,@args) = $db->contained_in_query($start,$stop)
 Function: create SQL fragment that selects the desired features by range
 Returns : a list containing the query and bind arguments
 Args    : the start and stop of a range, inclusive
 Status  : protected

This method is called by make_features_byrange_where_part() to construct the
part of the select WHERE section that selects a set of features
entirely enclosed by a range. It returns a multi-element list in which
the first element is the SQL fragment and subsequent elements are bind
values

=cut

# find features that are completely contained within a range
sub contained_in_query {
  my $self = shift;
  my ($start,$stop) = @_;
  my $query    = qq(fdata.fstart<=? AND fdata.fstop>=?);
  return wantarray ? ($query,$start,$stop) : $self->dbi_quote($query,$start,$stop);
}

=head2 types_query

 Title   : types_query
 Usage   : ($query,@args) = $db->types_query($types)
 Function: create SQL fragment that selects the desired features by type
 Returns : a list containing the query and bind arguments
 Args    : an array reference containing the types
 Status  : protected

This method is called by make_features_byrange_where_part() to construct the
part of the select WHERE section that selects a set of features based
on their type. It returns a multi-element list in which the first
element is the SQL fragment and subsequent elements are bind values.
The argument is an array reference containing zero or more
[$method,$source] pairs.

=cut

# generate the fragment of SQL responsible for searching for
# features with particular types and methods
sub types_query {
  my $self = shift;
  my $types = shift;

  my @method_queries;
  my @args;
  for my $type (@$types) {
    my ($method,$source) = @$type;
    my $meth_query = $self->exact_match('fmethod',$method) if defined $method && length $method;
    my $src_query  = $self->exact_match('fsource',$source) if defined $source && length $source;
    my @pair;
    if (defined $method && length $method) {
      push @pair,$self->exact_match('fmethod',$method);
      push @args,$method;
    }
    if (defined $source && length $source) {
      push @pair,$self->exact_match('fsource',$source);
      push @args,$source;
    }
    push @method_queries,"(" . join(' AND ',@pair) .")" if @pair;
  }
  my $query = " (".join(' OR ',@method_queries).")\n" if @method_queries;
  return wantarray ? ($query,@args) : $self->dbi_quote($query,@args);
}

=head2 make_types_select_part

 Title   : make_types_select_part
 Usage   : ($string,@args) = $db->make_types_select_part(@args)
 Function: create the select portion of the SQL for fetching features type list
 Returns : query string and bind arguments
 Args    : see below
 Status  : protected

This method is called by get_types() to generate the query fragment
and bind arguments for the SELECT part of the query that retrieves
lists of feature types.  The four positional arguments are as follows:

 $refseq      reference sequence name
 $start       start of region
 $stop        end of region
 $want_count  true to return the count of this feature type

If $want_count is false, the SQL fragment returned must produce a list
of feature types in the format (method, source).

If $want_count is true, the returned fragment must produce a list of
feature types in the format (method, source, count).

=cut

#------------------------- support for the types() query ------------------------
sub make_types_select_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  my $query = $want_count ? 'ftype.fmethod,ftype.fsource,count(fdata.ftypeid)'
                          : 'fmethod,fsource';
  return $query;
}

=head2 make_types_from_part

 Title   : make_types_from_part
 Usage   : ($string,@args) = $db->make_types_from_part(@args)
 Function: create the FROM portion of the SQL for fetching features type lists
 Returns : query string and bind arguments
 Args    : see below
 Status  : protected

This method is called by get_types() to generate the query fragment
and bind arguments for the FROM part of the query that retrieves lists
of feature types.  The four positional arguments are as follows:

 $refseq      reference sequence name
 $start       start of region
 $stop        end of region
 $want_count  true to return the count of this feature type

If $want_count is false, the SQL fragment returned must produce a list
of feature types in the format (method, source).

If $want_count is true, the returned fragment must produce a list of
feature types in the format (method, source, count).

=cut

sub make_types_from_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  my $query = defined($srcseq) || $want_count ? 'fdata,ftype' : 'ftype';
  return $query;
}

=head2 make_types_join_part

 Title   : make_types_join_part
 Usage   : ($string,@args) = $db->make_types_join_part(@args)
 Function: create the JOIN portion of the SQL for fetching features type lists
 Returns : query string and bind arguments
 Args    : see below
 Status  : protected

This method is called by get_types() to generate the query fragment
and bind arguments for the JOIN part of the query that retrieves lists
of feature types.  The four positional arguments are as follows:

 $refseq      reference sequence name
 $start       start of region
 $stop        end of region
 $want_count  true to return the count of this feature type

=cut

sub make_types_join_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  my $query = defined($srcseq) || $want_count ? 'fdata.ftypeid=ftype.ftypeid'
                                              : '';
  return $query || 1;
}

=head2 make_types_where_part

 Title   : make_types_where_part
 Usage   : ($string,@args) = $db->make_types_where_part(@args)
 Function: create the WHERE portion of the SQL for fetching features type lists
 Returns : query string and bind arguments
 Args    : see below
 Status  : protected

This method is called by get_types() to generate the query fragment
and bind arguments for the WHERE part of the query that retrieves
lists of feature types.  The four positional arguments are as follows:

 $refseq      reference sequence name
 $start       start of region
 $stop        end of region
 $want_count  true to return the count of this feature type

=cut

sub make_types_where_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count,$typelist) = @_;
  my (@query,@args);
  if (defined($srcseq)) {
    push @query,'fdata.fref=?';
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
  my $query = @query ? join(' AND ',@query) : '1';
  return wantarray ? ($query,@args) : $self->dbi_quote($query,@args);
}

=head2 make_types_group_part

 Title   : make_types_group_part
 Usage   : ($string,@args) = $db->make_types_group_part(@args)
 Function: create the GROUP BY portion of the SQL for fetching features type lists
 Returns : query string and bind arguments
 Args    : see below
 Status  : protected

This method is called by get_types() to generate the query fragment
and bind arguments for the GROUP BY part of the query that retrieves
lists of feature types.  The four positional arguments are as follows:

 $refseq      reference sequence name
 $start       start of region
 $stop        end of region
 $want_count  true to return the count of this feature type

=cut

sub make_types_group_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  return unless $srcseq or $want_count;
  return 'ftype.ftypeid';
}

=head2 do_query

 Title   : do_query
 Usage   : $sth = $db->do_query($query,@args)
 Function: perform a DBI query
 Returns : a statement handler
 Args    : query string and list of bind arguments
 Status  : Public

This method performs a DBI prepare() and execute(), returning a
statement handle.  You will typically call fetch() of fetchrow_array()
on the statement handle.  The parsed statement handle is cached for
later use.

This method differs from that used by the Bio::DB::GFF::Adaptor::dbi
parent in that it sets the dbi::mysql b<mysql_use_result> attribute,
which is much more memory efficient.

=cut

sub do_query {
  my $self = shift;
  my ($query,@args) = @_;
  warn $self->dbi_quote($query,@args),"\n" if $self->debug;

  my $fast_queries = $self->fast_queries;
  if (!defined $fast_queries) {  # user hasn't specified
    $fast_queries = $IN_ITERATOR ? 0 : 1;
  }
  warn "fast queries turned ",($fast_queries ? 'on' : 'off'),"\n" if $self->debug;

  my $sth = $self->{sth}{$query} ||= $self->features_db->prepare($query,
								 {mysql_use_result=>$fast_queries})
    || $self->throw("Couldn't prepare query $query:\n ".DBI->errstr."\n");
  $sth->execute(@args)
    || $self->throw("Couldn't execute query $query:\n ".DBI->errstr."\n");
  $sth;
}

sub get_features_iterator {
  my $self = shift;
  local $IN_ITERATOR = 1;
  $self->SUPER::get_features_iterator(@_);
}


################################ loading and initialization ##################################

=head2 tables

 Title   : tables
 Usage   : @tables = $db->tables
 Function: return list of tables that belong to this module
 Returns : list of tables
 Args    : none
 Status  : protected

This method lists the tables known to the module, namely qw(fdata fref
fgroup ftype fdna fnote).

=cut

# return list of tables that "belong" to us.
sub tables {
  qw(fdata fref fgroup ftype fdna fnote);
}

=head2 schema

 Title   : schema
 Usage   : $schema = $db->schema
 Function: return the CREATE script for the schema
 Returns : a string
 Args    : none
 Status  : protected

This method returns a string containing the various CREATE statements
needed to initialize the database tables.

=cut

sub schema {
  return <<END;
create table fdata (
    fid	         int not null  auto_increment,
    fref         varchar(20)    not null,
    fstart       int unsigned   not null,
    fstop        int unsigned   not null,
    ftypeid      int not null,
    fscore        float,
    fstrand       enum('+','-'),
    fphase        enum('0','1','2'),
    gid          int not null,
    ftarget_start int unsigned,
    ftarget_stop  int unsigned,
    primary key(fid),
    unique index(fref,fstart,fstop,ftypeid,gid),
    index(ftypeid),
    index(gid)
);

create table fgroup (
    gid	    int not null  auto_increment,
    gclass  varchar(20),
    gname   varchar(100),
    primary key(gid),
    unique(gclass,gname)
);

create table ftype (
    ftypeid      int not null   auto_increment,
    fmethod       varchar(30) not null,
    fsource       varchar(30),
    primary key(ftypeid),
    index(fmethod),
    index(fsource),
    unique ftype (fmethod,fsource)
);

create table fdna (
    fref          varchar(20) not null,
    fdna          longblob not null,
    primary key(fref)
);


create table fnote (
    fid      int not null,
    fnote    text,
    index(fid)
);

END
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

  my $lookup_type = $dbh->prepare('SELECT ftypeid FROM ftype WHERE fmethod=? AND fsource=?');
  my $insert_type = $dbh->prepare('INSERT INTO ftype (fmethod,fsource) VALUES (?,?)');

  my $lookup_group = $dbh->prepare('SELECT gid FROM fgroup WHERE gname=? AND gclass=?');
  my $insert_group = $dbh->prepare('INSERT INTO fgroup (gname,gclass) VALUES (?,?)');

  my $insert_note  = $dbh->prepare('INSERT INTO fnote (fid,fnote) VALUES (?,?)');
  my $insert_data  = $dbh->prepare(<<END);
INSERT INTO fdata (fref,fstart,fstop,ftypeid,fscore,
		   fstrand,fphase,gid,ftarget_start,ftarget_stop)
       VALUES(?,?,?,?,?,?,?,?,?,?)
END
;

  $self->{load_stuff}{lookup_ftype}  = $lookup_type;
  $self->{load_stuff}{insert_ftype}  = $insert_type;
  $self->{load_stuff}{lookup_fgroup} = $lookup_group;
  $self->{load_stuff}{insert_fgroup} = $insert_group;
  $self->{load_stuff}{insert_fdata}  = $insert_data;
  $self->{load_stuff}{insert_fnote}  = $insert_note;
  $self->{load_stuff}{types}  = {};
  $self->{load_stuff}{groups} = {};
  $self->{load_stuff}{counter} = 0;
}

=head2 load_gff_line

 Title   : load_gff_line
 Usage   : $db->load_gff_line(@args)
 Function: called to load one parsed line of GFF
 Returns : true if successfully inserted
 Args    : see below
 Status  : protected

This method is called once per line of the GFF and passed a series of
parsed data items.  The items are:

 $ref          reference sequence
 $source       annotation source
 $method       annotation method
 $start        annotation start
 $stop         annotation stop
 $score        annotation score (may be undef)
 $strand       annotation strand (may be undef)
 $phase        annotation phase (may be undef)
 $group_class  class of annotation's group (may be undef)
 $group_name   ID of annotation's group (may be undef)
 $target_start start of target of a similarity hit
 $target_stop  stop of target of a similarity hit
 $notes        array reference of text items to be attached

=cut

sub load_gff_line {
  my $self = shift;
  my $gff = shift;

  my $s    = $self->{load_stuff};
  my $dbh  = $self->features_db;
  local $dbh->{PrintError} = 0;

  defined(my $typeid  = $self->get_table_id('ftype', $gff->{method} => $gff->{source})) or return;
  defined(my $groupid = $self->get_table_id('fgroup',$gff->{gname}  => $gff->{gclass})) or return;

  my $result = $s->{insert_data}->execute($gff->{ref},
					  $gff->{start},$gff->{stop},
					  $typeid,
					  $gff->{score},$gff->{strand},$gff->{phase},
					  $groupid,
					  $gff->{tstart},$gff->{tstop});

  warn $dbh->errstr,"\n" and return unless $result;

  my $fid = $dbh->{mysql_insertid};
  $s->{insert_fnote}->execute($fid,$_) foreach @{$gff->{notes}||[]};

  if ( (++$s->{counter} % 1000) == 0) {
    print STDERR "$s->{counter} records loaded...";
    print STDERR -t STDOUT && !$ENV{EMACS} ? "\r" : "\n";
  }

  $fid;
}

=head2 finish_load

 Title   : finish_load
 Usage   : $db->finish_load
 Function: called after load_gff_line()
 Returns : number of records loaded
 Args    : none
 Status  : protected

This method performs schema-specific cleanup after loading a set of
GFF records.  It finishes each of the statement handlers prepared by
setup_load().

=cut

sub finish_load {
  my $self = shift;

  my $dbh = $self->features_db or return;
  $dbh->do('UNLOCK TABLES') if $self->lock_on_load;

  foreach (qw(lookup_type insert_type lookup_group insert_group insert_data)) {
    next unless defined $self->{load_stuff}{$_};
    $self->{load_stuff}{$_}->finish;
  }

  my $counter = $self->{load_stuff}{counter};
  delete $self->{load_stuff};
  return $counter;
}

=head2 get_table_id

 Title   : get_table_id
 Usage   : $integer = $db->get_table_id($table,$id1,$id2)
 Function: get the ID of a group or type
 Returns : an integer ID or undef
 Args    : none
 Status  : private

This internal method is called by load_gff_line to look up the integer
ID of an existing feature type or group.  The arguments are the name
of the table, and two string identifiers.  For feature types, the
identifiers are the method and source.  For groups, the identifiers
are group name and class.

This method requires that a statement handler named i<lookup_$table>,
have been created previously by setup_load().  It is here to overcome
deficiencies in mysql's INSERT syntax.

=cut


# get the object ID from a named table
sub get_table_id {
  my $self  = shift;
  my $table = shift;
  my ($id1,$id2) = @_;
  my $s = $self->{load_stuff};
  my $dbh = $self->features_db;

  $id1 ||= '';
  $id2 ||= '';

  unless (defined($s->{$table}{$id1,$id2})) {

    if ( (my $result = $s->{"lookup_$table"}->execute($id1,$id2)) > 0) {
      $s->{$table}{$id1,$id2} = ($s->{"lookup_$table"}->fetchrow_array)[0];
    } else {
      $s->{"insert_$table"}->execute($id1,$id2)
	&& ($s->{$table}{$id1,$id2} = $dbh->{mysql_insertid});
    }
  }

  my $id = $s->{$table}{$id1,$id2};
  unless (defined $id) {
    warn "No $table id for $id1:$id2 ",$dbh->errstr," Record skipped.\n";
    return;
  }
  $id;
}

1;
