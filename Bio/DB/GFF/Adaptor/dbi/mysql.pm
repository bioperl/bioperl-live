=head1 NAME

Bio::DB::GFF::Adaptor::dbi::mysql -- Database adaptor for a specific mysql schema

=head1 SYNOPSIS

See L<Bio::DB::GFF>

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
"Note" and a group value containing the text of the note.

Columns are:

    fid      feature ID (integer)
    fnote    text of the note (text)

The fdata.fid column joins with fnote.fid.

=back



=cut

package Bio::DB::GFF::Adaptor::dbi::mysql;

# a simple mysql adaptor
use strict;
use Bio::DB::GFF::Adaptor::dbi;
use vars qw($VERSION @ISA);
@ISA = qw(Bio::DB::GFF::Adaptor::dbi);
$VERSION = '0.20';

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

sub make_dna_query {
  my $self = shift;
  my ($name,$start,$stop,$class) = @_;
  # simply ignore class for now
  $self->do_query('SELECT substring(fdna.fdna,?,?) FROM fdna WHERE fref=?',$name,$start,$stop);
}

# given sequence name, return (reference,start,stop,strand)
sub make_abscoord_query {
  my $self = shift;
  my ($name,$class) = @_;
  $self->do_query(GETSEQCOORDS,$name,$class);
}

sub make_features_select_part {
  my $self = shift;
  return <<END;
fref,fstart,fstop,fsource,fmethod,fscore,fstrand,fphase,gclass,gname,ftarget_start,ftarget_stop,fid
END
}

sub make_features_from_part {
  my $self = shift;
  return "fdata,ftype,fgroup\n";
}

sub make_features_join_part {
  my $self = shift;
  return <<END;
  fgroup.gid = fdata.gid 
  AND ftype.ftypeid = fdata.ftypeid
END
}

sub make_features_order_by_part {
  my $self = shift;
  return "fdata.gid";
}

# IMPORTANT NOTE: THE MYSQL SCHEMA IGNORES THE SEQUENCE CLASS
# THIS SHOULD BE FIXED
sub refseq_query {
  my $self = shift;
  my ($refseq,$refclass) = @_;
  my $query = "fdata.fref=?";
  return wantarray ? ($query,$refseq) : $self->dbi_quote($query,$refseq);
}

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

# find features that overlap a given range
sub overlap_query {
  my $self = shift;
  my ($start,$stop) = @_;

  my $query    = qq(fdata.fstop>=? AND fdata.fstart<=?);
  return wantarray ? ($query,$start,$stop) : $self->dbi_quote($query,$start,$stop);
}

# find features that are completely contained within a range
sub contains_query {
  my $self = shift;
  my ($start,$stop) = @_;
  my $query    = qq(fdata.fstart>=? AND fdata.fstop<=?);
  return wantarray ? ($query,$start,$stop) : $self->dbi_quote($query,$start,$stop);
}

# find features that are completely contained within a range
sub contained_in_query {
  my $self = shift;
  my ($start,$stop) = @_;
  my $query    = qq(fdata.fstart<=? AND fdata.fstop>=?);
  return wantarray ? ($query,$start,$stop) : $self->dbi_quote($query,$start,$stop);
}

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

#------------------------- support for the types() query ------------------------
sub make_types_select_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  my $query = $want_count ? 'ftype.fmethod,ftype.fsource,count(fdata.ftypeid)'
                          : 'fmethod,fsource';
  return $query;
}

sub make_types_from_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  my $query = defined($srcseq) || $want_count ? 'fdata,ftype' : 'ftype';
  return $query;
}

sub make_types_join_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  my $query = defined($srcseq) || $want_count ? 'fdata.ftypeid=ftype.ftypeid'
                                              : '';
  return $query || 1;
}

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

sub make_types_group_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  return unless $srcseq or $want_count;
  return 'ftype.ftypeid';
}

# override this method in order to set the mysql_use_result attribute, which is an obscure
# but extremely powerful optimization for both performance and memory.
sub do_query {
  my $self = shift;
  my ($query,@args) = @_;
  warn $self->dbi_quote($query,@args),"\n" if $self->debug;
  my $sth = $self->{sth}{$query} ||= $self->features_db->prepare($query,{mysql_use_result=>1})
    || $self->throw("Couldn't prepare query $query:\n ".DBI->errstr."\n");
  $sth->execute(@args)
    || $self->throw("Couldn't execute query $query:\n ".DBI->errstr."\n");
  $sth;
}

################################ loading and initialization ##################################
# return list of tables that "belong" to us.
sub tables {
  qw(fdata fref fgroup ftype fdna fnote);
}

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
    index(fref,fstart,fstop,ftypeid),
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

# get the object ID from a named table
sub get_table_id {
  my $self  = shift;
  my $table = shift;
  my ($id1,$id2) = @_;
  my $s = $self->{load_stuff};
  my $dbh = $self->features_db;

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

sub finish_load {
  my $self = shift;

  my $dbh = $self->features_db or return;
  $dbh->do('UNLOCK TABLES') if $self->lock_on_load;

  $self->{load_stuff}{$_}->finish
    foreach qw(lookup_type insert_type lookup_group insert_group insert_data);

  my $counter = $self->{load_stuff}{counter};
  delete $self->{load_stuff};
  return $counter;
}

1;
