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
fgroup.gid = fdata.gid AND ftype.ftypeid = fdata.ftypeid
END
}

# IMPORTANT NOTE: THE MYSQL SCHEMA IGNORES THE SEQUENCE CLASS
sub refseq_query {
  my $self = shift;
  my ($srcseq,$refclass) = @_;
  my $query = "fdata.fref = ?\n";
  return wantarray ? ($query,$srcseq) : $self->dbi_quote($query,$srcseq);
}


# find features that overlap a given range
sub overlap_query {
  my $self = shift;
  my ($start,$stop) = @_;

  my $query    = qq(fdata.fstop>=? AND fdata.fstart<=?);
  return wantarray ? ($query,$start,$stop) : $self->dbi_quote($query,$start,$stop);
}

# find features that are completely contained within a range
sub range_query {
  my $self = shift;
  my ($start,$stop) = @_;
  my $query    = qq(fdata.fstart>=? AND fdata.fstop<=?);
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
    my $meth_query = $self->string_match('fmethod',$method) if defined $method && length $method;
    my $src_query  = $self->string_match('fsource',$source) if defined $source && length $source;
    my @pair;
    if (defined $method && length $method) {
      push @pair,$self->string_match('fmethod',$method);
      push @args,$method;
    }
    if (defined $source && length $source) {
      push @pair,$self->string_match('fsource',$source);
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
  my ($srcseq,$start,$stop,$want_count) = @_;
  my ($query,@args);
  if (defined($srcseq)) {
    $query .= 'fdata.fref=?';
    push @args,$srcseq;
    if (defined $start or defined $stop) {
      $start = 1           unless defined $start;
      $stop  = MAX_SEGMENT unless defined $stop;
      my ($q,@a) = $self->overlap_query($start,$stop);
      $query .= " AND ($q)";
      push @args,@a;
    }
  } else {
    $query = '1';
  }
  return wantarray ? ($query,@args) : $self->dbi_quote($query,@args);
}

sub make_types_group_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  return unless $srcseq or $want_count;
  return 'ftype.ftypeid';
}

################################ loading and initialization ##################################
# return list of tables that "belong" to us.
sub tables {
  qw(fdata fgroup ftype fdna)
}

sub schema {
  return <<END;
create table fdata (
    fid	         int not null auto_increment,
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
    gid	    int not null auto_increment,
    gclass  varchar(20),
    gname   varchar(100),
    primary key(gid),
    unique(gclass,gname)
);

create table ftype (
    ftypeid      int not null auto_increment,
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
END
}

sub setup_load {
  my $self      = shift;

  my $dbh = $self->features_db;

  # for the paranoid....
  #  $dbh->do("LOCK TABLES fdata WRITE, ftype WRITE, fgroup WRITE");

  my $lookup_type = $dbh->prepare('SELECT ftypeid FROM ftype WHERE fmethod=? AND fsource=?');
  my $insert_type = $dbh->prepare('INSERT INTO ftype (fmethod,fsource) VALUES (?,?)');

  my $lookup_group = $dbh->prepare('SELECT gid FROM fgroup WHERE gclass=? AND gname=?');
  my $insert_group = $dbh->prepare('INSERT INTO fgroup (gclass,gname) VALUES (?,?)');

  my $insert_data  = $dbh->prepare(<<END);
INSERT INTO fdata (fref,fstart,fstop,ftypeid,fscore,
		   fstrand,fphase,gid,ftarget_start,ftarget_stop)
       VALUES(?,?,?,?,?,?,?,?,?,?)
END
;

  $self->{load_stuff}{lookup_type}  = $lookup_type;
  $self->{load_stuff}{insert_type}  = $insert_type;
  $self->{load_stuff}{lookup_group} = $lookup_group;
  $self->{load_stuff}{insert_group} = $insert_group;
  $self->{load_stuff}{insert_data}  = $insert_data;
  $self->{load_stuff}{types}  = {};
  $self->{load_stuff}{groups} = {};
  $self->{load_stuff}{counter} = 0;
}

sub load_gff_line {
  my $self = shift;
  my ($ref,$source,$method,$start,$stop,
      $score,$strand,$phase,
      $group_class,$group_name,$target_start,$target_stop,
      $notes) = @_;

  my $s    = $self->{load_stuff};
  my $dbh  = $self->features_db;
  local $dbh->{PrintError} = 0;

  # get the type ID
  my $key = "\L$method$;$source\E";
  unless ($s->{types}{$key}) {

    if ( (my $result = $s->{lookup_type}->execute($method,$source)) > 0) {
      $s->{types}{$key} = ($s->{lookup_type}->fetchrow_array)[0];
    } else {
      $s->{insert_type}->execute($method,$source)
	&& ($s->{types}{$key} = $dbh->{mysql_insertid});
    }
  }

  my $typeid = $s->{types}{$key};
  unless ($typeid) {
    warn "No typeid for $method:$source; ",$dbh->errstr," Record skipped.\n";
    next;
  }

  # and the group ID
  $key = "\L$group_class$;$group_name\E";
  unless ($s->{groups}{$key}) {

    if ((my $result = $s->{lookup_group}->execute($group_class,$group_name)) > 0) {
      $s->{groups}{$key} = ($s->{lookup_group}->fetchrow_array)[0];
    } else {
      $s->{insert_group}->execute($group_class,$group_name)
	&& ($s->{groups}{$key} = $dbh->{mysql_insertid});
    }
  }

  my $groupid = $s->{groups}{$key};
  unless ($groupid) {
    warn "No groupid for $group_class:$group_name; ",$dbh->errstr," Record skipped.\n";
    next;
  }

  my $result = $s->{insert_data}->execute($ref,$start,$stop,$typeid,
					  $score,$strand,$phase,$groupid,
					  $target_start,$target_stop);
  unless ($result) {
    warn $dbh->errstr,"\n";
    next;
  }

  return unless $result;

  if ( (++$s->{counter} % 1000) == 0) {
    print STDERR "$s->{counter} records loaded...";
    print STDERR -t STDOUT && !$ENV{EMACS} ? "\r" : "\n";
  }

  1;
}

sub finish_load {
  my $self = shift;
#  $dbh->do('UNLOCK TABLES');

  $self->{load_stuff}{$_}->finish
    foreach qw(lookup_type insert_type lookup_group insert_group insert_data);

  my $counter = $self->{load_stuff}{counter};
  delete $self->{load_stuff};
  return $counter;
}

1;
