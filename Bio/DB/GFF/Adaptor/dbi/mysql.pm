package Bio::DB::GFF::Adaptor::dbi::mysql;

# a simple mysql adaptor
require 5.6.0;
use strict;
use base 'Bio::DB::GFF::Adaptor::dbi';

use constant MAX_SEGMENT => 100_000_000;  # the largest a segment can get

use constant GETSEQCOORDS =><<END;
SELECT fref,
       min(fstart),
       max(fstop),
       fstrand
  FROM fdata,fgroup
  WHERE fgroup.gclass=?
    AND fgroup.gname=?
    AND fgroup.gid=fdata.gid
    GROUP BY fref,fstrand
END

sub make_dna_query {
  my $self = shift;
  $self->do_query('SELECT substring(fdna.fdna,?,?) FROM fdna WHERE fref=?',@_);
}

# given sequence name, return (reference,start,stop,strand)
sub make_abscoord_query {
  my $self = shift;
  my ($class,$name) = @_;
  return wantarray ? (GETSEQCOORDS,$class,$name) 
                   : $self->dbi_quote(GETSEQCOORDS,$class,$name);
}

sub make_features_select_part {
  my $self = shift;
  return <<END;
fstart,fstop,fmethod,fsource,fscore,fstrand,fphase,gclass,gname,ftarget_start,ftarget_stop
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
  my ($refseq,$refclass) = @_;
  my $query = "fdata.fref = ?\n";
  return wantarray ? ($query,$refseq) : $self->dbi_quote($query,$refseq);
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
    my $meth_query = $self->string_match('fmethod',$method) if defined $method;
    my $src_query  = $self->string_match('fsource',$source) if defined $source;
    my @pair;
    if (defined $method) {
      push @pair,$self->string_match('fmethod',$method);
      push @args,$method;
    }
    if (defined $source) {
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
  my ($refseq,$start,$stop,$want_count) = @_;
  my $query = $want_count ? 'ftype.fmethod,ftype.fsource,count(fdata.ftypeid)'
                          : 'fmethod,fsource';
  return $query;
}

sub make_types_from_part {
  my $self = shift;
  my ($refseq,$start,$stop,$want_count) = @_;
  my $query = defined($refseq) || $want_count ? 'fdata,ftype' : 'ftype';
  return $query;
}

sub make_types_join_part {
  my $self = shift;
  my ($refseq,$start,$stop,$want_count) = @_;
  my $query = defined($refseq) || $want_count ? 'fdata.ftypeid=ftype.ftypeid'
                                              : '';
  return $query || 1;
}

sub make_types_where_part {
  my $self = shift;
  my ($refseq,$start,$stop,$want_count) = @_;
  my ($query,@args);
  if (defined($refseq)) {
    $query .= 'fdata.fref=?';
    push @args,$refseq;
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
  my ($refseq,$start,$stop,$want_count) = @_;
  return unless $refseq or $want_count;
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

sub load_gff {
  my $self      = shift;

  my $dbh = $self->features_db;
  local $dbh->{PrintError} = 0;

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

  # local caches of type and group ids
  my (%types,%groups,$counter);

  while (<>) {
    my ($ref,$source,$method,$start,$stop,$score,$strand,$phase,$group) = split "\t";
    next if /^\#/;

    my ($group_class,$group_name,@rest) = split_group($group);

    # truncate the group name
    $group_name = substr($group_name,0,100) if length $group_name > 100;

    my $target_start = $rest[0];
    my $target_stop  = $rest[1];
    $group_class = 'Sequence' if $group_class eq 'Target';

    # get the type ID
    my $key = "\L$method$;$source\E";
    unless ($types{$key}) {

      if ( (my $result = $lookup_type->execute($method,$source)) > 0) {
	$types{$key} = ($lookup_type->fetchrow_array)[0];
      } else {
	$insert_type->execute($method,$source)
	  && ($types{$key} = $dbh->{mysql_insertid});
      }
    }

    my $typeid = $types{$key};
    unless ($typeid) {
      warn "No typeid for $method:$source; ",$dbh->errstr," Record skipped.\n";
      next;
    }

    # and the group ID
    $key = "\L$group_class$;$group_name\E";
    unless ($groups{$key}) {

      if ((my $result = $lookup_group->execute($group_class,$group_name)) > 0) {
	$groups{$key} = ($lookup_group->fetchrow_array)[0];
      } else {
	$insert_group->execute($group_class,$group_name)
	  && ($groups{$key} = $dbh->{mysql_insertid});
      }
    }

    my $groupid = $groups{$key};
    unless ($groupid) {
      warn "No groupid for $group_class:$group_name; ",$dbh->errstr," Record skipped.\n";
      next;
    }

    my $result = $insert_data->execute($ref,$start,$stop,$typeid,
				       $score,$strand,$phase,$groupid,
				       $target_start,$target_stop);
    unless ($result) {
      warn $dbh->errstr,"\n";
      next;
    }

    next unless $result;
    if ( (++$counter % 1000) == 0) {
      print STDERR "$counter records loaded...";
      print STDERR -t STDOUT && !$ENV{EMACS} ? "\r" : "\n";
    }
  }

  $_->finish foreach ($lookup_type,$insert_type,$lookup_group,$insert_group,$insert_data);

#  $dbh->do('UNLOCK TABLES');

  return $counter;
}

1;
