package Bio::DB::GFF::Adaptor::dbi::mysqlopt;

# an optimized mysql adaptor
# For speed, this module:
# 1) uses Bio::DB::Fasta to retrieve the sequence
# 2) uses a binning scheme to make short range queries very much faster
# 3) optionally returns Acedb objects, for compatibility with WormBase

use strict;
use Bio::DB::GFF::Adaptor::dbi::mysql;
use Bio::DB::Fasta;
use Bio::DB::GFF::Util::Binning;
use Bio::DB::GFF::Util::Rearrange;

use vars qw($VERSION @ISA);
@ISA = qw(Bio::DB::GFF::Adaptor::dbi::mysql);
$VERSION = 0.30;

# this is the largest that any reference sequence can be (100 megabases)
use constant MAX_BIN    => 100_000_000;

# this is the smallest bin (1 K)
use constant MIN_BIN    => 1000;

# size of range over which it is faster to force mysql to use the range for indexing
use constant STRAIGHT_JOIN_LIMIT => 200_000;

sub new {
  my $class = shift;
  my ($dna_db,$acedb,
      $minbin,$maxbin,
      $join_limit,$other) =  rearrange([
					[qw(DNADB DNA FASTA FASTA_DIR)],
					'ACEDB',
					'MINBIN',
					'MAXBIN',
					[qw(STRAIGHT_JOIN_LIMIT JOIN_LIMIT)],
				       ],@_);
  my $self = $class->SUPER::new($other);

  if ($dna_db) {
    if (!ref($dna_db)) {
      my $fasta_dir = $dna_db;
      $dna_db = Bio::DB::Fasta->new($fasta_dir);
      $dna_db or $class->throw("new(): Failed to create new Bio::DB::Fasta from files in $fasta_dir");
    } else {
      $dna_db->isa('Bio::DB::Fasta') or $class->throw("new(): $dna_db is not a Bio::DB::Fasta object");
    }
    $self->{dna_db} = $dna_db;
  }

  if ($acedb) {
    $acedb->isa('Ace') or $class->throw("$acedb is not an acedb accessor object");
    $self->{acedb} = $acedb;
  }

  $self->{minbin} = defined($minbin) ? $minbin : MIN_BIN;
  $self->{maxbin} = defined($maxbin) ? $maxbin : MAX_BIN;
  $self->{straight_join_limit} =  defined($join_limit) ? $join_limit : STRAIGHT_JOIN_LIMIT;

  $self;
}

sub dna_db      { shift->{dna_db}      }
sub acedb       { shift->{acedb}       }


# given sequence name, and optional (start,stop) give raw dna
sub get_dna {
  my $self = shift;
  my ($name,$start,$stop,$class) = @_;
  my $dna_db = $self->dna_db or return $self->SUPER::get_dna(@_);
  # in actuality, the class is simply ignored by Bio::DB::Fasta
  $dna_db->seq($name,$start,$stop,$class);
}

sub do_straight_join {
  my $self = shift;
  my($refseq,$start,$stop,$types) = @_;

  # Might try turning on and off straight join based on the number of types
  # specified, but this turns out to be very difficult indeed!

  # if a list of types has been specified, then it is almost always faster
  # to let the query optimizer figure it out.
  # the exception is when a type of "similarity" has been specified, in which
  # case the range query is better.  (yes, definitely a hack)
  # return 0 if defined($types) and @$types > 0 
  # and !grep {$_->[0] =~ /similarity/ } @$types;

  # if no types are specified then it is faster to do a range search, up to a point.
  return $refseq && abs($stop-$start) < $self->{straight_join_limit};
}


# find features that overlap a given range
sub overlap_query {
  my $self = shift;
  my ($start,$stop) = @_;

  my ($bq,@bargs)   = $self->bin_query($start,$stop);
  my ($iq,@iargs) = $self->SUPER::overlap_query($start,$stop);
  my $query = "($bq)\n\tAND $iq";
  my @args  = (@bargs,@iargs);

  return wantarray ? ($query,@args) : $self->dbi_quote($query,@args);
}

# find features that are completely contained within a range
sub range_query {
  my $self = shift;
  my ($start,$stop) = @_;
  my ($bq,@bargs)   = $self->bin_query($start,$stop);
  my ($iq,@iargs) = $self->SUPER::range_query($start,$stop);
  my $query = "($bq)\n\tAND $iq";
  my @args  = (@bargs,@iargs);
  return wantarray ? ($query,@args) : $self->dbi_quote($query,@args);
}

sub make_object {
  my $self = shift;
  my ($class,$name,$start,$stop) = @_;

  if (my $db = $self->acedb) {
    my $class = $db->class;

    # for Notes we just return a text, no database associated
    return $class->new(Text=>$name) if $class eq 'Note';

    # for homols, we create the indicated Protein or Sequence object
    # then generate a bogus Homology object (for future compatability??)
    return Ace::Sequence::Homol->new($class,$name,$db,$start,$stop) if defined $start;

    # General case:
    my $obj = $class->new($class=>$name,$self->db);

    return $obj if defined $obj;

    # Last resort, return a Text
    return $class->new(Text=>$name);
  }

  return $self->SUPER::make_object($class,$name,$start,$stop);
}

sub bin_query {
  my $self = shift;
  my ($start,$stop) = @_;
  my ($query,@args);

  $start = 0               unless defined($start);
  $stop  = $self->{maxbin} unless defined($stop);

  my @bins;
  my $tier = $self->{maxbin};
  while ($tier >= $self->{minbin}) {
    push @bins,'fbin between ? and ?';
    push @args,bin_bot($tier,$start),bin_top($tier,$stop);
    $tier /= 10;
  }
  $query = join("\n\t OR ",@bins);
  return wantarray ? ($query,@args)
                   : $self->dbi_quote($query,@args);
}

########################## loading and initialization  #####################

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
INSERT INTO fdata (fref,fstart,fstop,fbin,ftypeid,fscore,
		   fstrand,fphase,gid,ftarget_start,ftarget_stop)
       VALUES(?,?,?,?,?,?,?,?,?,?,?)
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

    my $bin =  bin($start,$stop,$self->{minbin});

    my $result = $insert_data->execute($ref,$start,$stop,$bin,$typeid,
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

sub split_group {
  local $_ = shift;
  my ($tag,$value) = /(\S+)\s+\"([^\"]+)\"/;
  $value =~ s/\\t/\t/g;
  $value =~ s/\\r/\r/g;
  my $class = $tag;
  if ($value =~ /^(\S+):(\S.*)$/) {
    $class = $1;
    $value = $2;
  }
  return ($class,$value) unless $tag eq 'Target';
  $value =~ s/^\"//;
  $value =~ s/\"$//;
  my($start,$end) = /(\d+) (\d+)/;
  return ($class,$value,$start,$end);
}

sub schema {
  return <<END;
create table fdata (
    fid	         int not null auto_increment,
    fref         varchar(20)    not null,
    fstart       int unsigned   not null,
    fstop        int unsigned   not null,
    fbin         decimal(20,6)  not null,
    ftypeid      int not null,
    fscore        float,
    fstrand       enum('+','-'),
    fphase        enum('0','1','2'),
    gid          int not null,
    ftarget_start int unsigned,
    ftarget_stop  int unsigned,
    primary key(fid),
    index(fref,fbin,fstart,fstop,ftypeid),
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
;
}



1;

__END__

