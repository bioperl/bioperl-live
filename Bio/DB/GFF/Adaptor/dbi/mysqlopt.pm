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
  my ($name,$class,$start,$stop) = @_;
  my $dna_db = $self->dna_db or return $self->SUPER::get_dna(@_);
  # in actuality, the class is simply ignored by Bio::DB::Fasta
  $dna_db->seq($name,$start,$stop,$class);
}

sub do_straight_join {
  my $self = shift;
  my($srcseq,$start,$stop,$types) = @_;

  # Might try turning on and off straight join based on the number of types
  # specified, but this turns out to be very difficult indeed!

  # if a list of types has been specified, then it is almost always faster
  # to let the query optimizer figure it out.
  # the exception is when a type of "similarity" has been specified, in which
  # case the range query is better.  (yes, definitely a hack)
  # return 0 if defined($types) and @$types > 0 
  # and !grep {$_->[0] =~ /similarity/ } @$types;

  # if no types are specified then it is faster to do a range search, up to a point.
  return $srcseq && defined($start) && defined($stop) && abs($stop-$start) < $self->{straight_join_limit};
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
sub contains_query {
  my $self = shift;
  my ($start,$stop) = @_;
  my ($bq,@bargs)   = $self->bin_query($start,$stop);
  my ($iq,@iargs) = $self->SUPER::contains_query($start,$stop);
  my $query = "($bq)\n\tAND $iq";
  my @args  = (@bargs,@iargs);
  return wantarray ? ($query,@args) : $self->dbi_quote($query,@args);
}

# find features that are completely contained within a range
sub contained_in_query {
  my $self = shift;
  my ($start,$stop) = @_;
  my ($bq,@bargs)   = $self->bin_query($start,$stop);
  my ($iq,@iargs) = $self->SUPER::contained_in_query($start,$stop);
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
    my $obj = $class->new($class=>$name,$self->acedb);

    return $obj if defined $obj;

    # Last resort, return a Text
    return $class->new(Text=>$name);
  }

  return $self->SUPER::make_object($name,$class,$start,$stop);
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

# exactly like mysql.pm setup_load except that it inserts the bin field
sub setup_load {
  my $self      = shift;

  $self->SUPER::setup_load;
  my $dbh = $self->features_db;

  my $insert_note  = $dbh->prepare('INSERT INTO fnote (fid,fnote) VALUES (?,?)');
  my $insert_data  = $dbh->prepare(<<END);
INSERT INTO fdata (fref,fstart,fstop,fbin,ftypeid,fscore,
		   fstrand,fphase,gid,ftarget_start,ftarget_stop)
       VALUES(?,?,?,?,?,?,?,?,?,?,?)
END
;

  $self->{load_stuff}{insert_data}  = $insert_data;
  $self->{load_stuff}{insert_note}  = $insert_note;
}

sub load_gff_line {
  my $self      = shift;
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

  my $group_id;   # undef to start with

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

  unless ($group_id = $s->{groups}{$key}) {
    warn "No groupid for $group_class:$group_name; ",$dbh->errstr," Record skipped.\n";
    next;
  }

  my $bin =  bin($start,$stop,$self->{minbin});

  my $result = $s->{insert_data}->execute($ref,$start,$stop,$bin,$typeid,
					  $score,$strand,$phase,$group_id,
					  $target_start,$target_stop);
  unless ($result) {
    warn $dbh->errstr,"\n";
    next;
  }

  # add notes
  my $fid = $dbh->{mysql_insertid};
  foreach (@$notes) {
    $s->{insert_note}->execute($fid,$_);
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
  $self->{load_stuff}{insert_note}->finish;
  $self->SUPER::finish_load;
}

sub tables {
  qw(fdata fgroup ftype fdna fnote)
}

sub schema {
  return <<END;
create table fdata (
    fid	                int not null auto_increment,
    fref                varchar(20)    not null,
    fstart              int unsigned   not null,
    fstop               int unsigned   not null,
    fbin                double(20,6)  not null,
    ftypeid             int not null,
    fscore              float,
    fstrand             enum('+','-'),
    fphase              enum('0','1','2'),
    gid                 int not null,
    ftarget_start       int unsigned,
    ftarget_stop        int unsigned,
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

create table fnote (
    fid      int not null,
    fnote    text,
    index(fid)
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

