=head1 NAME

Bio::DB::GFF::Adaptor::dbi::mysqlopt -- Optimized Bio::DB::GFF adaptor for mysql

=head1 SYNOPSIS

See L<Bio::DB::GFF>

=head1 DESCRIPTION

This adaptor is similar to Bio::DB::GFF::Adaptor::mysqlopt, except
that it implements several optimizations:

=over 4

=item 1. Binning

It uses a hierarchical binning scheme to dramatically accelerate
feature queries that use positional information.

=item 2. DNA fetching

Because mysql is slow when fetching substrings out of large text
BLOBs, this adaptor uses Bio::DB::Fasta to fetch DNA segments rapidly.
out of FASTA files.

=item 3. An ACEDB interface

Features can be linked to ACEDB objects, allowing this module to be
used as a replacement for the Ace::Sequence module.

=back

The schema is identical to Bio::DB::GFF::Adaptor::dbi, except for the
fdata table:

    fid	           feature ID (integer)
    fref           reference sequence name (string)
    fstart         start position relative to reference (integer)
    fstop          stop postion relative to reference (integer)
    fbin           bin containing this feature (float)
    ftypeid        feature type ID (integer)
    fscore         feature score (float); may be null
    fstrand        strand; one of "+" or "-"; may be null
    fphase         phase; one of 0, 1 or 2; may be null
    gid            group ID (integer)
    ftarget_start  for similarity features, the target start position (integer)
    ftarget_stop   for similarity features, the target stop position (integer)

The only difference is the "fbin" field, which indicates the interval
in which the feature is contained.  This module uses a hierarchical
set of bins, the smallest of which are 1 kb, and the largest is 100
megabases.

In the call to initialize() you can set the following options:

  -minbin        minimum value to use for binning

  -maxbin        maximum value to use for binning

  -straight_join_limit
                 size of range over which it is faster to force mysql to use the range for indexing

-minbin and -maxbin indicate the minimum and maximum sizes of
features, and are important for range query optimization.  They are
set at reasonable values -- in particular, the maximum bin size is set
to 100 megabases. Do not change them unless you know what you are
doing.

=cut

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
$VERSION = 0.80;

# this is the largest that any reference sequence can be (100 megabases)
use constant MAX_BIN    => 100_000_000;

# this is the smallest bin (1 K)
use constant MIN_BIN    => 1000;

# size of range over which it is faster to force mysql to use the range for indexing
use constant STRAIGHT_JOIN_LIMIT => 200_000;

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

  -fasta         path to a directory containing FASTA files for this database
                    (e.g. "/usr/local/share/fasta")

  -acedb         an acedb URL to use when converting features into ACEDB
                    objects (e.g. sace://localhost:2005)

  -user          username for authentication

  -pass          the password for authentication

  -minbin        minimum value to use for binning

  -maxbin        maximum value to use for binning

The path indicated by -fasta must be writable by the current process.
This is needed in order to build an index of the fasta files.

-minbin and -maxbin indicate the minimum and maximum sizes of
features, and are important for range query optimization.  They are
set at reasonable values -- in particular, the maximum bin size is set
to 100 megabases. Do not change them unless you know what you are
doing.

=cut


sub new {
  my $class = shift;
  my ($dna_db,$acedb,$other) =  rearrange([
					[qw(DNADB DNA FASTA FASTA_DIR)],
					'ACEDB'
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

  $self;
}

sub dna_db      { shift->{dna_db}      }
sub acedb       { shift->{acedb}       }


=head2 freshen_ace

 Title   : freshen
 Usage   : $flag = Bio::DB::GFF->freshen_ace;
 Function: Refresh internal acedb handle
 Returns : flag if correctly freshened
 Args    : none
 Status  : Public

ACeDB has an annoying way of timing out, leaving dangling database
handles.  This method will invoke the ACeDB reopen() method, which
causes dangling handles to be refreshed.  It has no effect if you are
not using ACeDB to create ACeDB objects.

=cut

sub freshen_ace {
  my $acedb = shift->acedb or return;
  $acedb->reopen();
}

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
  my($srcseq,$class,$start,$stop,$types) = @_;

  # Might try turning on and off straight join based on the number of types
  # specified, but this turns out to be very difficult indeed!

  # if a list of types has been specified, then it is almost always faster
  # to let the query optimizer figure it out.
  # the exception is when a type of "similarity" has been specified, in which
  # case the range query is better.  (yes, definitely a hack)
  # return 0 if defined($types) and @$types > 0 
  # and !grep {$_->[0] =~ /similarity/ } @$types;

  # if no types are specified then it is faster to do a range search, up to a point.
  return $srcseq && defined($start) && defined($stop) && 
    abs($stop-$start) < $self->straight_join_limit;
}


# find features that overlap a given range
sub overlap_query {
  my $self = shift;
  my ($start,$stop) = @_;

  my ($bq,@bargs)   = $self->bin_query($start,$stop);
  my ($iq,@iargs) = $self->SUPER::overlap_query($start,$stop);
  my $query = "($bq)\n\tAND $iq";
  my @args  = (@bargs,@iargs);

  return wantarray ? ($query,@args) : $self->dbh->dbi_quote($query,@args);
}

# find features that are completely contained within a range
sub contains_query {
  my $self = shift;
  my ($start,$stop) = @_;
  my ($bq,@bargs)   = $self->bin_query($start,$stop,undef,bin($start,$stop,$self->min_bin));
  my ($iq,@iargs)   = $self->SUPER::contains_query($start,$stop);
  my $query = "($bq)\n\tAND $iq";
  my @args  = (@bargs,@iargs);
  return wantarray ? ($query,@args) : $self->dbh->dbi_quote($query,@args);
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

# find features that are completely contained within a range
sub contained_in_query {
  my $self = shift;
  my ($start,$stop) = @_;
  my ($bq,@bargs)   = $self->bin_query($start,$stop,abs($stop-$start)+1,undef);
  my ($iq,@iargs)   = $self->SUPER::contained_in_query($start,$stop);
  my $query = "($bq)\n\tAND $iq";
  my @args  = (@bargs,@iargs);
  return wantarray ? ($query,@args) : $self->dbh->dbi_quote($query,@args);
}

sub make_object {
  my $self = shift;
  my ($class,$name,$start,$stop) = @_;

  if (my $db = $self->acedb) {

    # for Notes we just return a text, no database associated
    return $class->new(Text=>$name) if $class eq 'Note';

    # for homols, we create the indicated Protein or Sequence object
    # then generate a bogus Homology object (for future compatability??)
    if ($start ne '') {
      require Ace::Sequence::Homol;
      return Ace::Sequence::Homol->new($class,$name,$db,$start,$stop);
    }

    # General case:
    my $obj = $db->class->new($class=>$name,$self->acedb);

    return $obj if defined $obj;

    # Last resort, return a Text
    return $class->new(Text=>$name);
  }

  return $self->SUPER::make_object($class,$name,$start,$stop);
}

# IMPORTANT NOTE:
# WHETHER OR NOT THIS WORKS IS CRITICALLY DEPENDENT ON THE RELATIVE MAGNITUDE OF THE
sub make_features_from_part {
  my $self = shift;
  my $sparse = shift;
  my $options = shift || {};
  my $index = $sparse ? ' USE INDEX(ftypeid)': '';
  return $options->{attributes} ? "fdata${index},ftype,fgroup,fattribute,fattribute_to_feature\n"
                                : "fdata${index},ftype,fgroup\n";
}

sub bin_query {
  my $self = shift;
  my ($start,$stop,$minbin,$maxbin) = @_;
  my ($query,@args);

  $start = 0               unless defined($start);
  $stop  = $self->meta('max_bin') unless defined($stop);

  my @bins;
  $minbin = defined $minbin ? $minbin : $self->min_bin;
  $maxbin = defined $maxbin ? $maxbin : $self->max_bin;
  my $tier = $maxbin;
  while ($tier >= $minbin) {
    my ($tier_start,$tier_stop) = (bin_bot($tier,$start),bin_top($tier,$stop));
    if ($tier_start == $tier_stop) {
      push @bins,'fbin=?';
      push @args,$tier_start;
    } else {
      push @bins,'fbin between ? and ?';
      push @args,($tier_start,$tier_stop);
    }
    $tier /= 10;
  }
  $query = join("\n\t OR ",@bins);
  return wantarray ? ($query,@args)
                   : $self->dbh->dbi_quote($query,@args);
}

########################## loading and initialization  #####################

# exactly like mysql.pm setup_load except that it inserts the bin field
sub setup_load {
  my $self      = shift;

  $self->SUPER::setup_load;
  my $dbh = $self->features_db;

  my $insert_data  = $dbh->prepare(<<END);
REPLACE INTO fdata (fref,fstart,fstop,fbin,ftypeid,fscore,
		   fstrand,fphase,gid,ftarget_start,ftarget_stop)
       VALUES(?,?,?,?,?,?,?,?,?,?,?)
END
;

  my @tables = map {"$_ WRITE"} $self->tables;
  my $lock_tables = join ', ',@tables;
#   $dbh->do("LOCK TABLES $lock_tables");

  $self->{load_stuff}{sth}{insert_data}  = $insert_data;
}

sub load_gff_line {
  my $self = shift;
  my $gff = shift;

  my $s    = $self->{load_stuff};
  my $dbh  = $self->features_db;
  local $dbh->{PrintError} = 0;

  defined(my $typeid  = $self->get_table_id('ftype', $gff->{method} => $gff->{source})) or return;
  defined(my $groupid = $self->get_table_id('fgroup',$gff->{gname}  => $gff->{gclass})) or return;

  my $bin =  bin($gff->{start},$gff->{stop},$self->min_bin);
  my $result = $s->{sth}{insert_data}->execute($gff->{ref},
					       $gff->{start},$gff->{stop},$bin,
					       $typeid,
					       $gff->{score},$gff->{strand},$gff->{phase},
					       $groupid,
					       $gff->{tstart},$gff->{tstop});

  warn $dbh->errstr,"\n" and return unless $result;

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

sub schema {
  my $self = shift;
  my $schema = $self->SUPER::schema;
  $schema->{fdata} = q{
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
    gid                 int not null,
    ftarget_start       int unsigned,
    ftarget_stop        int unsigned,
    primary key(fid),
    unique index(fref,fbin,fstart,fstop,ftypeid,gid),
    index(ftypeid),
    index(gid)
)
		      };
  $schema;
}

1;

__END__

