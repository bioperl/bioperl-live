package Bio::DB::GFF::Adaptor::dbi::mysqlopt;

# an optimized mysql adaptor
# For speed, this module:
# 1) uses Bio::DB::Fasta to retrieve the sequence
# 2) uses a binning scheme to make short range queries very much faster
# 3) optionally returns Acedb objects, for compatibility with WormBase

require 5.6.0;
use strict;
use Carp 'croak';
use base 'Bio::DB::GFF::Adaptor::dbi::mysql';
use Bio::DB::Fasta;
use Bio::DB::GFF::Util::Binning;
use Bio::DB::GFF::Util::Rearrange;

# this is the largest that any reference sequence can be (100 megabases)
use constant MAX_BIN    => 100_000_000;

# this is the smallest bin (1 K)
use constant MIN_BIN    => 1000;

# size of range over which it is faster to force mysql to use the range for indexing
use constant STRAIGHT_JOIN_LIMIT => 200_000;

sub new {
  my $class = shift;
  my ($features_db,$dna_db,$acedb,
      $minbin,$maxbin,
      $join_limit) =  rearrange([
				 [qw(FEATUREDB DB DSN)],
				 [qw(DNADB DNA FASTA FASTA_DIR)],
				 'ACEDB',
				 'MINBIN',
				 'MAXBIN',
				 [qw(STRAIGHT_JOIN_LIMIT JOIN_LIMIT)],
				],@_);
  my $self = $class->SUPER::new($features_db);

  if ($dna_db) {
    if (!ref($dna_db)) {
      my $fasta_dir = $dna_db;
      $dna_db = Bio::DB::Fasta->new($fasta_dir);
      $dna_db or croak __PACKAGE__."->new(): Failed to create new Bio::DB::Fasta from files in $fasta_dir";
    } else {
      $dna_db->isa('Bio::DB::Fasta') or croak __PACKAGE__."->new(): $dna_db is not a Bio::DB::Fasta object";
    }
    $self->{dna_db} = $dna_db;
  }

  if ($acedb) {
    $acedb->isa('Ace') or croak "$acedb is not an acedb accessor object";
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
  my ($name,$start,$stop) = @_;
  my $dna_db = $self->dna_db or return $self->SUPER::get_dna(@_);
  $dna_db->seq($name,$start,$stop);
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
  my ($query,@args)   = $self->bin_query($start,$stop);
  my ($iquery,@iargs) = $self->SUPER::range_query($start,$stop);
  $query .= "\n\t AND $iquery";
  push @args,@iargs;
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
  my $query = join("\n\t OR ",@bins);
  return wantarray ? ($query,@args)
                   : $self->dbi_quote($query,@args);
}

1;
