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

sub refseq_query {
  my $self = shift;
  my $ref = shift;
  my $query = "fdata.fref = ?\n";
  return wantarray ? ($query,$ref) : $self->dbi_quote($query,$ref);
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

1;
