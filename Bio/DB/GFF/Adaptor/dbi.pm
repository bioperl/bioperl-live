package Bio::DB::GFF::Adaptor::dbi;

# base class for dbi-based implementations
use strict;

use DBI;
use Bio::DB::GFF;
use Bio::DB::GFF::Util::Rearrange; # for rearrange()
use vars qw($VERSION @ISA);

@ISA =  qw(Bio::DB::GFF);
$VERSION = '0.20';

use constant MAX_SEGMENT => 100_000_000;  # the largest a segment can get

# Create a new Bio::DB::GFF::Adaptor::dbi object
sub new {
  my $class = shift;
  my ($features_db,$username,$auth) = rearrange([
						 [qw(FEATUREDB DB DSN)],
						 [qw(USERNAME USER)],
						 [qw(PASSWORD PASS)]
						],@_);

  $features_db  || $class->throw("new(): Provide a data source or DBI database");

  if (!ref($features_db)) {
    my $dsn = $features_db;
    my @args;
    push @args,$username if defined $username;
    push @args,$auth     if defined $username && defined $auth;
    $features_db = DBI->connect($dsn,@args)
      || $class->throw("new(): Failed to connect to $dsn: ".DBI->errstr);
  } else {
    $features_db->isa('DBI::db') 
      || $class->throw("new(): $features_db is not a DBI handle");
  }

  # fill in object
  return bless {
		features_db => $features_db
	       },$class;
}

sub features_db { shift->{features_db} }

sub do_query {
  my $self = shift;
  my ($query,@args) = @_;
  warn $self->dbi_quote($query,@args),"\n" if $self->debug;
  my $sth = $self->{sth}{$query} ||= $self->features_db->prepare($query)
    || $self->throw("Couldn't prepare query $query:\n ".DBI->errstr."\n");
  $sth->execute(@args)
    || $self->throw("Couldn't execute query $query:\n ".DBI->errstr."\n");
  $sth;
}

# given sequence name, and optional (start,stop) give raw dna
sub get_dna {
  my $self = shift;
  my ($name,$class,$start,$stop) = @_;
  my ($offset,$length);
  if ($stop > $start) {
    $offset = $stop - 1;
    $length = $stop - $start + 1;
  } elsif ($start > $stop) {
    $offset = $stop - 1;
    $length = $start - $start + 1;
  } else {
    return;   # zero length == empty string
  }

  my $sth = $self->dna_query($name,$start,$stop,$class);
  my @row = $sth->fetchrow_array;
  $sth->finish;

  my $dna = $row[0];
  if ($stop < $start) {
    $dna = reverse $dna;
    $dna =~ tr/gatcGATC/ctagCTAG/;
  }
  $dna;
}

sub dna_query {
  shift->throw("make_dna_query(): must be implemented by a subclass");
}

# given sequence name, return (reference,start,stop,strand)
sub get_abscoords {
  my $self = shift;
  my ($name,$class)  = @_;

  my ($query,@args)  = $self->make_abscoord_query($name,$class);
  my $sth            = $self->do_query($query,@args);

  my @result;
  while ( my @row = $sth->fetchrow_array) {
    push @result,\@row
  }
  $sth->finish;

  if (@result == 0) {
    $self->error("$name not found in database");
    return;
  } elsif (@result > 1) {
    $self->error("$name present more than once in database");
    return;
  } else {
    return @{$result[0]};
  }
}

# Given sequence name, range, and optional filter, retrieve list of all
# features.  Passes features through callback.
sub get_features {
  my $self = shift;
  my ($isrange,$srcseq,$class,$start,$stop,$types,$callback) = @_;
  $callback || $self->throw('must provide a callback argument');

  my $sth = $self->range_or_overlap($isrange,$srcseq,$class,$start,$stop,$types,$class) or return;

  my $count = 0;
  while (my @row = $sth->fetchrow_array) {
    $callback->(@row);
    $count++;
  }
  $sth->finish;
  return $count;
}

sub get_features_iterator {
  my $self = shift;
  my ($isrange,$srcseq,$class,$start,$stop,$types,$callback) = @_;
  $callback || $self->throw('must provide a callback argument');

  my $sth = $self->range_or_overlap($isrange,$srcseq,$class,$start,$stop,$types,$class) or return;

  return Bio::DB::GFF::Adaptor::dbi::iterator->new($sth,$callback);
}

sub get_types {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  my $straight      = $self->do_straight_join($srcseq,$start,$stop,[]) ? 'straight_join' : '';
  my ($select,@args1) = $self->make_types_select_part($srcseq,$start,$stop,$want_count);
  my ($from,@args2)   = $self->make_types_from_part($srcseq,$start,$stop,$want_count);
  my ($join,@args3)   = $self->make_types_join_part($srcseq,$start,$stop,$want_count);
  my ($where,@args4)  = $self->make_types_where_part($srcseq,$start,$stop,$want_count);
  my ($group,@args5)  = $self->make_types_group_part($srcseq,$start,$stop,$want_count);

  my $query = "SELECT $straight $select FROM $from WHERE $join AND $where";
  $query   .= " GROUP BY $group" if $group;
  my @args  = (@args1,@args2,@args3,@args4,@args5);
  my $sth = $self->do_query($query,@args) or return;

  my (%result,%obj);
  while (my ($method,$source,$count) = $sth->fetchrow_array) {
    my $type = Bio::DB::GFF::Typename->new($method,$source);
    $result{$type} = $count;
    $obj{$type} = $type;
  }
  return $want_count ? %result : values %obj;
}

# THE FOLLOWING ROUTINES ALL PERTAIN TO RANGE AND TYPE QUERIES
# this is what will need to change if the structure of the GFF table is altered
sub range_or_overlap {
  my $self = shift;
  my($isrange,$srcseq,$class,$start,$stop,$types) = @_;

  my $dbh = $self->features_db;

  # NOTE: straight_join is necessary in some database to force the right index to be used.
  my $straight      = $self->do_straight_join($srcseq,$start,$stop,$types) ? 'straight_join' : '';
  my $select        = $self->make_features_select_part;
  my $from          = $self->make_features_from_part;
  my $join          = $self->make_features_join_part;
  my ($where,@args) = $self->make_features_where_part($isrange,$srcseq,$class,
						      $start,$stop,$types,$class);
  my $query         = "SELECT $straight $select FROM $from WHERE $join AND $where";

  my $sth = $self->do_query($query,@args);
  $sth;
}

sub do_straight_join { 0 }  # false by default

sub make_features_select_part {
  shift->throw("make_features_select_part(): must be implemented by subclass");
}

sub make_features_from_part {
  shift->throw("make_features_from_part(): must be implemented by subclass");
}

sub make_features_join_part {
  shift->throw("make_features_join_part(): must be implemented by subclass");
}

sub make_features_where_part {
  my $self = shift;
  my($isrange,$srcseq,$class,$start,$stop,$types) = @_;
  my @query;
  my @args;


  if ($srcseq) {
    my ($q,@a) = $self->srcseq_query($srcseq,$class);
    push @query,$q;
    push @args,@a;
  }

  if (defined $start or defined $stop) {
    $start = 0               unless defined($start);
    $stop  = MAX_SEGMENT     unless defined($stop);

    my ($range_query,@range_args) = $isrange ? $self->range_query($start,$stop) 
                                             : $self->overlap_query($start,$stop);
    push @query,$range_query;
    push @args,@range_args;
  }

  if (defined $types && @$types) {
    my ($type_query,@type_args) = $self->types_query($types);
    push @query,$type_query;
    push @args,@type_args;
  }

  my $query = join "\n\tAND ",@query;
  return wantarray ? ($query,@args) : $self->dbi_quote($query,@args);
}

# generate the fragment of SQL responsible for returning the
# reference sequence, start, stop and strand given a sequence class
# and name.
sub make_abscoord_query {
  my $self = shift;
  my ($seq_name,$seq_class) = @_;
  $self->throw("make_abscoord_query(): must be implemented by subclass");
  # in scalar context, return a query string.
  # in array context, return a query string and bind arguments
}


sub srcseq_query {
  my $self = shift;
  my ($srcseq,$refclass) = @_;
  $self->throw("srcseq_query(): must be implemented by subclass");
  # in scalar context, return a query string.
  # in array context, return a query string and bind arguments
}

# generate the fragment of SQL responsible for searching for
# features that overlap a given range
sub overlap_query {
  my $self = shift;
  my ($start,$stop) = @_;
  $self->throw("overlap_query(): must be implemented by subclass");
  # in scalar context, return a query string.
  # in array context, return a query string and bind arguments
}

# generate the fragment of SQL responsible for searching for
# features that are completely contained within a range
sub range_query {
  my $self = shift;
  my ($start,$stop) = @_;
  $self->throw("range_query(): must be implemented by subclass");
  # in scalar context, return a query string.
  # in array context, return a query string and bind arguments
}

# generate the fragment of SQL responsible for searching for
# features with particular types and methods
sub types_query {
  my $self  = shift;
  my $types = shift;  # array ref
  $self->throw("types_query(): must be implemented by subclass");
  # in scalar context, return a query string.
  # in array context, return a query string and bind arguments
}

sub make_types_select_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  $self->throw("make_types_select_part(): must be implemented by subclass");
}

sub make_types_from_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  $self->throw("make_types_from_part(): must be implemented by subclass");
}

sub make_types_join_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  $self->throw("make_types_join_part(): must be implemented by subclass");
}

sub make_types_where_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  $self->throw("make_types_where_part(): must be implemented by subclass");
}

sub make_types_group_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  $self->throw("make_types_group_part(): must be implemented by subclass");
}

sub string_match {
  my $self           = shift;
  my ($field,$value) = @_;
  return qq($field = ?) if $value =~ /^[!@%&a-zA-Z0-9_\'\" ~-]+$/;
  return qq($field REGEXP ?);
}

sub dbi_quote {
  my $self = shift;
  my ($query,@args) = @_;
  my $dbi = $self->features_db;
  $query =~ s/\?/$dbi->quote(shift @args)/eg;
  $query;
}

########################## loading and initialization  #####################

# Create the schema from scratch.
# You will need create privileges for this.
sub do_initialize {
  my $self = shift;
  my $drop_all = shift;
  $self->drop_all if $drop_all;

  my $dbh = $self->features_db;
  my @statements = split "\n\n",$self->schema;
  foreach (@statements) {
    s/;.*\Z//s;
    return unless $dbh->do($_);
  }
  1;
}

# Drop all the GFF tables -- dangerous!
sub drop_all {
  my $self = shift;
  my $dbh = $self->features_db;
  local $dbh->{PrintError} = 0;
  foreach ($self->tables) {
    $dbh->do("drop table $_");
  }
}

# return list of tables that "belong" to us.
sub tables {
  shift->throw("tables(): must be implemented by subclass");
}

sub DESTROY {
  my $self = shift;
  $self->features_db->disconnect if defined $self->features_db;
}

package Bio::DB::GFF::Adaptor::dbi::iterator;

sub new {
  my $class = shift;
  my ($sth,$callback) = @_;
  return bless [$sth,$callback],$class;
}

sub next_feature {
  my $self = shift;
  return unless $self->[0];
  if (my @row = $self->[0]->fetchrow_array) {
    return $self->[1]->(@row);
  } else {
    $self->[0]->finish;
    undef $self->[0];
  }
}

1;
__END__
# Below is stub documentation for your module. You better edit it!

=head1 NAME

Ace::Sequence::Mysql - Perl extension for blah blah blah

=head1 SYNOPSIS

  use Ace::Sequence::Mysql;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for Ace::Sequence::Mysql, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.


=head1 AUTHOR

A. U. Thor, a.u.thor@a.galaxy.far.far.away

=head1 SEE ALSO

perl(1).

=cut
