=head1 NAME

Bio::DB::GFF::Adaptor::dbi -- Database adaptor for DBI (SQL) databases

=head1 SYNOPSIS

See L<Bio::DB::GFF>

=head1 DESCRIPTION

This is the base class for DBI-based adaptors.  It does everything
except generating the text of the queries to be used.  See the section
QUERIES TO IMPLEMENT for the list of queries that must be implemented.

=cut

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

=head2 new

 Title   : new
 Usage   : $db = Bio::DB::GFF->new(-adaptor => 'dbi:some_subclass',@args)
 Function: create a new adaptor
 Returns : a Bio::DB::GFF object
 Args    : see below
 Status  : Public

This is the constructor for the adaptor.  It is called automatically
by Bio::DB::GFF->new.  In addition to arguments that are common among
all adaptors, the following class-specific arguments are recgonized:

  Argument       Description
  --------       -----------

  -dsn           the DBI data source, e.g. 'dbi:mysql:ens0040'

  -user          username for authentication

  -pass          the password for authentication

=cut

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

=head2 features_db

 Title   : features_db
 Usage   : $dbh = $db->features_db
 Function: get database handle
 Returns : a DBI handle
 Args    : none
 Status  : Public

=cut

sub features_db { shift->{features_db} }


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

=cut

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

=head2 get_dna

 Title   : get_dna
 Usage   : $string = $db->get_dna($name,$class,$start,$stop)
 Function: get DNA string
 Returns : a string
 Args    : name, class, start and stop of desired segment
 Status  : Public

This method performs the low-level fetch of a DNA substring given its
name, class and the desired range.  It is actually a front end to the
abstract method make_dna_query(), which it calls after some argument
consistency checking.

=cut

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

  my $sth = $self->make_dna_query($name,$class,$start,$stop);
  my @row = $sth->fetchrow_array;
  $sth->finish;

  my $dna = $row[0];
  if ($stop < $start) {
    $dna = reverse $dna;
    $dna =~ tr/gatcGATC/ctagCTAG/;
  }
  $dna;
}

=head2 get_abscoords

 Title   : get_abscoords
 Usage   : ($refseq,$refclass,$start,$stop,$strand) = $db->get_abscoords($name,$class)
 Function: get absolute coordinates for landmark
 Returns : five-element list containing reference sequence name, class, start, stop and strand
 Args    : name and class of desired landmark
 Status  : Public

This method performs the low-level resolution of a landmark into a
reference sequence and position.

=cut

# given sequence name, return (reference,start,stop,strand)
sub get_abscoords {
  my $self = shift;
  my ($name,$class)  = @_;

  my $sth = $self->make_abscoord_query($name,$class);

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

=head2 get_features

 Title   : get_features
 Usage   : $db->get_features(@args)
 Function: retrieve features from the database
 Returns : number of features retrieved
 Args    : see below
 Status  : Public

This is the low-level method that is called to retrieve GFF lines from
the database.  It is responsible for retrieving features that satisfy
range and feature type criteria, and passing the GFF fields to a
callback subroutine.

The seven arguments are as follows:

  $isrange    true if this is a range query, false if an overlap query
              In a range query all features are completely contained between
              start and stop (inclusive).  An overlap query requests any feature
              that overlaps start and stop.

  $refseq     The reference sequence for this range.  This may be undef if no
              range restriction is desired.

  $refclass   The class of the reference sequence.  This can be ignored by
              database that do not recognize different types of reference
              sequence.

  $start      Start of range.  May be undef.

  $stop       Stop of range.  May be undef.

  $types      An array reference containing a list of [method,source] pairs.  If
              passed, only features that match method and/or source are requested.
              May be empty or undef.

  $callback   A code reference.  As features are retrieved, they are
              passed to this callback.

The callback expects thirteen arguments, which can be parsed directly out of GFF files:

  $refseq       The reference sequence
  $start        Start position
  $stop         Stop position
  $source       Feature source
  $method       Feature method
  $score        Feature score
  $strand       Feature strand
  $phase        Feature phase
  $group_class  Group class (may be undef)
  $group_name   Group name  (may be undef)
  $target_start Start of homology hit (target coordinates; may be undef)
  $target_stop  Stop of homology hit (target coordinates; may be undef)
  $feature_id   Unique feature ID (may be undef)

The group class, group name, target start, and target stop are all
variants of the GFF "group" field, and are optional. The target start
and stop fields are used by similarity match lines, and indicate the
coordinates of the match in the target sequence.

The feature ID, if provided, is a unique identifier of the feature
line.  The module does not depend on this ID in any way, but it is
available via Bio::DB::GFF->id() if wanted.  In the dbi::mysql and
dbi::mysqlopt adaptor, the ID is a unique row ID.  In the acedb
adaptor it is not used.

Internally, get_features() is a front end for range_or_overlap().  The
latter method constructs the query and executes it.  get_features()
calls fetchrow_array() to recover the fields and passes them to the
callback.

=cut

# Given sequence name, range, and optional filter, retrieve list of
# all features.  Passes features through callback.
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

=head2 get_features_iterator

 Title   : get_features_iterator
 Usage   : $db->get_features_iterator(@args)
 Function: get an iterator on a features query
 Returns : a Bio::SeqIO object
 Args    : as per get_features()
 Status  : Public

This method takes the same arguments as get_features(), but returns an
iterator that can be used to fetch features sequentially, as per
Bio::SeqIO.

Internally, this method is simply a front end to range_or_overlap().
The latter method constructs and executes the query, returning a
statement handle. This routine passes the statement handle to the
constructor for the iterator, along with the callback.

=cut

sub get_features_iterator {
  my $self = shift;
  my ($isrange,$srcseq,$class,$start,$stop,$types,$callback) = @_;
  $callback || $self->throw('must provide a callback argument');

  my $sth = $self->range_or_overlap($isrange,$srcseq,$class,$start,$stop,$types,$class) or return;

  return Bio::DB::GFF::Adaptor::dbi::iterator->new($sth,$callback);
}

=head2 get_types

 Title   : get_types
 Usage   : $db->get_types($refseq,$refclass,$start,$stop,$count)
 Function: get list of types
 Returns : a list of Bio::DB::GFF::Typename objects
 Args    : see below
 Status  : Public

This method is responsible for fetching the list of feature type names
from the database.  The query may be limited to a particular range, in
which case the range is indicated by a landmark sequence name and
class and its subrange, if any.  These arguments may be undef if it is
desired to retrieve all feature types in the database (which may be a
slow operation in some implementations).

If the $count flag is false, the method returns a simple list of
vBio::DB::GFF::Typename objects.  If $count is true, the method returns
a list of $name=>$count pairs, where $count indicates the number of
times this feature occurs in the range.

=cut

sub get_types {
  my $self = shift;
  my ($srcseq,$class,$start,$stop,$want_count) = @_;
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

=head2 range_or_overlap

 Title   : range_or_overlap
 Usage   : $db->range_or_overlap($isrange,$refseq,$refclass,$start,$stop,$types)
 Function: create statement handle for range/overlap queries
 Returns : a DBI statement handle
 Args    : see below
 Status  : Protected

This method constructs the statement handle for this module's central
query: given a range and/or a list of feature types, fetch their GFF
records.

The six positional arguments are as follows:

  Argument               Description

  $isrange               A flag indicating that this is a range.
			 query.  Otherwise an overlap query is
			 assumed.

  $refseq		 The reference sequence name (undef if no range).

  $refclass		 The reference sequence class (undef if no range).

  $start		 The start of the range (undef if none).

  $stop                  The stop of the range (undef if none).

  $types                 Array ref containing zero or feature types in the
			 format [method,source].

If successful, this method returns a statement handle.  The handle is
expected to return the fields described for get_features().

Internally, range_or_overlap() makes calls to the following methods,
each of which is expected to be overridden in subclasses:

  $select        = $self->make_features_select_part;
  $from          = $self->make_features_from_part;
  $join          = $self->make_features_join_part;
  ($where,@args) = $self->make_features_where_part($isrange,$srcseq,$class,
						   $start,$stop,$types,$class);

The query that is constructed looks like this:

  SELECT $select FROM $from WHERE $join AND $where

The arguments that are returned from make_features_where_part() are
passed to the statement handler's execute() method.

range_or_overlap() also calls a do_straight_join() method, described
below.  If this method returns true, then the keyword "straight_join"
is inserted right after SELECT.

=cut

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

=head2 make_features_select_part

 Title   : make_features_select_part
 Usage   : $string = $db->make_features_select_part()
 Function: make select part of the features query
 Returns : a string
 Args    : none
 Status  : Abstract

This abstract method creates the part of the features query that
immediately follows the SELECT keyword.  See
Bio::DB::Adaptor::dbi::mysql for an example.

=cut

sub make_features_select_part {
  shift->throw("make_features_select_part(): must be implemented by subclass");
}

=head2 make_features_from_part

 Title   : make_features_from_part
 Usage   : $string = $db->make_features_from_part()
 Function: make from part of the features query
 Returns : a string
 Args    : none
 Status  : Abstract

This abstract method creates the part of the features query that
immediately follows the FROM keyword.  See
Bio::DB::Adaptor::dbi::mysql for an example.

=cut

sub make_features_from_part {
  shift->throw("make_features_from_part(): must be implemented by subclass");
}

=head2 make_features_join_part

 Title   : make_features_join_part
 Usage   : $string = $db->make_features_join_part()
 Function: make join part of the features query
 Returns : a string
 Args    : none
 Status  : Abstract

This abstract method creates the part of the features query that
immediately follows the WHERE keyword.  It is combined with the output
of make_feautres_where_part() to form the full WHERE clause.  If you
do not need to join, return "1".  See Bio::DB::Adaptor::dbi::mysql for
an example.

=cut

sub make_features_join_part {
  shift->throw("make_features_join_part(): must be implemented by subclass");
}

=head2 make_features_where_part

 Title   : make_features_where_part
 Usage   : ($string,@args) =
     $db->make_features_select_part($isrange,$refseq,$class,$start,$stop,$types)
 Function: make where part of the features query
 Returns : the list ($query,@bind_args)
 Args    : see below
 Status  : Protected

This method creates the part of the features query that immediately
follows the WHERE keyword and is ANDed with the string returned by
make_features_join_part().

The six positional arguments are a flag indicating whether to perform
a range search or an overlap search, the reference sequence, class,
start and stop, all of which define an optional range to search in,
and an array reference containing a list [$method,$souce] pairs.

The method result is a multi-element list containing the query string
and the list of runtime arguments to bind to it with the execute()
method.

This method's job is to clean up arguments and perform consistency
checking.  The real work is done by the following abstract methods:

  Method             Description

  refseq_query()     Return the query string needed to match the reference
		     sequence.

  range_query()	     Return the query string needed to find all features contained
		     within a range.

  overlap_query()    Return the query string needed to find all features that overlap
		     a range.

See Bio::DB::Adaptor::dbi::mysql for an example of how this works.

=cut

sub make_features_where_part {
  my $self = shift;
  my($isrange,$refseq,$class,$start,$stop,$types) = @_;
  my @query;
  my @args;


  if ($refseq) {
    my ($q,@a) = $self->refseq_query($refseq,$class);
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

=head2 do_straight_join

 Title   : do_straight_join
 Usage   : $boolean = $db->do_straight_join($refseq,$start,$stop,$types)
 Function: optimization flag
 Returns : a flag
 Args    : see range_or_overlap()
 Status  : Protected

This subroutine, called by range_or_overlap() returns a boolean flag.
If true, range_or_overlap() will perform a straight join, which can be
used to optimize certain SQL queries.  The four arguments correspond
to similarly-named arguments passed to range_or_overlap().

=cut

sub do_straight_join { 0 }  # false by default

=head1 QUERIES TO IMPLEMENT

The following astract methods either return DBI statement handles or
fragments of SQL.  They must be implemented by subclasses of this
module.  See Bio::DB::GFF::Adaptor::dbi::mysql for examples.

=head2 make_dna_query

 Title   : make_dna_query
 Usage   : $sth = $db->make_dna_query($name,$class,$start,$stop);
 Function: create query that returns raw DNA sequence from database
 Returns : a DBI statement handle
 Args    : reference sequence name and class, and a range
 Status  : Abstract

The statement handler should return rows containing just one field,
the extracted DNA string.

=cut

sub make_dna_query {
  shift->throw("make_dna_query(): must be implemented by a subclass");
}

=head2 make_abscoord_query

 Title   : make_abscoord_query
 Usage   : $sth = $db->make_abscoord_query($name,$class);
 Function: create query that finds the reference sequence coordinates given a landmark & classa
 Returns : a DBI statement handle
 Args    : name and class of landmark
 Status  : Abstract

The statement handler should return rows containing five fields:

  1. reference sequence name
  2. reference sequence class
  3. start position
  4. stop position
  5. strand ("+" or "-")

If the database does not recognize different classes of reference
sequence, return "Sequence" as the class.

=cut

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

=head2 refseq_query

 Title   : refseq_query
 Usage   : ($query,@args) = $db->refseq_query($name,$class)
 Function: create SQL fragment that selects the desired reference sequence
 Returns : a list containing the query and bind arguments
 Args    : reference sequence name and class
 Status  : Abstract

This method is called by make_features_where_part() to construct the
part of the select WHERE section that selects a particular reference
sequence.  It returns a mult-element list in which the first element
is the SQL fragment and subsequent elements are bind values.  For
example:

  sub refseq_query {
     my ($name,$class) = @_;
     return ('gff.refseq=? AND gff.refclass=?',
	     $name,$class);
  }

=cut

sub refseq_query {
  my $self = shift;
  my ($refseq,$refclass) = @_;
  $self->throw("refseq_query(): must be implemented by subclass");
  # in scalar context, return a query string.
  # in array context, return a query string and bind arguments
}


=head2 overlap_query

 Title   : overlap_query
 Usage   : ($query,@args) = $db->overlap_query($start,$stop)
 Function: create SQL fragment that selects the desired features by range
 Returns : a list containing the query and bind arguments
 Args    : the start and stop of a range, inclusive
 Status  : Abstract

This method is called by make_features_where_part() to construct the
part of the select WHERE section that selects a set of features that
overlap a range. It returns a multi-element list in which the first
element is the SQL fragment and subsequent elements are bind values.
For example:

  sub overlap_query {
     my ($start,$stop) = @_;
     return ('gff.stop>=? AND gff.start<=?',
	     $start,$stop);
  }

=cut

# generate the fragment of SQL responsible for searching for
# features that overlap a given range
sub overlap_query {
  my $self = shift;
  my ($start,$stop) = @_;
  $self->throw("overlap_query(): must be implemented by subclass");
  # in scalar context, return a query string.
  # in array context, return a query string and bind arguments
}

=head2 range_query

 Title   : range_query
 Usage   : ($query,@args) = $db->range_query($start,$stop)
 Function: create SQL fragment that selects the desired features by range
 Returns : a list containing the query and bind arguments
 Args    : the start and stop of a range, inclusive
 Status  : Abstract

This method is called by make_features_where_part() to construct the
part of the select WHERE section that selects a set of features
entirely enclosed by a range. It returns a multi-element list in which
the first element is the SQL fragment and subsequent elements are bind
values.  For example:

  sub range_query {
     my ($start,$stop) = @_;
     return ('gff.start>=? AND gff.stop<=?',
	     $start,$stop);
  }

=cut

# generate the fragment of SQL responsible for searching for
# features that are completely contained within a range
sub range_query {
  my $self = shift;
  my ($start,$stop) = @_;
  $self->throw("range_query(): must be implemented by subclass");
  # in scalar context, return a query string.
  # in array context, return a query string and bind arguments
}

=head2 types_query

 Title   : types_query
 Usage   : ($query,@args) = $db->types_query($types)
 Function: create SQL fragment that selects the desired features by type
 Returns : a list containing the query and bind arguments
 Args    : an array reference containing the types
 Status  : Abstract

This method is called by make_features_where_part() to construct the
part of the select WHERE section that selects a set of features based
on their type. It returns a multi-element list in which the first
element is the SQL fragment and subsequent elements are bind values.
The argument is an array reference containing zero or more
[$method,$source] pairs.  See Bio::DB::GFF::Adaptor::dbi::mysql for
examples.

=cut

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
use Bio::SeqIO;
use vars '@ISA';
@ISA = 'Bio::SeqIO';

*next_seq = \&next_feature;

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
