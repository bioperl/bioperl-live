=head1 NAME

Bio::DB::GFF::Adaptor::dbi -- Database adaptor for DBI (SQL) databases

=head1 SYNOPSIS

See L<Bio::DB::GFF>

=head1 DESCRIPTION

This is the base class for DBI-based adaptors.  It does everything
except generating the text of the queries to be used.  See the section
QUERIES TO IMPLEMENT for the list of methods that must be implemented.

=cut

package Bio::DB::GFF::Adaptor::dbi;

# base class for dbi-based implementations
use strict;

use DBI;
use Bio::DB::GFF;
use Bio::DB::GFF::Util::Rearrange; # for rearrange()
use Bio::DB::GFF::Adaptor::dbi::iterator;
use Bio::DB::GFF::Adaptor::dbi::caching_handle;
use vars qw($VERSION @ISA);

@ISA =  qw(Bio::DB::GFF);
$VERSION = '0.30';

use constant MAX_SEGMENT => 100_000_000;  # the largest a segment can get

=head2 new

 Title   : new
 Usage   : $db = Bio::DB::GFF->new(@args)
 Function: create a new adaptor
 Returns : a Bio::DB::GFF object
 Args    : see below
 Status  : Public

This is the constructor for the adaptor.  It is called automatically
by Bio::DB::GFF-E<gt>new.  In addition to arguments that are common among
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
    push @args,$auth     if defined $auth;
    $features_db = Bio::DB::GFF::Adaptor::dbi::caching_handle->new($dsn,@args)
      || $class->throw("new(): Failed to connect to $dsn: "
		       . Bio::DB::GFF::Adaptor::dbi::caching_handle->errstr);
  } else {
    $features_db->isa('DBI::db') 
      || $class->throw("new(): $features_db is not a DBI handle");
  }

  # fill in object
  return bless {
		features_db => $features_db
	       },$class;
}

sub debug {
  my $self = shift;
  $self->features_db->debug(@_);
  $self->SUPER::debug(@_);
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
sub dbh         { shift->{features_db} }

=head2 get_dna

 Title   : get_dna
 Usage   : $string = $db->get_dna($name,$start,$stop,$class)
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
  my ($name,$start,$stop,$class) = @_;
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
 Returns : an array ref -- see below
 Args    : name and class of desired landmark
 Status  : Public

This method performs the low-level resolution of a landmark into a
reference sequence and position.

The result is an array ref, each element of which is a five-element
list containing reference sequence name, class, start, stop and strand.

=cut

# given sequence name, return (reference,start,stop,strand)
sub get_abscoords {
  my $self = shift;
  my ($name,$class,$refseq)  = @_;

  my $sth = $self->make_abscoord_query($name,$class,$refseq);

  my @result;
  while ( my @row = $sth->fetchrow_array) {
    push @result,\@row
  }
  $sth->finish;

  if (@result == 0) {
    $self->error("$name not found in database");
    return;
  } else {
    return \@result;
  }
}

=head2 get_features

 Title   : get_features
 Usage   : $db->get_features($search,$options,$callback)
 Function: retrieve features from the database
 Returns : number of features retrieved
 Args    : see below
 Status  : Public

This is the low-level method that is called to retrieve GFF lines from
the database.  It is responsible for retrieving features that satisfy
range and feature type criteria, and passing the GFF fields to a
callback subroutine.

See the manual page for Bio::DB::GFF for the interpretation of the
arguments and how the information retrieved by get_features is passed
to the callback for processing.

Internally, get_features() is a front end for range_query().  The
latter method constructs the query and executes it.  get_features()
calls fetchrow_array() to recover the fields and passes them to the
callback.

=cut

# Given sequence name, range, and optional filter, retrieve list of
# all features.  Passes features through callback.
sub get_features {
  my $self = shift;
  my ($search,$options,$callback) = @_;
  $callback || $self->throw('must provide a callback argument');

  my $sth = $self->range_query(@{$search}{qw(rangetype
					     refseq
					     refclass
					     start
					     stop
					     types) },
			       @{$options}{qw(
					      sparse
					      sort_by_group
					      ATTRIBUTES)}) or return;

  my $count = 0;
  while (my @row = $sth->fetchrow_array) {
    $callback->(@row);
    $count++;
  }
  $sth->finish;
  return $count;
}

=head2 classes

 Title   : classes
 Usage   : $db->classes
 Function: return list of landmark classes in database
 Returns : a list of classes
 Args    : none
 Status  : public

This routine returns the list of reference classes known to the
database, or empty if classes are not used by the database.  Classes
are distinct from types, being essentially qualifiers on the reference
namespaces.

NOTE: In the current mysql-based schema, this query takes a while to
run due to the classes not being normalized.

=cut

sub classes {
  my $self = shift;
  my ($query,@args) = $self->make_classes_query or return;
  my $sth           = $self->dbh->do_query($query,@args);
  my @classes;
  while (my ($c) = $sth->fetchrow_array) {
     push @classes,$c;
  }
  @classes;
}

=head2 make_classes_query

 Title   : make_classes_query
 Usage   : ($query,@args) = $db->make_classes_query
 Function: return query fragment for generating list of reference classes
 Returns : a query and args
 Args    : none
 Status  : public

=cut

sub make_classes_query {
  my $self = shift;
  return;
}

=head2 _feature_by_name

 Title   : _feature_by_name
 Usage   : $db->get_features_by_name($name,$class,$callback)
 Function: get a list of features by name and class
 Returns : count of number of features retrieved
 Args    : name of feature, class of feature, and a callback
 Status  : protected

This method is used internally.  The callback arguments are those used
by make_feature().  Internally, it invokes the following abstract procedures:

 make_features_select_part
 make_features_from_part
 make_features_by_name_where_part
 make_features_join_part

=cut

sub _feature_by_name {
  my $self = shift;
  my ($class,$name,$callback) = @_;
  $callback || $self->throw('must provide a callback argument');

  my $select         = $self->make_features_select_part;
  my $from           = $self->make_features_from_part;
  my ($where,@args)  = $self->make_features_by_name_where_part($class,$name);
  my $join           = $self->make_features_join_part;
  my $query  = "SELECT $select FROM $from WHERE $where AND $join";
  my $sth    = $self->dbh->do_query($query,@args);

  my $count = 0;
  while (my @row = $sth->fetchrow_array) {
    $callback->(@row);
    $count++;
  }
  $sth->finish;
  return $count;
}

=head2 _feature_by_id

 Title   : _feature_by_id
 Usage   : $db->_feature_by_id($ids,$type,$callback)
 Function: get a list of features by ID
 Returns : count of number of features retrieved
 Args    : arrayref containing list of IDs to fetch and a callback
 Status  : protected

This method is used internally.  The $type selector is one of
"feature" or "group".  The callback arguments are those used by
make_feature().  Internally, it invokes the following abstract
procedures:

 make_features_select_part
 make_features_from_part
 make_features_by_id_where_part
 make_features_join_part

=cut

sub _feature_by_id {
  my $self = shift;
  my ($ids,$type,$callback) = @_;
  $callback || $self->throw('must provide a callback argument');

  my $select         = $self->make_features_select_part;
  my $from           = $self->make_features_from_part;
  my ($where,@args)  = $type eq 'feature' ? $self->make_features_by_id_where_part($ids)
                                          : $self->make_features_by_gid_where_part($ids);
  my $join           = $self->make_features_join_part;
  my $query          = "SELECT $select FROM $from WHERE $where AND $join";
  my $sth            = $self->dbh->do_query($query,@args);

  my $count = 0;
  while (my @row = $sth->fetchrow_array) {
    $callback->(@row);
    $count++;
  }
  $sth->finish;
  return $count;
}

sub _feature_by_attribute {
  my $self = shift;
  my ($attributes,$callback) = @_;
  $callback || $self->throw('must provide a callback argument');

  my $select         = $self->make_features_select_part;
  my $from           = $self->make_features_from_part(undef,{attributes=>$attributes});
  my ($where,@args)  = $self->make_features_by_range_where_part('',{attributes=>$attributes});
  my $join           = $self->make_features_join_part({attributes=>$attributes});
  my $query          = "SELECT $select FROM $from WHERE $where AND $join";
  my $sth            = $self->dbh->do_query($query,@args);

  my $count = 0;
  while (my @row = $sth->fetchrow_array) {
    $callback->(@row);
    $count++;
  }
  $sth->finish;
  return $count;
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
a list of $name=E<gt>$count pairs, where $count indicates the number of
times this feature occurs in the range.

Internally, this method calls upon the following functions to generate
the SQL and its bind variables:

  ($q1,@args) = make_types_select_part(@args);
  ($q2,@args) = make_types_from_part(@args);
  ($q3,@args) = make_types_where_part(@args);
  ($q4,@args) = make_types_join_part(@args);
  ($q5,@args) = make_types_group_part(@args);

The components are then combined as follows:

  $query = "SELECT $q1 FROM $q2 WHERE $q3 AND $q4 GROUP BY $q5";

If any of the query fragments contain the ? bind variable, then the
same number of bind arguments must be provided in @args.  The
fragment-generating functions are described below.

=cut

sub get_types {
  my $self = shift;
  my ($srcseq,$class,$start,$stop,$want_count,$typelist) = @_;
  my $straight      = $self->do_straight_join($srcseq,$start,$stop,[]) ? 'straight_join' : '';
  my ($select,@args1) = $self->make_types_select_part($srcseq,$start,$stop,$want_count,$typelist);
  my ($from,@args2)   = $self->make_types_from_part($srcseq,$start,$stop,$want_count,$typelist);
  my ($join,@args3)   = $self->make_types_join_part($srcseq,$start,$stop,$want_count,$typelist);
  my ($where,@args4)  = $self->make_types_where_part($srcseq,$start,$stop,$want_count,$typelist);
  my ($group,@args5)  = $self->make_types_group_part($srcseq,$start,$stop,$want_count,$typelist);

  my $query = "SELECT $straight $select FROM $from WHERE $join AND $where";
  $query   .= " GROUP BY $group" if $group;
  my @args  = (@args1,@args2,@args3,@args4,@args5);
  my $sth = $self->dbh->do_query($query,@args) or return;

  my (%result,%obj);
  while (my ($method,$source,$count) = $sth->fetchrow_array) {
    my $type = Bio::DB::GFF::Typename->new($method,$source);
    $result{$type} = $count;
    $obj{$type} = $type;
  }
  return $want_count ? %result : values %obj;
}

=head2 range_query

 Title   : range_query
 Usage   : $db->range_query($range_type,$refseq,$refclass,$start,$stop,$types,$order_by_group)
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

  $order_by_group        A flag indicating that statement handler should group
                         the features by group id (handy for iterative fetches)

If successful, this method returns a statement handle.  The handle is
expected to return the fields described for get_features().

Internally, range_query() makes calls to the following methods,
each of which is expected to be overridden in subclasses:

  $select        = $self->make_features_select_part;
  $from          = $self->make_features_from_part;
  $join          = $self->make_features_join_part;
  ($where,@args) = $self->make_features_by_range_where_part($isrange,$srcseq,$class,
						           $start,$stop,$types,$class);

The query that is constructed looks like this:

  SELECT $select FROM $from WHERE $join AND $where

The arguments that are returned from make_features_by_range_where_part() are
passed to the statement handler's execute() method.

range_query() also calls a do_straight_join() method, described
below.  If this method returns true, then the keyword "straight_join"
is inserted right after SELECT.

=cut

sub range_query {
  my $self = shift;
  my($rangetype,$refseq,$class,$start,$stop,$types,$sparse,$order_by_group,$attributes) = @_;

  my $dbh = $self->features_db;

  # NOTE: straight_join is necessary in some database to force the right index to be used.
  my %a             = (refseq=>$refseq,class=>$class,start=>$start,stop=>$stop,types=>$types,attributes=>$attributes);
  my $straight      = $self->do_straight_join(\%a) ? 'straight_join' : '';
  my $select        = $self->make_features_select_part(\%a);
  my $from          = $self->make_features_from_part($sparse,\%a);
  my $join          = $self->make_features_join_part(\%a);
  my ($where,@args) = $self->make_features_by_range_where_part($rangetype,\%a);
  my ($group_by,@more_args) = $self->make_features_group_by_part(\%a);
  my $order_by      = $self->make_features_order_by_part(\%a) if $order_by_group;

  my $query         = "SELECT $straight $select FROM $from WHERE $join";
  $query           .= " AND $where" if $where;
  if ($group_by) {
    $query           .= " GROUP BY $group_by";
    push @args,@more_args;
  }
  $query           .= " ORDER BY $order_by" if $order_by;

  my $sth = $self->dbh->do_query($query,@args);
  $sth;
}

=head2 make_features_by_range_where_part

 Title   : make_features_by_range_where_part
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

#'

sub make_features_by_range_where_part {
  my $self = shift;
  my ($rangetype,$options) = @_;
  $options ||= {};
  my ($refseq,$class,$start,$stop,$types,$attributes) =
    @{$options}{qw(refseq class start stop types attributes)};

  my (@query,@args);

  if ($refseq) {
    my ($q,@a) = $self->refseq_query($refseq,$class);
    push @query,$q;
    push @args,@a;
  }

  if (defined $start or defined $stop) {
    $start = 0               unless defined($start);
    $stop  = MAX_SEGMENT     unless defined($stop);

    my ($range_query,@range_args) =
           $rangetype eq 'overlaps'     ? $self->overlap_query($start,$stop)
	 : $rangetype eq 'contains'     ? $self->contains_query($start,$stop)
         : $rangetype eq 'contained_in' ? $self->contained_in_query($start,$stop)
         : ();

    push @query,$range_query;
    push @args,@range_args;
  }

  if (defined $types && @$types) {
    my ($type_query,@type_args) = $self->types_query($types);
    push @query,$type_query;
    push @args,@type_args;
  }

  if ($attributes) {
    my ($attribute_query,@attribute_args) = $self->make_features_by_attribute_where_part($attributes);
    push @query,"($attribute_query)";
    push @args,@attribute_args;
  }

  my $query = join "\n\tAND ",@query;
  return wantarray ? ($query,@args) : $self->dbh->dbi_quote($query,@args);
}

=head2 do_straight_join

 Title   : do_straight_join
 Usage   : $boolean = $db->do_straight_join($refseq,$class,$start,$stop,$types)
 Function: optimization flag
 Returns : a flag
 Args    : see range_query()
 Status  : Protected

This subroutine, called by range_query() returns a boolean flag.
If true, range_query() will perform a straight join, which can be
used to optimize certain SQL queries.  The four arguments correspond
to similarly-named arguments passed to range_query().

=cut

sub do_straight_join { 0 }  # false by default

=head2 string_match

 Title   : string_match
 Usage   : $string = $db->string_match($field,$value)
 Function: create a SQL fragment for performing exact or regexp string matching
 Returns : query string
 Args    : the table field and match value
 Status  : public

This method examines the passed value for meta characters.  If so it
produces a SQL fragment that performs a regular expression match.
Otherwise, it produces a fragment that performs an exact string match.

This method is not used in the module, but is available for use by
subclasses.

=cut

sub string_match {
  my $self           = shift;
  my ($field,$value) = @_;
  return qq($field = ?) if $value =~ /^[!@%&a-zA-Z0-9_\'\" ~-]+$/;
  return qq($field REGEXP ?);
}

=head2 exact_match

 Title   : exact_match
 Usage   : $string = $db->exact_match($field,$value)
 Function: create a SQL fragment for performing exact string matching
 Returns : query string
 Args    : the table field and match value
 Status  : public

This method produces the SQL fragment for matching a field name to a
constant string value.

=cut

sub exact_match {
  my $self           = shift;
  my ($field,$value) = @_;
  return qq($field = ?);
}

=head2 meta

 Title   : meta
 Usage   : $value = $db->meta($name [,$newval])
 Function: get or set a meta variable
 Returns : a string
 Args    : meta variable name and optionally value
 Status  : public

Get or set a named metavariable for the database.  Metavariables can
be used for database-specific settings.  This method calls two
class-specific methods which must be implemented:

  make_meta_get_query()   Returns a sql fragment which given a meta
                          parameter name, returns its value.  One bind
                          variable.
  make_meta_set_query()   Returns a sql fragment which takes two bind
                          arguments, the parameter name and its value


Don't make changes unless you know what you're doing!  It will affect the
persistent database.

=cut

sub meta {
  my $self = shift;
  my $param_name = uc shift;

  # getting
  if (@_) {
    my $value = shift;
    my $sql = $self->make_meta_set_query() or return;
    my $sth = $self->dbh->prepare_delayed($sql) 
              or $self->error("Can't prepare $sql: ",$self->dbh->errstr), return;
    $sth->execute($param_name,$value)
              or $self->error("Can't execute $sql: ",$self->dbh->errstr), return;
    $sth->finish;
    return $self->{meta}{$param_name} = $value;
  }

  elsif (exists $self->{meta}{$param_name}) {
    return $self->{meta}{$param_name};
  }

  else {
    undef $self->{meta}{$param_name};  # so that we don't check again
    my $sql = $self->make_meta_get_query() or return;
    my $sth  = $self->dbh->prepare_delayed($sql)
            or $self->error("Can't prepare $sql: ",$self->dbh->errstr), return;
    $sth->execute($param_name)
            or $self->error("Can't execute $sql: ",$sth->errstr),return;
    my ($value) = $sth->fetchrow_array;
    $sth->finish;
    return $self->{meta}{$param_name} = $value;
  }

}

=head2 make_meta_get_query

 Title   : make_meta_get_query
 Usage   : $sql = $db->make_meta_get_query
 Function: return SQL fragment for getting a meta parameter
 Returns : SQL fragment
 Args    : none
 Status  : public

By default this does nothing; meta parameters are not stored or
retrieved.

=cut

sub make_meta_get_query {
  return;
}

=head2 make_meta_set_query

 Title   : make_meta_set_query
 Usage   : $sql = $db->make_meta_set_query
 Function: return SQL fragment for setting a meta parameter
 Returns : SQL fragment
 Args    : none
 Status  : public

By default this does nothing; meta parameters are not stored or
retrieved.

=cut

sub make_meta_set_query {
  return;
}

=head2 get_features_iterator

 Title   : get_features_iterator
 Usage   : $iterator = $db->get_features_iterator($search,$options,$callback)
 Function: create an iterator on a features() query
 Returns : A Bio::DB::GFF::Adaptor::dbi::iterator object
 Args    : see get_features()
 Status  : public

This method is similar to get_features(), except that it returns an
iterator across the query.  See
L<Bio::DB::GFF::Adaptor::dbi::iterator>.

=cut

sub get_features_iterator {
  my $self = shift;
  my ($search,$options,$callback) = @_;
  $callback || $self->throw('must provide a callback argument');
  my $sth = $self->range_query(@{$search}{qw(rangetype
					     refseq
					     refclass
					     start
					     stop
					     types)},
			       @{$options}{qw(
					      sparse
					      sort_by_group
					      ATTRIBUTES)}) or return;
  return Bio::DB::GFF::Adaptor::dbi::iterator->new($sth,$callback);
}

########################## loading and initialization  #####################

=head2 do_initialize

 Title   : do_initialize
 Usage   : $success = $db->do_initialize($drop_all)
 Function: initialize the database
 Returns : a boolean indicating the success of the operation
 Args    : a boolean indicating whether to delete existing data
 Status  : protected

This method will load the schema into the database.  If $drop_all is
true, then any existing data in the tables known to the schema will be
deleted.

Internally, this method calls schema() to get the schema data.

=cut

# Create the schema from scratch.
# You will need create privileges for this.
sub do_initialize {
  my $self = shift;
  my $erase = shift;
  $self->drop_all if $erase;

  my $dbh = $self->features_db;
  my $schema = $self->schema;
  foreach (values %$schema) {
    $dbh->do($_) || warn $dbh->errstr;
  }
  1;
}

=head2 drop_all

 Title   : drop_all
 Usage   : $db->drop_all
 Function: empty the database
 Returns : void
 Args    : none
 Status  : protected

This method drops the tables known to this module.  Internally it
calls the abstract tables() method.

=cut

# Drop all the GFF tables -- dangerous!
sub drop_all {
  my $self = shift;
  my $dbh = $self->features_db;
  local $dbh->{PrintError} = 0;
  foreach ($self->tables) {
    $dbh->do("drop table $_");
  }
}

=head1 QUERIES TO IMPLEMENT

The following astract methods either return DBI statement handles or
fragments of SQL.  They must be implemented by subclasses of this
module.  See Bio::DB::GFF::Adaptor::dbi::mysql for examples.

=head2 make_features_select_part

 Title   : make_features_select_part
 Usage   : $string = $db->make_features_select_part()
 Function: make select part of the features query
 Returns : a string
 Args    : none
 Status  : Abstract

This abstract method creates the part of the features query that
immediately follows the SELECT keyword.

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
immediately follows the FROM keyword.

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
do not need to join, return "1".

=cut

sub make_features_join_part {
  shift->throw("make_features_join_part(): must be implemented by subclass");
}

=head2 make_features_order_by_part

 Title   : make_features_order_by_part
 Usage   : ($query,@args) = $db->make_features_order_by_part()
 Function: make the ORDER BY part of the features() query
 Returns : a SQL fragment and bind arguments, if any
 Args    : none
 Status  : Abstract

This abstract method creates the part of the features query that
immediately follows the ORDER BY part of the query issued by
features() and related methods.  Return undef if it is not necessary
to order the returned features in some way.

=cut

sub make_features_order_by_part {
  my $self = shift;
  $self->throw("please implement a make_features_order_by_part() method");
}

=head2 make_features_by_name_where_part

 Title   : make_features_by_name_where_part
 Usage   : $db->make_features_by_name_where_part($class,$name)
 Function: create the SQL fragment needed to select a feature by its group name & class
 Returns : a SQL fragment and bind arguments
 Args    : class and name of the feature group
 Status  : Protected

=cut

sub make_features_by_name_where_part {
  my $self = shift;
  my ($class,$name) = @_;
  shift->throw('make_features_by_name_where_part(): must be implemented by subclass');
}

=head2 make_features_by_id_where_part

 Title   : make_features_by_id_where_part
 Usage   : $db->make_features_by_id_where_part($ids)
 Function: create the SQL fragment needed to select a set of features by their ids
 Returns : a SQL fragment and bind arguments
 Args    : arrayref of IDs
 Status  : Protected

=cut

sub make_features_by_id_where_part {
  my $self = shift;
  my $ids  = shift;
  shift->throw('make_features_by_id_where_part(): must be implemented by subclass');
}

sub make_features_by_attribute_where_part {
  my $self = shift;
  my $ids  = shift;
  shift->throw('make_features_by_attribute_where_part(): must be implemented by subclass');
}

=head2 make_features_by_gid_where_part

 Title   : make_features_by_gid_where_part
 Usage   : $db->make_features_by_gid_where_part($ids)
 Function: create the SQL fragment needed to select a set of features by their group ids
 Returns : a SQL fragment and bind arguments
 Args    : arrayref of IDs
 Status  : Protected

=cut

sub make_features_by_gid_where_part {
  my $self = shift;
  my $ids  = shift;
  shift->throw('make_features_by_gid_where_part(): must be implemented by subclass');
}

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

This method is called by make_features_by_range_where_part() to construct the
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

This method is called by make_features_by_range_where_part() to construct the
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

=head2 contains_query

 Title   : contains_query
 Usage   : ($query,@args) = $db->contains_query($start,$stop)
 Function: create SQL fragment that selects the desired features by range
 Returns : a list containing the query and bind arguments
 Args    : the start and stop of a range, inclusive
 Status  : Abstract

This method is called by make_features_by_range_where_part() to construct the
part of the select WHERE section that selects a set of features
entirely enclosed by a range. It returns a multi-element list in which
the first element is the SQL fragment and subsequent elements are bind
values.  For example:

  sub contains_query {
     my ($start,$stop) = @_;
     return ('gff.start>=? AND gff.stop<=?',
	     $start,$stop);
  }

=cut

# generate the fragment of SQL responsible for searching for
# features that are completely contained within a range
sub contains_query {
  my $self = shift;
  my ($start,$stop) = @_;
  $self->throw("contains_query(): must be implemented by subclass");
  # in scalar context, return a query string.
  # in array context, return a query string and bind arguments
}

=head2 contained_in_query

 Title   : contained_in_query
 Usage   : ($query,@args) = $db->contained_in_query($start,$stop)
 Function: create SQL fragment that selects the desired features by range
 Returns : a list containing the query and bind arguments
 Args    : the start and stop of a range, inclusive
 Status  : Abstract

This method is called by make_features_by_range_where_part() to construct the
part of the select WHERE section that selects a set of features
entirely enclosed by a range. It returns a multi-element list in which
the first element is the SQL fragment and subsequent elements are bind
values.  For example:

  sub contained_in_query {
     my ($start,$stop) = @_;
     return ('gff.start<=? AND gff.stop>=?',
	     $start,$stop);
  }

=cut

# generate the fragment of SQL responsible for searching for
# features that are completely contained within a range
sub contained_in_query {
  my $self = shift;
  my ($start,$stop) = @_;
  $self->throw("contained_in_query(): must be implemented by subclass");
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

This method is called by make_features_by_range_where_part() to construct the
part of the select WHERE section that selects a set of features based
on their type. It returns a multi-element list in which the first
element is the SQL fragment and subsequent elements are bind values.
The argument is an array reference containing zero or more
[$method,$source] pairs.

=cut

#'

# generate the fragment of SQL responsible for searching for
# features with particular types and methods
sub types_query {
  my $self  = shift;
  my $types = shift;  # array ref
  $self->throw("types_query(): must be implemented by subclass");
  # in scalar context, return a query string.
  # in array context, return a query string and bind arguments
}

=head2 make_types_select_part

 Title   : make_types_select_part
 Usage   : ($string,@args) = $db->make_types_select_part(@args)
 Function: create the select portion of the SQL for fetching features type list
 Returns : query string and bind arguments
 Args    : see below
 Status  : Abstract

This abstract method is called by get_types() to generate the query
fragment and bind arguments for the SELECT part of the query that
retrieves lists of feature types.  The four positional arguments are
as follows:

 $refseq      reference sequence name
 $start       start of region
 $stop        end of region
 $want_count  true to return the count of this feature type

If $want_count is false, the SQL fragment returned must produce a list
of feature types in the format (method, source).

If $want_count is true, the returned fragment must produce a list of
feature types in the format (method, source, count).

=cut

sub make_types_select_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  $self->throw("make_types_select_part(): must be implemented by subclass");
}

=head2 make_types_from_part

 Title   : make_types_from_part
 Usage   : ($string,@args) = $db->make_types_from_part(@args)
 Function: create the FROM portion of the SQL for fetching features type lists
 Returns : query string and bind arguments
 Args    : see below
 Status  : Abstract

This abstract method is called by get_types() to generate the query
fragment and bind arguments for the FROM part of the query that
retrieves lists of feature types.  The four positional arguments are
as follows:

 $refseq      reference sequence name
 $start       start of region
 $stop        end of region
 $want_count  true to return the count of this feature type

If $want_count is false, the SQL fragment returned must produce a list
of feature types in the format (method, source).

If $want_count is true, the returned fragment must produce a list of
feature types in the format (method, source, count).

=cut

sub make_types_from_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  $self->throw("make_types_from_part(): must be implemented by subclass");
}

=head2 make_types_join_part

 Title   : make_types_join_part
 Usage   : ($string,@args) = $db->make_types_join_part(@args)
 Function: create the JOIN portion of the SQL for fetching features type lists
 Returns : query string and bind arguments
 Args    : see below
 Status  : Abstract

This abstract method is called by get_types() to generate the query
fragment and bind arguments for the JOIN part of the query that
retrieves lists of feature types.  The four positional arguments are
as follows:

 $refseq      reference sequence name
 $start       start of region
 $stop        end of region
 $want_count  true to return the count of this feature type

=cut

sub make_types_join_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  $self->throw("make_types_join_part(): must be implemented by subclass");
}

=head2 make_types_where_part

 Title   : make_types_where_part
 Usage   : ($string,@args) = $db->make_types_where_part(@args)
 Function: create the WHERE portion of the SQL for fetching features type lists
 Returns : query string and bind arguments
 Args    : see below
 Status  : Abstract

This abstract method is called by get_types() to generate the query
fragment and bind arguments for the WHERE part of the query that
retrieves lists of feature types.  The four positional arguments are
as follows:

 $refseq      reference sequence name
 $start       start of region
 $stop        end of region
 $want_count  true to return the count of this feature type

=cut

sub make_types_where_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  $self->throw("make_types_where_part(): must be implemented by subclass");
}

=head2 make_types_group_part

 Title   : make_types_group_part
 Usage   : ($string,@args) = $db->make_types_group_part(@args)
 Function: create the GROUP BY portion of the SQL for fetching features type lists
 Returns : query string and bind arguments
 Args    : see below
 Status  : Abstract

This abstract method is called by get_types() to generate the query
fragment and bind arguments for the GROUP BY part of the query that
retrieves lists of feature types.  The four positional arguments are
as follows:

 $refseq      reference sequence name
 $start       start of region
 $stop        end of region
 $want_count  true to return the count of this feature type

This method may return an empty list if grouping is not required.
Note, however, that grouping will always be required if a type count
is requested.

=cut

sub make_types_group_part {
  my $self = shift;
  my ($srcseq,$start,$stop,$want_count) = @_;
  $self->throw("make_types_group_part(): must be implemented by subclass");
}

=head2 tables

 Title   : tables
 Usage   : @tables = $db->tables
 Function: return list of tables that belong to this module
 Returns : list of tables
 Args    : none
 Status  : protected

This method lists the tables known to the module.

=cut

# return list of tables that "belong" to us.
sub tables {
  my $schema = shift->schema;
  return keys %$schema;
}

=head2 schema

 Title   : schema
 Usage   : $schema = $db->schema
 Function: return the CREATE script for the schema
 Returns : a hashref
 Args    : none
 Status  : abstract

This method returns an array ref containing the various CREATE
statements needed to initialize the database tables.  The keys are the
table names, and the values are strings containing the appropriate
CREATE statement.

=cut

sub schema {
  shift->throw("The schema() method must be implemented by subclass");
}

=head2 DESTROY

 Title   : DESTROY
 Usage   : $db->DESTROY
 Function: disconnect database at destruct time
 Returns : void
 Args    : none
 Status  : protected

This is the destructor for the class.

=cut

sub DESTROY {
  my $self = shift;
  $self->features_db->disconnect if defined $self->features_db;
}

################## query cache ##################


1;

__END__

=head1 BUGS

Schemas need work to support multiple hierarchical groups.

=head1 SEE ALSO

L<Bio::DB::GFF>, L<bioperl>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

