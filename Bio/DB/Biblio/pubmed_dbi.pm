package Bio::DB::Biblio::pubmed_dbi;
use DBI;
use Storable qw(freeze thaw);
use Compress::LZF qw(:compress);
use Convert::UU qw(uudecode uuencode);
use Text::Abbrev;
use Storable qw(freeze thaw);

use vars qw(@ISA @AUTO_ATTRIBUTES @OTHER_ATTRIBUTES %SYNONYMS);
use vars qw($WHAT %WHAT $CREATE_MODE);
use strict;
use Bio::Root::AutoClass;
use Bio::Biblio;
@ISA = qw(Bio::Root::AutoClass Bio::Biblio);
$WHAT='refs';
%WHAT=abbrev qw(pubmed_ids pmids ids references refs ref);
$CREATE_MODE=0666;

BEGIN {
  @AUTO_ATTRIBUTES=qw(access dsn dbh user password collection read_only create
		     _needs_disconnect _collection_id _iteration_sth _iteration_end);
  @OTHER_ATTRIBUTES=qw(what);
  %SYNONYMS=(get_collection_id=>'collection_id');
  Bio::Root::AutoClass::declare(__PACKAGE__,\@AUTO_ATTRIBUTES,\%SYNONYMS);
}

sub _init_self {
  my($self,$class,$args)=@_;
  return unless $class eq __PACKAGE__; # to prevent subclasses from re-running this
  # This code adapted from Tie::DBI by Lincoln Stein
  my $dsn=$self->dsn;
  my $dbh=$self->dbh;
  $self->throw("Exactly one of -dsn or -dbh should be specified") unless $dsn || $dbh;
  $self->throw("Exactly one of -dsn or -dbh should be specified") if $dsn && $dbh;
  if ($dsn) {
    $dsn = "DBI:$dsn" unless $dsn=~ /^dbi/i;
    my($driver) = $dsn =~ /\w+:(\w+)/;
    # Try to establish connection with data source.
    $dbh = DBI->connect($dsn,$self->user,$self->password,
			{ AutoCommit=>1,
			  ChopBlanks=>1,
			  PrintError=>0,
			  Warn=>0,
			});
    $self->dsn($dsn);
    $self->dbh($dbh);
    $self->_needs_disconnect(1);
    $self->throw("Can't open $dsn: ".DBI->errstr) unless $dbh;
  } elsif ($dbh) {
    $dbh->{Warn} = $self->{WARN};
  }
  my $collection=$dbh->quote($self->collection);
  my $ret;
  my $create=$self->create;
  my @row=$dbh->selectrow_array
    (qq(select count(*) from collection where 
	collection_name = $collection));
  my $exists=$row[0];
  $self->throw("Collection ".$self->collection." does not exist") if $create eq 0 && $exists;
  $create=1 if !defined $create && !$exists; # create if doesn't exist
  
  if ($create) {
    if ($exists) {		# delete objects from existing collection
      $ret=$dbh->do
	(qq(delete object from collection, object where
	    collection_name = $collection and 
	    collection.collection_id = object.collection_id));
    } else {			# creat new collection
      $ret=$dbh->do
	(qq(insert into collection (collection_name)
	    values ($collection)));
    }
  }
  my @row=$dbh->selectrow_array
    (qq(select collection_id from collection where collection_name = $collection));
  my $collection_id=$row[0];
  $self->_collection_id($collection_id);
}
sub what {
  my $self=shift @_;
  my $what;
  if (@_) {
    $what=_what($_[0]) ||  $self->throw("Invalid what $_[0]");
    $self->{'WHAT'}=$what;
  } else {
    $what=$self->{'WHAT'} || $WHAT;
  }
  $what;
}
sub _what {
  my($what_param)=@_;
  my $what=$WHAT{lc($what_param)};
  $what='ids' if $what=~/pubmed_ids|pmids|ids/;
  $what='refs' if $what=~/ref/;
  $what;
}
sub collection_id {
  my $self=shift @_;
  $self->collection(@_);
}
sub count {
  my $self=shift @_;
  return $self->{'COUNT'}=$_[0] if @_;
  return $self->{'COUNT'} if defined $self->{'COUNT'};
  # else have to calculate it...
  my $dbh=$self->dbh;
  my $collection=$dbh->quote($self->collection);
  my @row=$dbh->selectrow_array
    (qq(select count(*) from object, collection where 
	collection_name = $collection and
	collection.collection_id = object.collection_id));
  my $count=$row[0];
  $self->{'COUNT'}=$count;
}
sub put {
  my($self,@objects)=@_;
  $self->throw("put called on collection opened for read-only") if $self->read_only;
  my $dbh=$self->dbh;
  my $collection_id=$self->_collection_id;
  for my $object (@objects) {
    my $id=$object->pmid;
    next unless defined $id;
    $id=$dbh->quote($id);
    my $value=$dbh->quote(uuencode(freeze($object)));
    my $ret=$dbh->do
      (qq(insert into object (collection_id, object_external_identifier, object)
	  values ($collection_id, $id, $value)));
  }
  $self->count(undef)		# cached version of count no longer valid
}
sub del {
  my($self,@args)=@_;
  $self->throw("del called on collection opened for read-only") if $self->read_only;
  my $dbh=$self->dbh;
  my $collection_id=$self->_collection_id;
  for my $arg (@args) {
    my $id=!ref($arg)? $arg: $arg->can('pmid')? $arg->pmid: $arg->can('id')? $arg->id: undef;
    next unless defined $id;
    $id=$dbh->quote($id);
    my $ret=$dbh->do
      (qq(delete from object where
	  object_external_identifier = $id and
	  collection_id = $collection_id));
  }
  $self->count(undef);		# cached version of count no longer valid
}
sub del_all {
  my($self,@args)=@_;
  $self->throw("del_all called on collection opened for read-only") if $self->read_only;
  my $dbh=$self->dbh;
  my $collection_id=$self->_collection_id;
  my $ret=$dbh->do (qq(delete from object where collection_id = $collection_id));
  $self->count(undef);		# cached version of count no longer valid
}
sub get_all {
  my $self=shift @_;
  my $what=@_? _what($_[0]) || $self->throw("Invalid what $_[0]"): $self->what;
  my $dbh=$self->dbh;
  my $collection_id=$self->_collection_id;
  my $column=$what eq 'refs'? 'object': 'object_external_identifier';
  my $rows=$dbh->selectall_arrayref 
    (qq(select $column from object where collection_id = $collection_id));
  my @results=map {$_->[0]} @$rows;
  @results=map {my $s=uudecode($_); thaw($s)} @results if $what eq 'refs';
  @results;
}
sub get_by_id {
  my($self,@ids)=@_;
  my $what=$self->what;
  my $dbh=$self->dbh;
  my $collection_id=$self->_collection_id;
  my $column=$what eq 'refs'? 'object': 'object_external_identifier';
  my @results;
  for my $id (@ids) {
    $id=$dbh->quote($id);
    my @row=$dbh->selectrow_array
      (qq(select $column from object where 
	  object_external_identifier = $id and
	  collection_id = $collection_id));
    next unless @row;
    if ($what eq 'refs') {
      my $freeze=uudecode($row[0]); # note uudecode returns status info, too, in list context
      push(@results,thaw($freeze));
    } else  {
      push(@results,$row[0]);
    }
  }
  wantarray? @results: $results[0];
}
sub get_by_ids {my $self=shift @_; $self->get_by_id(@_);}

sub get_next {
  my $self=shift @_;
  my $n=@_? $_[0]: 1;
  my $what=$self->what;
  my @results;
  unless ($self->_iteration_end) {
    my $sth=$self->_iteration_sth;
    unless ($sth) {
      my $dbh=$self->dbh;
      my $collection_id=$self->_collection_id;
      my $column=$what eq 'refs'? 'object': 'object_external_identifier';
      my $statement=(qq(select $column from object where collection_id = $collection_id));
      $sth=$dbh->prepare($statement);
      $sth->execute;
      $self->_iteration_sth($sth);
    }
    while ($n-->0) {
      my @row=$sth->fetchrow_array;
      $self->_iteration_end(1), last unless @row;
      if ($what eq 'refs') {
	my $freeze=uudecode($row[0]); # note uudecode returns status info, too, in list context
	push(@results,thaw($freeze));
      } else  {
	push(@results,$row[0]);
      }
    }
  }
  wantarray? @results: $results[0];
}
sub reset {
  my $self=shift @_;
  @_ and $self->what($_[0]);
  my $dbh=$self->dbh;
  $self->_iteration_end(0);
  $self->_iteration_sth(undef);
  $self;
}
sub get_more {my $self=shift @_; $self->get_next(@_);}
sub reset_retrieval {my $self=shift @_; $self->reset(@_);}


sub DESTROY {
  $_[0]->dbh->disconnect if defined $_[0]->dbh && $_[0]->_needs_disconnect;
}

# Code below is adapted from Bio::Biblio::IO::medlinexml by Martin Senger
sub fh {
  my $self = shift;
  my $class = ref($self) || $self;
  my $s = Symbol::gensym;
  tie $$s,$class,$self;
  return $s;
}
sub TIEHANDLE {
  my ($class,$val) = @_;
  return bless {'biblio' => $val}, $class;
}
sub READLINE {
  my $self = shift;
  return $self->{'biblio'}->get_next() unless wantarray;
  my (@list, $obj);
  push @list, $obj while $obj = $self->{'biblio'}->get_next();
  return @list;
}
sub PRINT {
  my $self = shift;
  $self->{'biblio'}->put(@_);
}
sub CLOSE {
  my $self=shift @_;
  my $dbh=$self->{'biblio'}->dbh;
  $self->{'biblio'}->dbh(undef);
}
sub UNTIE {}

# TIEHASH interface
# tie %h,Bio::DB::Biblio::pubmed_dbi;,..., args
sub TIEHASH {
  my $class = shift;
  my $self=$class->new(@_);
  $self;
}  
sub FETCH {
  my $self=shift @_;
  $self->get_by_id(@_);
}
sub STORE {
  my($self,$id,$object)=@_;
  $self->throw("Attempt to store into collection opened for read-only") if $self->read_only;
  $self->throw("Attempt to store object under wrong ID") unless $id eq $object->pmid;
  $self->put($object);
}
sub DELETE {
  my $self=shift @_;
  $self->throw("Attempt to delete from collection opened for read-only") if $self->read_only;
  $self->del(@_);
}
sub CLEAR {
  my $self=shift @_;
  $self->throw("Attempt to clear collection opened for read-only") if $self->read_only;
  $self->del_all;
}
sub EXISTS {
  my($self,$id)=@_;
  # dumb implementation -- should implement 'exists' method on object
  $self->get_by_id($id)? 1: 0;
}
sub FIRSTKEY {
  my $self=shift @_;
  $self->reset;
  my $old_what=$self->what;
  $self->what('ids');
  my $id=$self->get_next;
  $self->what($old_what);
  $id;
}  
sub NEXTKEY {
  my $self=shift @_;
  my $old_what=$self->what;
  $self->what('ids');
  my $id=$self->get_next;
  $self->what($old_what);
  $id;
}  

1;
