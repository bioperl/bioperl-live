package Bio::DB::Biblio::pubmed_gdbm;
use GDBM_File;
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
  @AUTO_ATTRIBUTES=qw(access dbm_module file read_only create
		     _db _iteration_end);
  @OTHER_ATTRIBUTES=qw(what);
  %SYNONYMS=(get_collection_id=>'collection_id');
  Bio::Root::AutoClass::declare(__PACKAGE__,\@AUTO_ATTRIBUTES,\%SYNONYMS);
}

sub _init_self {
  my($self,$class,$args)=@_;
  return unless $class eq __PACKAGE__; # to prevent subclasses from re-running this
  my %db;
  tie(%db,'GDBM_File',$self->file,$self->gdbm_flags,$self->create_mode) || 
    $self->throw("Cannot tie GDBM file ".$self->file."\n");
  $self->_db(\%db);
  $self->dbm_module('GDBM_File');
}
sub create_mode {
  my $self=shift @_;
  $self->{'CREATE_MODE'}=$_[0] if @_;
  defined $self->{'CREATE_MODE'}? $self->{'CREATE_MODE'}: $CREATE_MODE;
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
  $self->{'COLLECTION_ID'}=$_[0] if @_;
  defined $self->{'COLLECTION_ID'}? $self->{'COLLECTION_ID'}: $self->file;
}
sub get_collection_id {$_[0]->collection_id;}

sub count {
  my $self=shift @_;
  return $self->{'COUNT'}=$_[0] if @_;
  return $self->{'COUNT'} if defined $self->{'COUNT'};
  # else have to calculate it...
  my $db=$self->_db;
  $self->{'COUNT'}=scalar keys %$db;
}

sub put {
  my($self,@objects)=@_;
  $self->throw("put called on collection opened for read-only") if $self->read_only;
  my $db=$self->_db;
  for my $object (@objects) {
    my $id=$object->pmid;
    next unless defined $id;
    $db->{"$id"}=freeze($object);
  }
  $self->count(undef)		# cached version of count no longer valid
}
sub del {
  my($self,@args)=@_;
  $self->throw("del called on collection opened for read-only") if $self->read_only;
  my $db=$self->_db;
  for my $arg (@args) {
    my $id=!ref($arg)? $arg: $arg->can('pmid')? $arg->pmid: $arg->can('id')? $arg->id: undef;
    next unless defined $id;
    next unless exists $db->{"$id"};
    delete $db->{"$id"};
  }
  $self->count(undef);		# cached version of count no longer valid
}
sub del_all {
  my($self,@args)=@_;
  $self->throw("del_all called on collection opened for read-only") if $self->read_only;
  my $db=$self->_db;
  %$db=();
  $self->count(undef);
}
sub get_all {
  my $self=shift @_;
  my $what=@_? _what($_[0]) || $self->throw("Invalid what $_[0]"): $self->what;
  my $db=$self->_db;
  scalar keys %$db;		# reset built-in iterator
  my @results;
  while (my($id,$freeze)=each %$db) {
    push(@results,thaw($freeze)) if $what eq 'refs';
    push(@results,$id) if $what eq 'ids';
  }
  @results;
}
sub get_by_id {
  my($self,@ids)=@_;
  my $what=$self->what;
  my $db=$self->_db;
  my @results;
  for my $id (@ids) {
    next unless exists $db->{"$id"};
    my $freeze=$db->{"$id"};
    push(@results,thaw($freeze)) if $what eq 'refs';
    push(@results,$id) if $what eq 'ids';
  }
  wantarray? @results: $results[0];
}
sub get_by_ids {my $self=shift @_; $self->get_by_id(@_);}

sub get_next {
  my $self=shift @_;
  my $n=@_? $_[0]: 1;
  my $db=$self->_db;
  my $what=$self->what;
  my @results;
  unless ($self->_iteration_end) {
    while ($n-->0) {
      my @entry=each %$db;
      $self->_iteration_end(1), last unless @entry;
      my($id,$freeze)=@entry;
      push(@results,thaw($freeze)) if $what eq 'refs';
      push(@results,$id) if $what eq 'ids';
    }
  }
  wantarray? @results: $results[0];
}
sub reset {
  my $self=shift @_;
  @_ and $self->what($_[0]);
  my $db=$self->_db;
  $self->_iteration_end(0);
  scalar keys %$db;		# reset built-in iterator
  $self;
}
sub get_more {my $self=shift @_; $self->get_next(@_);}
sub reset_retrieval {my $self=shift @_; $self->reset(@_);}

sub gdbm_flags {
  my($self)=@_;
  my $create_always=$self->create;
  my $create_sometimes=!defined $self->create;
  return GDBM_READER if $self->read_only;
  return GDBM_WRITER unless $self->read_only || $create_always || $create_sometimes;
  return GDBM_NEWDB if $create_always;
  return GDBM_WRCREAT if $create_sometimes;
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
  my $db=$self->{'biblio'}->_db;
  untie %$db;
  $self->{'biblio'}->_db(undef);
}
sub UNTIE {}

# TIEHASH interface
# tie %h,Bio::DB::Biblio::pubmed_gdbm,..., args
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
