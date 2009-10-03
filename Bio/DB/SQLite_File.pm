# $Id$
#
# BioPerl module for Bio::DB::SQLite_File
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Mark A. Jensen <maj@fortinbras.us>
#
# Copyright Mark A. Jensen
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::SQLite_File - A minimal DBM using SQLite

=head1 SYNOPSIS

 # tie a simple hash to a SQLite DB file
 my %db;
 tie(%db, 'Bio::DB::SQLite_File', 'my.db');

 # tie an array
 my @db;
 tie(@db, 'Bio::DB::SQLite_File', 'my.db');

 # tie to a tempfile
 tie(%db, 'Bio::DB::SQLite_File', undef);

 # get attributes of the tied object

 $SQLite_handle = (tied %db)->dbh;
 $db_file = (tied %db)->file;

=head2 SQL interface 
 # use as an option in AnyDBM_File
 @AnyDBM_File::ISA = qw( DB_File Bio::DB::SQLite_File SDBM );
 my %db;
 tie(%db, 'AnyDBM_File', 'my.db', $flags, $mode, @dbmargs)
 
=head1 DESCRIPTION

This module allows a hash or an array to be tied to a SQLite DB via
L<DBI> plus L<DBD::SQLite>, in a way that emulates some features of
DB_File. The class is suitable for L<Bio::DB::FileCache>, for which it
was written. In particular, this module offers another choice for
ActiveState users, who may find it difficult to get a working
L<DB_File> (based on Berkeley DB) installed, but can't failover to
SDBM due to its record length restrictions. Bio::DB::SQLite_File
requires L<DBD::SQLite>, which has SQLite built-in -- no external
application install required.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support

Please direct usage questions or support issues to the mailing list:

L<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and
reponsive experts will be able look at the problem and quickly
address it. Please include a thorough description of the problem
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Mark A. Jensen

Email maj -at- fortinbras -dot- us

=head1 CONTRIBUTORS

This code owes an intellectual debt to Lincoln Stein. Inelegancies and
bugs are mine.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# for pod discussion
#
# R_CURSOR points to a record, not between records
#
# interaction between the two APIs
# duplicate protection and handling


# Let the code begin...

package Bio::DB::SQLite_File;
use strict;
use warnings;

BEGIN {
    unless (eval "require DBD::SQLite; 1") {
 	Bio::Root::Root->throw( "SQLite_File requires DBD::SQLite" );
     }
}

our @EXPORT = qw( 
                 $DB_HASH $DB_BTREE $DB_RECNO 
                 R_DUP R_CURSOR R_FIRST R_LAST 
                 R_NEXT R_PREV R_IAFTER R_IBEFORE 
                 R_NOOVERWRITE R_SETCURSOR
                 O_CREAT O_RDWR O_RDONLY O_SVWST
                 );

our $DB_HASH = { 'type' => 'HASH' };
our $DB_BTREE = { 'type' => 'BINARY' };
our $DB_RECNO = { 'type' => 'RECNO' };

# constants hacked out of DB_File:
sub R_DUP { 32678 }
sub R_CURSOR { 27 }
sub R_FIRST { 7 }
sub R_LAST { 15 }
sub R_NEXT { 16 }
sub R_PREV { 23 }
sub R_IAFTER { 1 }
sub R_IBEFORE { 3 }
sub R_NOOVERWRITE { 20 }
sub R_SETCURSOR { -100 }

sub O_CREAT { 512 };
sub O_RDWR { 2 };
sub O_RDONLY  { 0 };
sub O_SVWST { 514 };

# exporter...
# AnyDBM_File performs a require, but Exporter only exports at compile time.
# so we have to kludge it...
# my ($pkg, $fn, $ln) = caller;

# no warnings; # avoid fscking "uninit" warns
# for (@EXPORT) {
#     m/^\$(.*)/ && do {
# 	eval "\$${pkg}::$1 = \$Bio::DB::SQLite_File::$1";
#     };
#     m/^[^\$@%]/ && do {
# 	eval "\*\{${pkg}::$_\} = \\&Bio::DB::SQLite_File::$_";
#     };
#     print STDERR $@ if $@;
# }
# use warnings;

# Object preamble - inherits from Bio::Root::Root

# my notes /maj
# FileCache.pm does the following using DB_File
# $cache object first looks for the sequence in the cache (with _get)
#   if found, returns the sequence
#   if not found, finds the sequence in the database, 
#    and inserts it in the cache (with _store), then returns it
#
# the _store method: freezes the object with Storable, then 
# puts it in the cacheDB by using the tied hash
# very simple: just ${type}_${id} => $frozen_objects
# 
# the _get method:
# gets the frozen object via the tied hash ($self->db->{"${type}_${id}"})
# thaws it out with Storable, and returns the real object.

# SO... FileCache.pm can be used almost verbatim, but the tied hash
# representing the cacheDB implemented in SQLite needs to be built.

# SO..this module provides the tied SQLite hash class.

# strange bug using t/data/taxdump/names.db
#
# with r16188, getting a crash in transfac_pro.t with "database locked" 
# error - couldn't debug or trace to it. This was a hash element set.
# 
# offending records were of name 'rhodotorula', about 8-10 entries
# with this name
# Removing all 'rhodotoruala' records beyond 3 allows the test to proceed.
#
# My intuition is that the sqlite hashing alg couldn't create a
# unique hashkey for all these entries.
#
# My fix will be to add an autoincrementing primary key to the table.
# That didn't work.
# 


use Bio::Root::Root;
use DBI;
use File::Temp qw( tempfile );
use Fcntl qw(O_CREAT O_RDWR O_RDONLY);
use base qw( Bio::Root::Root );

# these globals probably need to be associated with the object
# for multiple objects to be present.

our $AUTOKEY = 0;
# for providing DB_File seq functionality
our $AUTOPK = 0;

sub SEQIDX {
    my $self = shift;
    return $self->{SEQIDX} = [] if (!defined $self->{SEQIDX});
    return $self->{SEQIDX};
}

sub CURSOR {
    my $self = shift;
    $self->{CURSOR} = 0 if (!defined $self->{CURSOR});
    return \$self->{CURSOR};
}

sub TIEHASH {
    my $class = shift;
    my ($file, $flags, $mode, $index, $keep) = @_;
    my $self = {};
    bless($self, $class);
    $self->{ref} = 'HASH';
    my $fh;
    # db file handling
    if ($file) {
	# you'll love this...
	my $infix;
	my $setmode;
	for ($flags) {
	    !defined && do {
		$infix = ">";
		last;
	    };
	    $_ eq 'O_SVWST' && do { #bullsith kludge
		$_ = 514;
	    };
	    ($_ & O_CREAT) && do {
		$setmode = 1 if ! -e $file;
		$infix = (-e $file ? '<' : '>');
	    };
	    ($_ & O_RDWR) && do {
		$infix = '+'.($infix ? $infix : '<');
	    };
	    do { # O_RDONLY
		$infix = '<' unless $infix;
	    };
	}
	open($fh, $infix, $file) or $self->throw("Can't open db file: $!");
	chmod $mode, $file if $setmode;
	# if file explicitly specified, but keep is not, 
	# retain file at destroy...
	$keep = 1 if !defined $keep;
    }
    else {
	# if no file specified, use a temp file...
	($fh, $file) = tempfile();
	# if keep not explicitly specified, 
	# remove the tempfile at destroy...
	$keep = 0 if !defined $keep;
    }
    $self->file($file);
    $self->_fh($fh);
    $self->keep($keep);

    # create SQL statements
     my $hash_tbl = <<END;
    (
      id      blob,
      obj     blob not null,
      pk      integer primary key autoincrement
    );
END
    my $create_idx = <<END;
    CREATE INDEX IF NOT EXISTS id_idx ON hash ( id, pk );
END
    my $dbh = DBI->connect("dbi:SQLite:dbname=".$self->file,"","",
			   {RaiseError => 1, AutoCommit => 1});
    $self->dbh( $dbh );
    # rudimentary DB_File emulation:
    if (defined $index) {
	$self->throw("Index selector must be a hashref") unless ref($index) eq 'HASH';
	my $flags = $index->{flags} || 0;
	for ($index->{'type'}) {
	    !defined && do {
		$self->dbh->do("CREATE TABLE IF NOT EXISTS hash $hash_tbl");
		last;
	    };
	    $_ eq 'BINARY' && do {
		if ($flags & R_DUP ) {
		    $self->dup(1);
		    $self->dbh->do("CREATE TABLE IF NOT EXISTS hash $hash_tbl");
		    $self->dbh->do($create_idx);
		}
		else {
		    $self->dup(0);
		    $self->dbh->do("CREATE TABLE IF NOT EXISTS hash $hash_tbl");
		    $self->dbh->do($create_idx);
		}
	    };
	    $_ eq 'HASH' && do {
		$self->dbh->do("CREATE TABLE IF NOT EXISTS hash $hash_tbl");
		last;
	    };
	    $_ eq 'RECNO' && do {
		$self->throw("$DB_RECNO is not meaningful for tied hashes");
		last;
	    };

	}
    }
    else {
	$self->dbh->do("CREATE TABLE IF NOT EXISTS hash $hash_tbl");
    }
    return $self;
}

sub TIEARRAY {
    my $class = shift;
    my ($file, $flags, $mode, $index, $keep) = @_;
    my $self = {};
    bless($self, $class);
    $self->{ref} = 'ARRAY';
    my $fh;
    # db file handling
    if ($file) {
	my $infix;
	my $setmode;
	for ($flags) {
	    !defined && do {
		$infix = ">";
		last;
	    };
	    $_ eq 'O_SVWST' && do { #bullsith kludge
		$_ = 514;
	    };
	    ($_ & O_CREAT) && do {
		$setmode = 1 if ! -e $file;
		$infix = (-e $file ? '<' : '>');
	    };
	    ($_ & O_RDWR) && do {
		$infix = '+'.($infix ? $infix : '<');
	    };
	    do { # O_RDONLY
		$infix = '<' unless $infix;
	    };
	}
	open($fh, $infix, $file) or $self->throw("Can't open db file: $!");
	chmod $mode, $file if $setmode;
	# if file explicitly specified, but keep is not, 
	# retain file at destroy...
	$keep = 1 if !defined $keep;
    }
    else {
	# if no file specified, use a temp file...
	($fh, $file) = tempfile();
	# if keep not explicitly specified, 
	# remove the tempfile at destroy...
	$keep = 0 if !defined $keep;
    }
    $self->file($file);
    $self->_fh($fh);
    $self->keep($keep);
    
    my $arr_tbl = <<END;
    (
      id      integer primary key,
      obj     blob not null
    );
END
    
    my $create_idx = <<END;
    CREATE INDEX IF NOT EXISTS id_idx ON hash ( id );
END
    
    my $dbh = DBI->connect("dbi:SQLite:dbname=".$self->file,"","",
			   {RaiseError => 1, AutoCommit => 1});
    $self->dbh( $dbh );
    # rudimentary DB_File emulation:
    if (defined $index) {
	$self->throw("Index selector must be a hashref") unless ref($index) eq 'HASH';
	for ($index->{'type'}) {
	    $_ eq 'BINARY' && do {
		$self->dbh->disconnect;
		$self->throw("$DB_BTREE is not meaningful for a tied array");
		last;
	    };
	    $_ eq 'HASH' && do {
		$self->dbh->disconnect;
		$self->throw("$DB_HASH is not meaningful for a tied array");
		last;
	    };
	    $_ eq 'RECNO' && do {
		$self->dbh->do("CREATE TABLE IF NOT EXISTS hash $arr_tbl");
		$self->dbh->do($create_idx);
	    };
	    !defined && do {
		$self->dbh->do("CREATE TABLE IF NOT EXISTS hash $arr_tbl");
		last;
	    };
	}
    }
    else {
	$self->dbh->do("CREATE TABLE IF NOT EXISTS hash $arr_tbl");
	$self->dbh->do($create_idx);
    }

    $self->{len} = 0;
    return $self;
}



# common methods

sub FETCH {
    my $self = shift;
    my $key = shift;
    return unless $self->dbh;
    if (!$self->{ref} or $self->ref eq 'HASH') {
	$self->get_sth->execute($key); # fetches on column 'id'
    }
    elsif ($self->ref eq 'ARRAY') {
	if (defined ${$self->SEQIDX}[$key]) {
	    $self->get_sth->execute($self->get_idx($key));
	}
	else {
	    $self->_last_pk(undef);
	    return undef;
	}
    }
    else { # type not recognized
        $self->throw("tied type not recognized");
    }
    my $ret = $self->get_sth->fetch;
    if ($ret) {
	$self->_last_pk( $ret->[1] ); # store the returned pk
	$ret->[0] =~ s{<SQUOT>}{'}g;
	$ret->[0] =~ s{<DQUOT>}{"}g;
	return $ret->[0]; # always returns the object
    }
    else {
	$self->_last_pk( undef ); # fail in pk
	return $ret;
    }
}

sub STORE {
    my $self = shift;
    my ($key, $value) = @_;
    return unless $self->dbh;
    $value =~ s{'}{<SQUOT>}g;
    $value =~ s{"}{<DQUOT>}g;
    if ( !defined $self->{ref} or $self->ref eq 'HASH' ) {
	if ( $self->dup ) { # allowing duplicates
	    my $pk = $self->_get_pk;
	    $self->put_sth->execute($key, $value, $pk);
	    push @{$self->SEQIDX}, $pk;
	}
	else { # no duplicates...
	    #need to check if key is already present
	    if ( $self->EXISTS($key) ) {
		$self->upd_sth->execute( $value, $key, $self->_last_pk )
	    }
	    else {
		my $pk = $self->_get_pk;
		$self->put_sth->execute( $key, $value, $pk );
		push @{$self->SEQIDX}, $pk;
	    }
	}
    }
    elsif ( $self->ref eq 'ARRAY' ) {
	# need to check if the key is already present
	if (!defined ${$self->SEQIDX}[$key] ) {
	    $self->put_sth->execute($self->get_idx($key), $value);
	}
	else {
	    $self->upd_sth->execute($value,$self->get_idx($key));
	}
    }
    $value;
}

sub DELETE {
    my $self = shift;
    my $key = shift;
    return unless $self->dbh;
    my $oldval;
    if (!$self->ref or $self->ref eq 'HASH') {
	return unless $self->get_sth->execute($key);
	my $ret = $self->get_sth->fetch;
	$oldval = $ret->[0];
	$self->dbh->do("DELETE FROM hash WHERE id = '$key'");
	# update the sequential side
	for (@{$self->SEQIDX}) { # not very efficient.
	    if ($_ == $ret->[1]) {
		undef $_;
		last;
	    }
	}
    }
    elsif ($self->ref eq 'ARRAY') {
	my $SEQIDX = $self->SEQIDX;
	if ($$SEQIDX[$key]) {
	    $oldval = $self->FETCH($$SEQIDX[$key]);
	    $self->dbh->do("DELETE FROM hash WHERE id = '$$SEQIDX[$key]'");
	    $self->rm_idx($key);
	}
    }
    else {
	$self->throw( "tied type not recognized" );
    }
    return $oldval;
}

sub EXISTS {
    my $self = shift;
    my $key = shift;
    return unless $self->dbh;
    if (!$self->ref or $self->ref eq 'HASH') {
	return defined $self->FETCH($key) ? 1 : 0;
    }
    elsif ($self->ref eq 'ARRAY') {
	if (defined(${$self->SEQIDX}[$key])) {
	    $self->_last_pk(${$self->SEQIDX}[$key]);
	    return 1;
	}
	else {
	    $self->_last_pk(undef);
	    return 0;
	}
    }
    else {
	$self->throw("tied type not recognized");
    }
}

sub CLEAR {
    my $self = shift;
    return unless $self->dbh;
    $self->dbh->do("DELETE FROM hash");
    @{$self->SEQIDX} = ();
    return 1;
}

# hash methods

sub FIRSTKEY {
    my $self = shift;
    return unless $self->dbh;
    return if ($self->{ref} and $self->ref ne 'HASH');
    my $ids = $self->dbh->selectall_arrayref("SELECT id FROM hash");
    return unless $ids;
    return $self->_keys($ids);
}

sub NEXTKEY {
    my $self = shift;
    my $lastkey = shift;
    return unless $self->dbh;
    return if ($self->{ref} and $self->ref ne 'HASH');
    return $self->_keys;
}

# array methods

sub FETCHSIZE {
    my $self = shift;
    return unless $self->dbh;
    return if (!$self->{ref} or $self->ref ne 'ARRAY');
    $self->len;
}

sub STORESIZE {
    my $self = shift;
    my $count = shift;
    return unless $self->dbh;
    return if (!$self->ref or $self->ref ne 'ARRAY');
    if ($count > $self->len) {
	foreach ($count - $self->len .. $count) {
	    $self->STORE($_, '');
	}
    }
    elsif ($count < $self->len) {
	foreach (0 .. $self->len - $count - 2) {
	    $self->POP();
	}
    }
}

# EXTEND is no-op
sub EXTEND {
    my $self = shift;
    my $count = shift;
    return;
}

sub POP {
    my $self = shift;
    return unless $self->dbh;
    return if (!$self->{ref} or $self->ref ne 'ARRAY');
    $self->get_sth->execute($self->get_idx($self->len-1));
    my $ret = $self->get_sth->fetch;
    $self->dbh->do("DELETE FROM hash WHERE id = ".$self->get_idx($self->len-1));
    # bookkeeping
    $self->rm_idx($self->len-1);
    return defined $ret ? $ret->[0] : $ret;
}

sub PUSH {
    my $self = shift;
    my @values = @_;
    return unless $self->dbh;
    return if (!$self->{ref} or $self->ref ne 'ARRAY');
    my $ret = @values;
    my $beg = $self->len;
    my $end = $self->len + @values - 1;
    for my $i ($beg..$end) {
	$self->put_sth->execute($self->get_idx($i), shift @values);
    }
    return $ret;
}

sub SHIFT {
    my $self = shift;
    return unless $self->dbh;
    return if (!$self->{ref} or $self->ref ne 'ARRAY');
    $self->get_sth->execute( $self->get_idx(0) );
    my $ret = $self->get_sth->fetch;
    $self->dbh->do("DELETE FROM hash WHERE id = ".$self->get_idx(0));
    # bookkeeping
    $self->shift_idx;
    return defined $ret ? $ret->[0] : $ret;
}

sub UNSHIFT {
    my $self = shift;
    my @values = @_;
    return if (!$self->{ref} or $self->ref ne 'ARRAY');
    my $n = @values;
    return unless $self->dbh;
    for ($self->unshift_idx($n)) {
	$self->put_sth->execute($_,shift @values);
    }
    return $n;
}

sub SPLICE {
    my $self = shift;
    return if !$self->dbh;
    return if (!$self->{ref} or $self->ref ne 'ARRAY');
    die "Won't do splice yet."
}

# destructors

sub UNTIE {
    my $self = shift;
    my $count = shift;
    $self->throw( "untie attempted while $count inner references still exist" ) if ($count);}

sub DESTROY {
    my $self = shift;
    # finish and destroy stmt handles
    for ( qw( put_sth put_seq_sth get_sth get_seq_sth upd_sth upd_seq_sth ) ) {
	$self->$_->finish if $self->{$_};
	undef $self->{$_};
    }
    # disconnect
    $self->dbh->disconnect;
    # remove file if nec
    $self->_fh->close() if $self->_fh;
    unlink $self->file if (!$self->keep && $self->_fh);
    $self->throw("SQLite_File unlink issue: $!") if $!;
    undef $self;
    1;
}


=head2 SQL interface

=head2 dbh

 Title   : dbh
 Usage   : $obj->dbh($newval)
 Function: database handle
 Example : 
 Returns : value of dbh (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub dbh {
    my $self = shift;
    
    return $self->{'dbh'} = shift if @_;
    return $self->{'dbh'};
}

=head2 get_sth

 Title   : get_sth
 Usage   : $obj->get_sth($new_handle)
 Function: returns a select statement handle
           prepares a new one if necessary
 Example : 
 Returns : value of get_sth (a statement handle)
 Args    : none

=cut

sub get_sth {
    my $self = shift;
    $self->throw ("No active database handle") unless $self->dbh;
    # prepare a new stmt if dne or if the column requested is different...
    if (!$self->{'get_sth'}) {
	if (!$self->{ref} or $self->ref eq 'HASH') {
	    $self->{'get_sth'} =  $self->dbh->prepare("SELECT obj, pk FROM hash WHERE id = ?");
	}
	elsif ($self->ref eq 'ARRAY') {
	    $self->{'get_sth'} =  $self->dbh->prepare("SELECT obj, id FROM hash WHERE id = ?");
	}
	else {
	    $self->throw("tied type not recognized");
	}
    }
    return $self->{'get_sth'};
}


=head2 put_sth

 Title   : put_sth
 Usage   : $obj->put_sth($new_handle)
 Function: returns an insert statement handle
           prepares a new one if necessary
 Example : 
 Returns : value of put_sth (a statement handle)
 Args    : none

=cut

sub put_sth {
    my $self = shift;
    $self->throw ("No active database handle") unless $self->dbh;
    if (!$self->{'put_sth'}) {
	if (!$self->{ref} or $self->ref eq 'HASH') {
	    $self->{'put_sth'} = $self->dbh->prepare("INSERT INTO hash (id, obj, pk) VALUES ( ?, ?, ? )");
	}
	elsif ($self->{ref} eq 'ARRAY') {
	    $self->{'put_sth'} = $self->dbh->prepare("INSERT INTO hash (id, obj) VALUES ( ?, ?)");
	}
	else {
	    $self->throw("tied type not recognized");
	}
    }
    return $self->{'put_sth'};
}

=head2 upd_sth

 Title   : upd_sth
 Usage   : $obj->upd_sth($new_handle)
 Function: returns an update statement handle
           prepares a new one if necessary
 Example : 
 Returns : value of upd_sth (a statement handle)
 Args    : none

=cut

sub upd_sth {
    my $self = shift;
    $self->throw ("No active database handle") unless $self->dbh;
    if (!$self->{'upd_sth'}) {
	if (!$self->{ref} or $self->ref eq 'HASH') {
	    $self->{'upd_sth'} = $self->dbh->prepare("UPDATE hash SET obj = ? WHERE id = ? AND pk = ?");
	}
	elsif ($self->ref eq 'ARRAY') {
	    $self->{'upd_sth'} = $self->dbh->prepare("UPDATE hash SET obj = ? WHERE id = ?");
	}
	else {
	    $self->throw("tied type not recognized");
	}
    }
    return $self->{'upd_sth'};
}

=head2 get_seq_sth

 Title   : get_seq_sth
 Usage   : $obj->get_seq_sth($new_handle)
 Function: returns a select statement handle
           prepares a new one if necessary
 Example : 
 Returns : value of get_seq_sth (a statement handle)
 Args    : none

=cut

sub get_seq_sth {
    my $self = shift;
    $self->throw ("No active database handle") unless $self->dbh;
    # prepare a new stmt if dne or if the column requested is different...
    if (!$self->{'get_seq_sth'}) {
	if (!$self->{ref} or $self->ref eq 'HASH') {
	    $self->{'get_seq_sth'} =  $self->dbh->prepare("SELECT id, obj FROM hash WHERE pk = ?");
	}
	elsif ($self->ref eq 'ARRAY') {
	    $self->{'get_seq_sth'} =  $self->dbh->prepare("SELECT id, obj FROM hash WHERE id = ?");
	}
	else {
	    $self->throw("tied type not recognized");
	}
    }
    return $self->{'get_seq_sth'};
}


=head2 put_seq_sth

 Title   : put_seq_sth
 Usage   : $obj->put_seq_sth($new_handle)
 Function: returns an insert statement handle
           prepares a new one if necessary
 Example : 
 Returns : value of put_seq_sth (a statement handle)
 Args    : none

=cut

sub put_seq_sth {
    my $self = shift;
    $self->throw ("No active database handle") unless $self->dbh;
    if (!$self->{'put_seq_sth'}) {
	if (!$self->{ref} or $self->ref eq 'HASH') {
	    $self->{'put_seq_sth'} = $self->dbh->prepare("INSERT INTO hash (id, obj, pk) VALUES ( ?, ?, ? )");
	}
	elsif ($self->{ref} eq 'ARRAY') {
	    $self->{'put_seq_sth'} = $self->dbh->prepare("INSERT INTO hash (obj, id) VALUES ( ?, ?)");
	}
	else {
	    $self->throw("tied type not recognized");
	}
    }
    return $self->{'put_seq_sth'};
}

=head2 upd_seq_sth

 Title   : upd_seq_sth
 Usage   : $obj->upd_seq_sth($new_handle)
 Function: returns an update statement handle
           prepares a new one if necessary
 Example : 
 Returns : value of upd_seq_sth (a statement handle)
 Args    : none

=cut

sub upd_seq_sth {
    my $self = shift;
    $self->throw ("No active database handle") unless $self->dbh;
    if (!$self->{'upd_seq_sth'}) {
	if (!$self->{ref} or $self->ref eq 'HASH') {
	    $self->{'upd_seq_sth'} = $self->dbh->prepare("UPDATE hash SET id = ?, obj = ? WHERE pk = ?");
	}
	elsif ($self->ref eq 'ARRAY') {
	    $self->{'upd_seq_sth'} = $self->dbh->prepare("UPDATE hash SET obj = ? WHERE id = ?");
	}
	else {
	    $self->throw("tied type not recognized");
	}
    }
    return $self->{'upd_seq_sth'};
}

=head2 Attribute accessors

=head2 file

 Title   : file
 Usage   : $obj->file($newval)
 Function: filename for the SQLite db
 Example : 
 Returns : value of file (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub file {
    my $self = shift;
    
    return $self->{'file'} = shift if @_;
    return $self->{'file'};
}

=head2 _fh

 Title   : _fh
 Usage   : $obj->_fh($newval)
 Function: holds the temp file handle
 Example : 
 Returns : value of _fh (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub _fh {
    my $self = shift;
    
    return $self->{'_fh'} = shift if @_;
    return $self->{'_fh'};
}

=head2 keep

 Title   : keep
 Usage   : $obj->keep($newval)
 Function: flag allows preservation of db file when set
 Example : 
 Returns : value of keep (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub keep {
    my $self = shift;
    
    return $self->{'keep'} = shift if @_;
    return $self->{'keep'};
}

# HASH or ARRAY?
sub ref {
    my $self = shift;
    return $self->{ref};
}

=head2 _keys

 Title   : _keys
 Usage   : internal
 Function: points to a hash to make iterating easy and fun
 Example : 
 Returns : value of _keys (a hashref)
 Args    : on set, an arrayref of scalar keys

=cut

sub _keys {
    my $self = shift;
    my $load = shift;
    if ($load) {
	$self->{'_keys'} = {};
	@{$self->{'_keys'}}{ map {$_->[0]} @$load } = (undef) x @$load;
	my $a = keys %{$self->{'_keys'}}; #reset each
    }
    return each %{$self->{'_keys'}};
}

=head2 Array object helper functions

=cut

sub inc {
    return ++(shift->{len});
}

sub dec {
    my $self = shift;
    return $self->{len} ? --($self->{len}) : 0;
}

sub len {
    scalar @{shift->SEQIDX};
}

sub get_idx {
    my $self = shift;
    my $index = shift;
    my $SEQIDX = $self->SEQIDX;
    return $$SEQIDX[$index] if defined $$SEQIDX[$index];
    push @$SEQIDX, $AUTOKEY;
    $$SEQIDX[$index] = $AUTOKEY++;
}

sub shift_idx {
    my $self = shift;
    return shift( @{$self->SEQIDX} );
}

# returns the set of new db ids to use
sub unshift_idx {
    my $self = shift;
    my $n = shift;
    my @new;
    push(@new, $AUTOKEY++) for (0..$n-1);
    unshift @{$self->SEQIDX}, @new;
    return @new;
}

sub rm_idx {
    my $self = shift;
    my $index = shift;
    unless (delete ${$self->SEQIDX}[$index]) {
	warn "Element $index did not exist";
    }
}
=head2 BDB API Emulation: random access

=head2 get()

 Title   : get
 Usage   : $X->get($key, $value, $flags)
 Function: see DB_File
 Returns : 
 Args    : 

=cut

sub get {
    my $self = shift;
    my ($key, $value, $flags) = @_;
    return unless $self->dbh;
    $_[1] = ($self->ref eq 'ARRAY' ? $self->FETCH(${$self->SEQIDX}[$key]) : $self->FETCH($key));
    return 0 if defined $_[1];
    return 1;
}

=head2 put()

 Title   : put
 Usage   : $X->put($key, $value, $flags)
 Function: see DB_File
 Returns : 
 Args    : 

=cut

sub put {
    my $self = shift;
    my ($key, $value, $flags) = @_;
    return unless $self->dbh;
    my $SEQIDX = $self->SEQIDX;
    my $CURSOR = $self->CURSOR;
    my ($status, $pk, @parms);
    for ($flags) {
	(!defined || $_ == R_SETCURSOR) && do { # put or upd
	    if ($self->dup) { # just make a new one
		$pk = $self->_get_pk;
		@parms = ($self->ref eq 'ARRAY' ? ($value, $pk) : ($key, $value, $pk));
		$status = !$self->put_seq_sth->execute(@parms);
		unless ($status) {
		    push @$SEQIDX, $pk;
		    $$CURSOR = $#$SEQIDX if $_ == R_SETCURSOR;
		}
	    }
	    else {

		$self->FETCH($key);
		$pk = $self->_last_pk || $self->_get_pk;
		@parms = ($self->ref eq 'ARRAY' ? ($value, $pk) : ($key, $value, $pk));
		if ($self->_last_pk) {
		    $status = !$self->upd_seq_sth->execute(@parms);
		}
		else {
		    $status = !$self->put_seq_sth->execute(@parms);
		    push @$SEQIDX, $pk;
		}
		$_ && do { # R_SETCURSOR
		    if ( $pk = $$SEQIDX[-1] ) {
			$$CURSOR = $#$SEQIDX;
		    }
		    else {
			for (my $i=0; $i<@$SEQIDX; $i++) {
			    $$CURSOR = $i;
			    last if $$SEQIDX[$i] == $pk;
			}
		    };
		}
	    }
	    last;
	};
	$_ == R_IAFTER && do {
	    $self->_wring_SEQIDX unless $$SEQIDX[$$CURSOR];
	    # duplicate protect
	    return 1 unless ($self->ref eq 'ARRAY') || $self->dup || !$self->EXISTS($key);
	    $pk = $self->_get_pk;
	    if ($$CURSOR == $#$SEQIDX) {
		push @$SEQIDX, $pk;
	    }
	    else {
		splice(@$SEQIDX,$$CURSOR,0,$pk);
	    }
	    @parms = ($self->ref eq 'ARRAY' ? ($value, $pk) : ($key, $value, $pk));
	    $status = !$self->put_seq_sth->execute(@parms);
	    $_[0] = $$CURSOR+1 if !$status;
	    last;
	};
	$_ == R_IBEFORE && do {
	    $self->_wring_SEQIDX unless $$SEQIDX[$$CURSOR];
	    # duplicate protect
	    return 1 unless ($self->ref eq 'ARRAY') || $self->dup || !$self->EXISTS($key);
	    $pk = $self->_get_pk;
	    if ($$CURSOR) {
		splice(@$SEQIDX,$$CURSOR-1,0,$pk);
	    }
	    else {
		unshift(@$SEQIDX, $pk);
	    }
	    @parms = ($self->ref eq 'ARRAY' ? ($value, $pk) : ($key, $value, $pk));
	    $status = !$self->put_seq_sth->execute(@parms);
	    $_[0] = $$CURSOR if !$status;
	    $$CURSOR++; # preserve cursor
	    last;
	};
	$_ == R_CURSOR && do { # upd only
	    $self->_wring_SEQIDX unless $$SEQIDX[$$CURSOR];
	    # duplicate protect
	    return 1 unless ($self->ref eq 'ARRAY') || $self->dup || !$self->EXISTS($key);
	    $pk = $$SEQIDX[$$CURSOR];
	    @parms = ($self->ref eq 'ARRAY' ? ($value, $pk) : ($key, $value, $pk));
	    $status = !$self->upd_seq_sth->execute(@parms);
	    last;
	};
	$_ == R_NOOVERWRITE && do { # put only/add to the "end"
	    #will create a duplicate if $self->dup is set!
	    return 1 unless ($self->ref eq 'ARRAY') || $self->dup || !$self->EXISTS($key);
	    $pk = $self->_get_pk;
	    push @$SEQIDX, $pk;
	    @parms = ($self->ref eq 'ARRAY' ? ($value, $pk) : ($key, $value, $pk));
	    $self->put_seq_sth->execute(@parms);
	    last;
	};

    }
    use warnings;
    return $status;
}

=head2 del()

 Title   : del
 Usage   : 
 Function: as in DB_file
 Returns : 
 Args    : 

=cut

sub del {
    my $self = shift;
    my ($key, $flags) = @_;
    return unless $self->dbh;
    my $SEQIDX = $self->SEQIDX;
    my $CURSOR = $self->CURSOR;
    my $status;
    if ($flags eq R_CURSOR) {
	_wring_SEQIDX($self->SEQIDX) unless $$SEQIDX[$$CURSOR];
	my $pk = $$SEQIDX[$$CURSOR];
	my $col = ($self->ref eq 'ARRAY' ? 'id' : 'pk');
	$status = $self->dbh->do("DELETE FROM hash WHERE $col = $pk");
	if ($status) { # successful delete
	    $$SEQIDX[$$CURSOR] = undef;
	    $self->_wring_SEQIDX;
	}
	1;
    }
    else {
	# delete all matches
	$status = $self->DELETE($key);
	1;
    }
    return 0 if $status;
    return 1;
}

=head2 BDB API Emulation : sequential access

=head2 seq()

 Title   : seq
 Usage   : 
 Function: as in DB_File
 Returns : 
 Args    : 

=cut

sub seq {
    my $self = shift;
    my ($key, $value, $flags) = @_;
    return 1 unless $flags;
    # to modify $key, set $_[0]
    # to modify $value, set $_[1]
    my $SEQIDX = $self->SEQIDX;
    my $CURSOR = $self->CURSOR;
    $self->_wring_SEQIDX unless defined $$SEQIDX[$$CURSOR];
    for ($flags) {
	$_ eq R_CURSOR && do {
	    last;
	};
	$_ eq R_FIRST && do {
	    $$CURSOR  = 0;
	    $self->_wring_SEQIDX() unless defined $$SEQIDX[$$CURSOR];
	    last;
	};
	$_ eq R_LAST && do {
	    $$CURSOR = $#$SEQIDX;
	    # the ff is necessary: the original cursor may have been
	    # defined, but not the new one
	    $self->_wring_SEQIDX() unless defined $$SEQIDX[$$CURSOR];
	    last;
	};
	$_ eq R_NEXT && do {
	    return 1 if ($$CURSOR >= $#$SEQIDX);
	    ($$CURSOR)++;
	    # the ff is necessary: the original cursor may have been
	    # defined, but not the new one
	    $self-> _wring_SEQIDX() unless defined $$SEQIDX[$$CURSOR];
	    last;
	};
	$_ eq R_PREV && do {
	    return 1 if $$CURSOR == 0;
	    ($$CURSOR)--;
	    # the ff is necessary: the original cursor may have been
	    # defined, but not the new one
	    $self->_wring_SEQIDX() unless defined $$SEQIDX[$$CURSOR];
	    last;
	};
    }
    # get by pk, set key and value.
    $self->get_seq_sth->execute($$SEQIDX[$$CURSOR]);
    my $ret = $self->get_seq_sth->fetch;
    ($_[0], $_[1]) = (($self->ref eq 'ARRAY' ? $$CURSOR : $$ret[0]), $$ret[1]);
    return 0;
}

# remove undefs from @SEQIDX
# taking care of the cursor
sub _wring_SEQIDX {
    my $self = shift;
    my $SEQIDX = $self->SEQIDX;
    my $CURSOR = $self->CURSOR;
    my ($i, $j, @a);
    $j = 0;
    for $i (0..$#$SEQIDX) {
	if (defined $$SEQIDX[$i]) {
	    $$CURSOR = $j if $$CURSOR == $i;
	    $a[$j++] = $$SEQIDX[$i];
	}
	else {
	    $$CURSOR = $i+1 if $$CURSOR == $i;
	}
    }
    @$SEQIDX = @a;
    return;
}
    

=head2 _get_pk()

 Title   : _get_pk
 Usage   : 
 Function: provide an unused primary key integer for seq access
 Returns : scalar int
 Args    : none

=cut

sub _get_pk {
    my $self = shift;
    # do the primary key auditing for the cursor functions...
    return ++$AUTOPK;
}

=head2 _last_pk

 Title   : _last_pk
 Usage   : $obj->_last_pk($newval)
 Function: the primary key integer returned on the last FETCH
 Example : 
 Returns : value of _last_pk (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub _last_pk {
    my $self = shift;
    
    return $self->{'_last_pk'} = shift if @_;
    return $self->{'_last_pk'};
}

=head2 BDB API Emulation: dup

=head2 dup

 Title   : dup
 Usage   : $obj->dup($newval)
 Function: flag to indicate whether duplicate keys are handled
           (compare R_DUP of DB_File)
 Example : 
 Returns : value of dup (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub dup {
    my $self = shift;
    return $self->{'dup'} = shift if @_;
    return $self->{'dup'};
}

=head2 get_dup()

 Title   : get_dup
 Usage   : 
 Function: as in DB_File
 Returns : 
 Args    : 

=cut

sub get_dup {
    my $self = shift;
    my ($key, $aa) = @_;
    return unless $self->dbh;
    unless ($self->dup) {
	warn("DB not created in dup context; ignoring");
	return;
    }
    $self->get_sth->execute($key);
    my $ret = $self->get_sth->fetchall_arrayref;
    return scalar @$ret unless wantarray;
    my @ret = map {$_->[0]} @$ret;
    if (!$aa) {
	return @ret;
    }
    else {
	my %h;
	$h{$_}++ for @ret;
	return %h;
    }
}

=head2 find_dup()

 Title   : find_dup
 Usage   : 
 Function: as in DB_File
 Returns : 
 Args    : 
 Note    : no cursor functionality
=cut

sub find_dup {
    my $self = shift;
    my ($key, $value) = @_;
    return unless $self->dbh;
    unless ($self->dup) {
	warn("DB not created in dup context; ignoring");
	return;
    }
    $self->get_sth->execute($key);
    my $ret = $self->get_sth->fetchall_arrayref;
    return 0 if grep(/^$value$/, map {$_->[0]} @$ret);
    return 1;
}

=head2 del_dup()

 Title   : del_dup
 Usage   : 
 Function: as in DB_File
 Returns : 
 Args    : 

=cut

sub del_dup {
    my $self = shift;
    my ($key, $value) = @_;
    return unless $self->dbh;
    unless ($self->dup) {
	warn("DB not created in dup context; ignoring");
	return;
    }
    if ($self->dbh->do("DELETE FROM hash WHERE id = '$key' AND obj = '$value'")) {
	return 0; # success
    }
    else {
	return 1; # fail
    }
}

1;
