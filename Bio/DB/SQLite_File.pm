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

# Let the code begin...

package Bio::DB::SQLite_File;
use strict;
use warnings;



BEGIN {
    unless (eval "require DBD::SQLite; 1") {
 	Bio::Root::Root->throw( "SQLite_File requires DBD::SQLite" );
     }
    use Exporter;
    our @ISA = qw( Exporter );
    our @EXPORT = qw( $DB_HASH $DB_BTREE $DB_RECNO R_DUP);
}

our $DB_HASH = { 'type' => 'HASH' };
our $DB_BTREE = { 'type' => 'BINARY' };
our $DB_RECNO = { 'type' => 'RECNO' };
sub R_DUP { 32678 }

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

use Bio::Root::Root;
use DBI;
use File::Temp qw( tempfile );
use Fcntl qw(O_CREAT O_RDWR O_RDONLY);
use base qw( Bio::Root::Root );
our %KEYTBL; # for accounting for pop/push/shift/unshift without
            # modifying the db
our $AUTOKEY = 0;

sub TIEHASH {
    my $class = shift;
    my ($file, $flags, $mode, $index, $keep) = @_;
    my $self = {};
    bless($self, $class);
    $self->{ref} = 'HASH';
    my $infix;
    my $fh;
    # db file handling
    if ($file) {
	# you'll love this...
	for ($flags) {
	    !defined && do {
		$infix = ">";
		last;
	    };
	    ($_ & O_CREAT) && do {
		$infix = (-e $file ? '<' : '>');
	    };
	    ($_ & O_RDWR) && do {
		$infix = '+'.($infix ? $infix : '>');
	    };
	    do { # O_RDONLY
		$infix = '<' unless $infix;
	    };
	}
	open($fh, $infix, $file) or $self->throw("Can't open db file: $!");
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
      obj     blob not null
    );
END
    my $hash_tbl_dup = <<END;
    ( id      blob,
      obj     blob,
      dup     integer primary key autoincrement
    );
END
    my $create_idx = <<END;
    CREATE INDEX IF NOT EXISTS id_idx ON hash ( id );
END
    my $create_idx_dup = <<END;
    CREATE INDEX IF NOT EXISTS id_idx ON hash ( id, dup );
END
    my $dbh = DBI->connect("dbi:SQLite:dbname=".$self->file,"","",
			   {RaiseError => 1, AutoCommit => 1});
    $self->dbh( $dbh );
    # rudimentary DB_File emulation:
    if (defined $index) {
	$self->throw("Index selector must be a hashref") unless ref($index) eq 'HASH';
	for ($index->{'type'}) {
	    $_ eq 'BINARY' && do {
		if ($index->{flags} eq R_DUP()) {
		    $self->dup(1);
		    $self->dbh->do("CREATE TABLE IF NOT EXISTS hash $hash_tbl_dup");
		    $self->dbh->do($create_idx_dup);
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
	    !defined && do {
		$self->dbh->do("CREATE TABLE IF NOT EXISTS hash $hash_tbl");
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
	for ($flags) {
	    !defined && do {
		$infix = ">";
		last;
	    };
	    ($_ & O_CREAT) && do {
		$infix = (-e $file ? '<' : '>');
	    };
	    ($_ & O_RDWR) && do {
		$infix = '+'.($infix ? $infix : '>');
	    };
	    do { # O_RDONLY
		$infix = '<' unless $infix;
	    };
	}
	open($fh, $infix, $file) or $self->throw("Can't open db file: $!");
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
    if (!$self->{ref} or $self->{ref} eq 'HASH') {
	$self->get_sth->execute($key);
    }
    elsif ($self->{ref} eq 'ARRAY') {
	$self->get_sth->execute($self->get_key($key));
    }
    else { # type not recognized
	return; 
    }
    my $ret = $self->get_sth->fetch;
    return $ret unless $ret;
    return $ret->[0];
}

sub STORE {
    my $self = shift;
    my ($key, $value) = @_;
    return unless $self->dbh;
    if ( !defined $self->{ref} or $self->{ref} eq 'HASH' ) {
	$value =~ s{'}{::}g;
	$self->put_sth->execute($key, $value);
    }
    elsif ( $self->{ref} eq 'ARRAY' ) {
	# need to check if the key is already present
	if ($self->is_empty($self->get_key($key))) {
	    $self->put_sth->execute($self->get_key($key), $value);
	}
	else {
	    # escape quotes (or SQL will barf)
	    $value =~ s{'}{::}g;
	    my $sth = $self->dbh->do("UPDATE hash SET obj = '$value' WHERE id = ".$self->get_key($key));

	}
	$self->inc;
    }
    $value;
}

sub DELETE {
    my $self = shift;
    my $key = shift;
    return unless $self->dbh;
    my $oldval = $self->FETCH($key);
    $self->dbh->do("DELETE FROM hash WHERE id = '$key'") if $oldval;
    return $oldval;
}

sub EXISTS {
    my $self = shift;
    my $key = shift;
    return unless $self->dbh;
    return $self->FETCH($key) ? 1 : 0;
}

# hash methods

sub CLEAR {
    my $self = shift;
    return unless $self->dbh;
    return if ($self->{ref} and $self->{ref} ne 'HASH');
    $self->dbh->do("DELETE FROM hash");
    return 1;
}

sub FIRSTKEY {
    my $self = shift;
    return unless $self->dbh;
    return if ($self->{ref} and $self->{ref} ne 'HASH');
    my $ids = $self->dbh->selectall_arrayref("SELECT id FROM hash");
    return unless $ids;
    return $self->_keys($ids);
}

sub NEXTKEY {
    my $self = shift;
    my $lastkey = shift;
    return unless $self->dbh;
    return if ($self->{ref} and $self->{ref} ne 'HASH');
    return $self->_keys;
}

# array methods

sub FETCHSIZE {
    my $self = shift;
    return unless $self->dbh;
    return if (!$self->{ref} or $self->{ref} ne 'ARRAY');
    $self->len;
}

sub STORESIZE {
    my $self = shift;
    my $count = shift;
    return unless $self->dbh;
    return if (!$self->{ref} or $self->{ref} ne 'ARRAY');
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
    return if (!$self->{ref} or $self->{ref} ne 'ARRAY');
    $self->get_sth->execute($self->get_key($self->len-1));
    my $ret = $self->get_sth->fetch;
    $self->dbh->do("DELETE FROM hash WHERE id = ".$self->get_key($self->len-1));
    # bookkeeping
    $self->rm_key($self->len-1);
    $self->dec;
    return defined $ret ? $ret->[0] : $ret;
}

sub PUSH {
    my $self = shift;
    my @values = @_;
    return unless $self->dbh;
    return if (!$self->{ref} or $self->{ref} ne 'ARRAY');
    my $ret = @values;
    my $beg = $self->len;
    my $end = $self->len + @values - 1;
    # if the bookkeeping has been right, 
    # none of these indices should be assigned
    # in %KEYTBL...
    for my $i ($beg..$end) {
	$self->put_sth->execute($self->get_key($i), shift @values);
    }
    $self->{len} += $ret;
    return $ret;
}

sub SHIFT {
    my $self = shift;
    return unless $self->dbh;
    return if (!$self->{ref} or $self->{ref} ne 'ARRAY');
    $self->get_sth->execute( $self->get_key(0) );
    my $ret = $self->get_sth->fetch;
    $self->dbh->do("DELETE FROM hash WHERE id = ".$self->get_key(0));
    # bookkeeping
    $self->shift_key;
    return defined $ret ? $ret->[0] : $ret;
}

sub UNSHIFT {
    my $self = shift;
    my @values = @_;
    return if (!$self->{ref} or $self->{ref} ne 'ARRAY');
    my $n = @values;
    return unless $self->dbh;
    for ($self->unshift_key($n)) {
	$self->put_sth->execute($_,shift @values);
    }
    $self->{len}+=$n;
    return $n;
}

sub SPLICE {
    my $self = shift;
    return if !$self->dbh;
    return if (!$self->{ref} or $self->{ref} ne 'ARRAY');
    die "Won't do splice yet."
}
# destructors

sub UNTIE {
    my $self = shift;
    my $count = shift;
    $self->throw( "untie attempted while $count inner references still exist" ) if ($count);}

sub DESTROY {
    my $self = shift;
    $self->put_sth->finish if $self->put_sth;
    $self->get_sth->finish if $self->get_sth;
    undef $self->{put_sth};
    undef $self->{get_sth};
    $self->dbh->disconnect;
    $self->_fh->close() if $self->_fh;
    unlink $self->file if (!$self->keep && $self->_fh);
    $self->throw("SQLite_File unlink issue: $!") if $!;
    1;
}

=head2 Bio::DB::SQLiteTIED

This module implements the tie; see L<perltie> for details. Other object methods follow; they are essentially internal, but can be accessed as follows:

 # %db is a tied hash...
 $SQLite_handle = (tied %db)->dbh;
 $db_file = (tied %db)->file;

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
    if (!$self->{'get_sth'}) {
	$self->{'get_sth'} =  $self->dbh->prepare("SELECT obj FROM hash WHERE id = ?");
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
	$self->{'put_sth'} = $self->dbh->prepare("INSERT INTO hash (id, obj) VALUES ( ?, ? )");
    }
    return $self->{'put_sth'};
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
    return $self->{_ref};
}

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

=head2 BDB API Emulation



=head2 get()

 Title   : get
 Usage   : $X->get($key,\$value)
 Function: see DB_File
 Returns : 
 Args    : 

=cut

sub get {
    my $self = shift;
    my ($key, $value) = @_;
    $$value = $self->FETCH($key);
    return 0 if defined $$value;
    return 1;
}



=head2 put()

 Title   : put
 Usage   : $X->put($key, $value)
 Function: see DB_File
 Returns : 
 Args    : 

=cut

sub put {
    my $self = shift;
    my ($key, $value) = @_;
    my $status = $self->STORE($key, $value);
    return 0 if $status;
    return 1;
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
    shift->{len};
}

sub get_key {
    my $self = shift;
    my $index = shift;
    return $KEYTBL{$index} if defined $KEYTBL{$index};
    $KEYTBL{$index} = $AUTOKEY++;
}

sub shift_key {
    my $self = shift;
    return if !$self->len;
    for my $i (0..$self->len-1) {
	$KEYTBL{$i} = $KEYTBL{$i+1}; # should undef the last elt
    }
    $self->dec;
    return;
}

# returns the set of new db ids to use
sub unshift_key {
    my $self = shift;
    my $n = shift; # number of new elements
    my (@new, @old);
    for my $i (0..$n-1) {
	push @new, $AUTOKEY++;
    }
    @old=@KEYTBL{ (0..$self->len-1) };
    @KEYTBL{ (0..$n-1) } = @new;
    @KEYTBL{ ($n..$self->len+$n-1) } = @old;
    return @new;
}

sub rm_key {
    my $self = shift;
    my $index = shift;
    unless (delete $KEYTBL{$index}) {
	warn "Element $index did not exist";
    }
}



=head2 is_empty()

 Title   : is_empty
 Usage   : 
 Function: check if the DB array slot is already occupied
           by examining the values of %KEYTBL
 Returns : 
 Args    : actual hash.id value

=cut

sub is_empty {
    my $self = shift;
    my $k = shift;
    # TODO: better to sort only when necessary
    my @ids = sort {$a <=> $b} values %KEYTBL;
    my $ret = 1;
    foreach (@ids) {
	if ( $k == $_ ) {
	    $ret = 0;
	    last;
	}
	elsif ( $k < $_ ) {
	    last;
	}
    }
    return $ret;
}


1;
