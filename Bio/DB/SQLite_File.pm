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
 
 # use as an option in AnyDBM_File
 @AnyDBM_File::ISA = qw( DB_File Bio::DB::SQLite_File SDBM );
 my %db;
 tie(%db, 'AnyDBM_File', 'my.db', $flags, $mode, @dbmargs)
 
=head1 DESCRIPTION

This is a simple wrapper around L<Bio::DB::SQLiteHASH>, also in this
file, that allows a hash to be tied to a SQLite DB via L<DBI> plus
L<DBD::SQLite>. The class does not currently have the functionality of
a L<DB_File>, but is suitable for L<Bio::DB::FileCache>, for which it
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
    use vars qw( $DB_HASH $DB_BTREE $DB_RECNO );
    use Exporter;
    our @ISA = qw( Exporter );
    our @EXPORT = qw( $DB_HASH $DB_BTREE $DB_RECNO );
}

# Object preamble - inherits from Bio::Root::Root
use Bio::Root::Root;

use base qw(Bio::Root::Root );

push @Bio::DB::SQLite_File::ISA, qw(Bio::DB::SQLiteTIED);
$DB_HASH = { 'type' => 'HASH' };
$DB_BTREE = { 'type' => 'BINARY' };
$DB_RECNO = { 'type' => 'RECNO' };

sub TIEHASH {
    my $class = shift;
    my ($file, $flags, $mode, $index, @args) = @_;
    my $keep = pop @args;
    my $self = TIEHASH($file, $flags, $index, $keep);
    $self->{_ref} = 'HASH';
    bless ($self, $class);
}

sub TIEARRAY {
    my $class = shift;
    my ($file, $flags, $mode, $index, @args) = @_;
    my $keep = pop @args;
    my $self = TIEARRAY($file, $flags, $index, $keep);
    $self->{_ref} = 'ARRAY';
    bless ($self, $class);
}

# HASH or ARRAY?
sub ref {
    my $self = shift;
    return $self->{_ref};
}
1;

package Bio::DB::SQLiteTIED;
use strict;
use warnings;

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

use DBI;
use File::Temp qw( tempfile );
use Fcntl;

my %KEYTBL; # for accounting for pop/push/shift/unshift without
            # modifying the db
my $AUTOKEY = 0;

sub TIEHASH {
    my $class = shift;
    my ($file, $flags, $index, $keep) = @_;
    my $self = {};
    bless($self, $class);
    my $infix;
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
    my $fh;
    # db file handling
    if ($file) {
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
    
     my $hash_tbl = <<END;
    (
      id      blob,
      obj     blob not null
    );
END

    my $create_idx = <<END;
    CREATE INDEX IF NOT EXISTS id_idx ON hash ( id );
END
    
    my $dbh = DBI->connect("dbi:SQLite:dbname=".$self->file,"","",
			   {RaiseError => 1, AutoCommit => 1});
    $self->dbh( $dbh );
    $self->dbh->do("CREATE TABLE IF NOT EXISTS hash $hash_tbl");
    # rudimentary DB_File emulation:
    if (defined $index) {
	$self->throw("Index selector must be a hashref") unless ref($index) eq 'HASH';
	for ($index->{'type'}) {
	    $_ eq 'BINARY' && do {
		$self->dbh->do($create_idx);
	    };
	    $_ eq 'HASH' && do {
		last;
	    };
	    $_ eq 'RECNO' && do {
		last;
	    };
	    !defined && do {
		last;
	    };
	}
    }
    return $self;
}

sub TIEARRAY {
    my $class = shift;
    my ($file, $flags, $index, $keep) = @_;
    my $self = {};
    bless($self, $class);
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
    my $fh;
    # db file handling
    if ($file) {
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
    $self->dbh->do("CREATE TABLE IF NOT EXISTS hash $hash_tbl");
    # rudimentary DB_File emulation:
    if (defined $index) {
	$self->throw("Index selector must be a hashref") unless ref($index) eq 'HASH';
	for ($index->{'type'}) {
	    $_ eq 'BINARY' && do {
		$self->dbh->do($create_idx);
	    };
	    $_ eq 'HASH' && do {
		last;
	    };
	    $_ eq 'RECNO' && do {
		$self->dbh->do($create_idx);
	    };
	    !defined && do {
		last;
	    };
	}
    }
    $self->{i} = 0;
    return $self;
}



# common methods

sub FETCH {
    my $self = shift;
    my $key = shift;
    return unless $self->dbh;
    if (!$self->get_sth) {
	$self->get_sth( $self->dbh->prepare("SELECT obj FROM hash WHERE id = ?") );
    }
    $self->get_sth->execute($key);
    my $ret = $self->get_sth->fetch;
    return $ret unless $ret;
    return $ret->[0];
}

sub STORE {
    my $self = shift;
    my ($key, $value) = @_;
    return unless $self->dbh;
    if (!$self->put_sth) {
	$self->put_sth(
	    $self->dbh->prepare("INSERT INTO hash (id, obj) VALUES ( ?, ? )")
	    );
    }
    if ( !defined $self->{ref} or $self->{ref} eq 'HASH' ) {
	$self->put_sth->execute($key, $value);
    }
    elsif ( $self->{ref} eq 'ARRAY' ) {
	$self->put_sth->execute($self->get_key($key), $value);
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
    $self->len - 1;
}

sub STORESIZE {
    my $self = shift;
    return unless $self->dbh;
    return if (!$self->{ref} or $self->{ref} ne 'ARRAY');
    $self->len;
}

sub POP {
    my $self = shift;
    return unless $self->dbh;
    return if (!$self->{ref} or $self->{ref} ne 'ARRAY');
    if (!$self->get_sth) {
	$self->get_sth( $self->dbh->prepare("SELECT obj FROM hash WHERE id = ?") );
    }
    $self->get_sth->execute($self->get_key($self->len-1));
    my $ret = $self->get_sth->fetch;
    $self->dbh->do("DELETE FROM hash WHERE id = ".$self->get_key($self->len-1));
    # bookkeeping
    $self->rm_key($self->len-1);
    $self->dec;
    return $ret;
}

sub PUSH {
    my $self = shift;
    my @values = @_;
    return unless $self->dbh;
    return if (!$self->{ref} or $self->{ref} ne 'ARRAY');
    my $ret = @values;
    if (!$self->put_sth) {
	$self->put_sth(
	    $self->dbh->prepare("INSERT INTO hash (id, obj) VALUES ( ?, ? )")
	    );
    }
    my $beg = $self->len;
    my $end = $self->len + @values - 1;
    # if the bookkeeping has been right, 
    # none of these indices should be assigned
    # in %KEYTBL...
    for my $i ($beg..$end) {
	$self->put_sth->execute($self->get_key($i), shift @values)
    }
    return $ret;
}

sub SHIFT {
    my $self = shift;
    return unless $self->dbh;
    return if (!$self->{ref} or $self->{ref} ne 'ARRAY');
    if (!$self->get_sth) {
	$self->get_sth( $self->dbh->prepare("SELECT obj FROM hash WHERE id = ?") );
    }
    my $ret = $self->get_sth->execute( $self->get_key(0) );
    $self->dbh->do("DELETE FROM hash WHERE id = ".$self->get_key(0));
    # bookkeeping
    $self->shift_key;
    return $ret;
}

sub UNSHIFT {
    my $self = shift;
    my @values = @_;
    return if (!$self->{ref} or $self->{ref} ne 'ARRAY');
    my $n = @values;
    return unless $self->dbh;
    if (!$self->put_sth) {
	$self->put_sth(
	    $self->dbh->prepare("INSERT INTO hash (id, obj) VALUES ( ?, ? )")
	    );
    }
    for ($self->unshift_key($n)) {
	$self->put_sth->execute($_,shift @values);
    }
    return $n;
}

sub SPLICE {
    my $self = shift;
    return if !$self->dbh;
    return if (!$self->{ref} or $self->{ref} ne 'ARRAY');
    die "Won't do splice yet."
}

sub inc {
    my $self = shift;
    return ++($self->{len});
}

sub dec {
    my $self = shift;
    return $self->{len} ? --($self->{len}) : 0;
}

sub len {
    my $self = shift;
    $self->{len};
}

sub get_key {
    my $self = shift;
    my $index = shift;
    return defined $KEYTBL{$index} ? $KEYTBL{$index} = $AUTOKEY++ : $KEYTBL{$index};
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
    $self->{len}+=$n;
    return @new;
}

sub rm_key {
    my $self = shift;
    my $index = shift;
    unless (delete $KEYTBL{$index}) {
	warn "Element $index did not exist";
    }
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
 Function: select statement handle
 Example : 
 Returns : value of get_sth (a statement handle)
 Args    : [optional] new handle

=cut

sub get_sth {
    my $self = shift;
    
    return $self->{'get_sth'} = shift if @_;
    return $self->{'get_sth'};
}


=head2 put_sth

 Title   : put_sth
 Usage   : $obj->put_sth($new_handle)
 Function: insert statement handle
 Example : 
 Returns : value of put_sth (a statement handle)
 Args    : [optional] new handle

=cut

sub put_sth {
    my $self = shift;
    
    return $self->{'put_sth'} = shift if @_;
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

1;
