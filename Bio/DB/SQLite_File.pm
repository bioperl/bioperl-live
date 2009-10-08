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
use vars qw( $AUTOLOAD ) ;

BEGIN {
    unless (eval "require DBD::SQLite; 1") {
 	Bio::Root::Root->throw( "SQLite_File requires DBD::SQLite" );
     }
    use Fcntl qw(O_CREAT O_RDWR O_RDONLY);
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
sub R_DUP () { 32678 }
sub R_CURSOR () { 27 }
sub R_FIRST () { 7 }
sub R_LAST () { 15 }
sub R_NEXT () { 16 }
sub R_PREV () { 23 }
sub R_IAFTER () { 1 }
sub R_IBEFORE () { 3 }
sub R_NOOVERWRITE () { 20 }
sub R_SETCURSOR () { -100 }

#sub O_CREAT  () { 512 };
#sub O_RDWR () { 2 };
#sub O_RDONLY () { 0 };
sub O_SVWST (){ O_CREAT() | O_RDWR() };

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use DBI qw(:sql_types);
use File::Temp qw( tempfile );
use base qw( Bio::Root::Root );

$Bio::DB::SQLite_File::MAXPEND = 250;

our $AUTOKEY = 0;
# for providing DB_File seq functionality
our $AUTOPK = 0;

# statement tables
our %STMT = (
    HASH => {
	put     => "INSERT INTO hash (id, obj, pk) VALUES ( ?, ?, ? )",
	put_seq => "INSERT INTO hash (id, obj, pk) VALUES ( ?, ?, ? )",
	get     => "SELECT obj, pk FROM hash WHERE id = ?",
	get_seq => "SELECT id, obj FROM hash WHERE pk = ?",
	upd     => "UPDATE hash SET obj = ? WHERE id = ? AND pk = ?",
	upd_seq => "UPDATE hash SET id = ?, obj = ? WHERE pk = ?",
	del     => "DELETE FROM hash WHERE id = ?",
	del_seq => "DELETE FROM hash WHERE pk = ?",
	del_dup => "DELETE FROM hash WHERE id = ? AND obj = ?",
	sel_dup => "SELECT pk FROM hash WHERE id = ? AND obj = ?",
	part_seq=> "SELECT id, obj, pk FROM hash WHERE id >= ? LIMIT 1"
    },
    ARRAY => {
	put     => "INSERT INTO hash (id, obj) VALUES ( ?, ?)",
	put_seq => "INSERT INTO hash (obj, id) VALUES ( ?, ?)",
	get     => "SELECT obj, id FROM hash WHERE id = ?",
	get_seq => "SELECT id, obj FROM hash WHERE id = ?",
	upd     => "UPDATE hash SET obj = ? WHERE id = ?",
	upd_seq => "UPDATE hash SET obj = ? WHERE id = ?",
	del     => "DELETE FROM hash WHERE id = ?",
	del_seq => "DELETE FROM hash WHERE id = ?"
    }
    );

	

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
    $index ||= $DB_HASH;
    $self->{index} = $index;
    $self->{pending} = 0;
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
      id      blob collate NOCASE,
      obj     blob not null,
      pk      integer primary key autoincrement
    );
END
    my $create_idx = <<END;
    CREATE INDEX IF NOT EXISTS id_idx ON hash ( id, pk );
END
    my $dbh = DBI->connect("dbi:SQLite:dbname=".$self->file,"","",
			   {RaiseError => 1, AutoCommit => 0});
    $self->dbh( $dbh );
#    $dbh->do("PRAGMA journal_mode = OFF");
    if (defined $index) {
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
    $self->commit(1);
    return $self;
}

sub TIEARRAY {
    my $class = shift;
    my ($file, $flags, $mode, $index, $keep) = @_;
    my $self = {};
    bless($self, $class);
    $self->{ref} = 'ARRAY';
    $index ||= $DB_RECNO;
    $self->throw("Arrays must be tied to type RECNO") unless 
	$index->{type} eq 'RECNO';
    $self->{index} = $index;
    $self->{pending} = 0;
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
			   {RaiseError => 1, AutoCommit => 0});
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
    $self->commit;
    $self->{len} = 0;
    return $self;
}



# common methods

sub FETCH {
    my $self = shift;
    my $key = shift;
    return unless $self->dbh;
    $self->commit;
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
    my ($pk, $sth);
    if ( !defined $self->{ref} or $self->ref eq 'HASH' ) {
	if ( $self->dup ) { # allowing duplicates
	    $pk = $self->_get_pk;
	    $sth = $self->put_sth;
	    $sth->bind_param(1,$key);
	    $sth->bind_param(2,$value, SQL_BLOB);
	    $sth->bind_param(3,$pk);
	    $self->put_sth->execute();
	    push @{$self->SEQIDX}, $pk;
	}
	else { # no duplicates...
	    #need to check if key is already present
	    if ( $self->EXISTS($key) ) {
		$sth = $self->upd_sth;
		$sth->bind_param(1,$value, SQL_BLOB);
		$sth->bind_param(2,$key);
		$sth->bind_param(3,$self->_last_pk);
		$sth->execute();
	    }
	    else {
		$pk = $self->_get_pk;
		$sth = $self->put_sth;
		$sth->bind_param(1,$key);
		$sth->bind_param(2,$value, SQL_BLOB);
		$sth->bind_param(3,$pk);
		$sth->execute();
		push @{$self->SEQIDX}, $pk;
	    }
	}
	$self->{_stale} = 1;
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
    ++$self->{pending};
    $value;
}

sub DELETE {
    my $self = shift;
    my $key = shift;
    return unless $self->dbh;
    $self->_reindex if ($self->index->{type} eq 'BINARY' and $self->_index_is_stale);
    my $oldval;
    if (!$self->ref or $self->ref eq 'HASH') {
	return unless $self->get_sth->execute($key);
	my $ret = $self->get_sth->fetch;
	$oldval = $ret->[0];
#	$self->dbh->do("DELETE FROM hash WHERE id = '$key'");
	$self->del_sth->execute($key); # del on id
	# update the sequential side
	if ($ret->[1]) {
	    delete ${$self->SEQIDX}[_find_idx($ret->[1],$self->SEQIDX)];
	}
    }
    elsif ($self->ref eq 'ARRAY') {
	my $SEQIDX = $self->SEQIDX;
	if ($$SEQIDX[$key]) {
	    $oldval = $self->FETCH($$SEQIDX[$key]);
#	    $self->dbh->do("DELETE FROM hash WHERE id = '$$SEQIDX[$key]'");
	    $self->del_sth->execute($$SEQIDX[$key]); # del on id
	    $self->rm_idx($key);
	}
    }
    else {
	$self->throw( "tied type not recognized" );
    }
    ++$self->{pending};
    return $oldval;
}

sub EXISTS {
    my $self = shift;
    my $key = shift;
    return unless $self->dbh;
    $self->commit;
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
    $self->dbh->commit;
    my $sth = $self->dbh->prepare("DELETE FROM hash");
    $sth->execute;
    $self->dbh->commit;
    @{$self->SEQIDX} = ();
    return 1;
}

# hash methods

sub FIRSTKEY {
    my $self = shift;
    return unless $self->dbh;
    $self->commit;
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
    $self->commit;
    return if (!$self->{ref} or $self->ref ne 'ARRAY');
    $self->get_sth->execute($self->get_idx($self->len-1));
    my $ret = $self->get_sth->fetch;
#    $self->dbh->do("DELETE FROM hash WHERE id = ".$self->get_idx($self->len-1));
    $self->del_sth->execute($self->get_idx($self->len-1));
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
    ++$self->{pending};
    return $ret;
}

sub SHIFT {
    my $self = shift;
    return unless $self->dbh;
    $self->commit;
    return if (!$self->{ref} or $self->ref ne 'ARRAY');
    $self->get_sth->execute( $self->get_idx(0) );
    my $ret = $self->get_sth->fetch;
#    $self->dbh->do("DELETE FROM hash WHERE id = ".$self->get_idx(0));
    $self->del_sth->execute($self->get_idx(0));
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
    ++$self->{pending};
    return $n;
}

sub SPLICE {
    my $self = shift;
    return if !$self->dbh;
    return if (!$self->{ref} or $self->ref ne 'ARRAY');
    die "Won't do splice yet.";
    ++$self->{pending};
}

# destructors

sub UNTIE {
    my $self = shift;
    my $count = shift;
    $self->throw( "untie attempted while $count inner references still exist" ) if ($count);}

sub DESTROY {
    my $self = shift;
    $self->dbh->commit; #'hard' commit
    my $tbl = $STMT{$self->ref};
    # finish and destroy stmt handles
    for ( keys %$tbl ) {
	$self->{$_."_sth"}->finish if $self->{$_."_sth"};
	undef $self->{$_."_sth"};
    }
    # disconnect
    $self->throw($self->dbh->errstr) unless $self->dbh->disconnect;
    $self->{dbh}->DESTROY;
    undef $self->{dbh};
    # remove file if nec
    $self->_fh->close() if $self->_fh;
    unlink $self->file if (!$self->keep && $self->_fh);
    $self->throw("SQLite_File unlink issue: $!") if $!;
    undef $self;
    1;
}


=head2 SQL interface : internal

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



=head2 sth()

 Title   : sth
 Usage   : $obj->sth($stmt_descriptor)
 Function: statement handle generator
 Returns : a prepared DBI statement handle
 Args    : scalar string (statement descriptor)
 Note    : calls such as $obj->put_sth are autoloaded through
           this method
=cut

sub sth {
    my $self = shift;
    my $desc = shift;
    $self->throw("No active database handle") unless $self->dbh;
    my $tbl = $STMT{$self->ref};
    unless ($tbl) {
	$self->throw("Tied type '".$self->ref."' not recognized");
    }
    if (!$self->{"${desc}_sth"}) {
	$self->throw("Statement descriptor '$desc' not recognized for type ".$self->ref) unless grep(/^$desc$/,keys %$tbl);
	$self->{"${desc}_sth"} = $self->dbh->prepare($tbl->{$desc});
    }
    return $self->{"${desc}_sth"};
}

# autoload statement handle getters

sub AUTOLOAD {
    my $self = shift;
    my @pth = split(/::/, $AUTOLOAD); 
    my $desc = $pth[-1];
    unless ($desc =~ /^(.*?)_sth$/) {
	$self->throw("Subroutine '$AUTOLOAD' is undefined in ".__PACKAGE__);
    }
    $desc = $1;
    unless (grep /^$desc$/, keys %{$STMT{$self->ref}}) {
	$self->throw("Statement accessor ${desc}_sth not defined for type ".$self->ref);
    }
    $self->sth($desc);
}
=head2 commit()

 Title   : commit
 Usage   : 
 Function: commit transactions
 Returns : 
 Args    : commit(1) forces, commit() commits when
           number of pending transactions > $MAXPEND

=cut

sub commit {

    my $self = shift;
    if (@_ or ($self->{pending} > $Bio::DB::SQLite_File::MAXPEND)) {
	$self->warn("commit failed") unless $self->dbh->commit();
	$self->{pending} = 0;
    }
    return 1;
}



=head2 pending()

 Title   : pending
 Usage   : $obj->pending
 Function: count of pending (uncommitted) transactions
 Returns : scalar int
 Args    : none (rdonly)

=cut

sub pending {
    shift->{pending};
}

=head2 trace()

 Title   : trace
 Usage   : 
 Function: invoke the DBI trace logging service
 Returns : 
 Args    : scalar int trace level

=cut

sub trace {
    my $self = shift;
    my $level = shift;
    return unless $self->dbh;
    $level ||= 3;
    $self->dbh->{TraceLevel} = $level;
    $self->dbh->trace;
    return $level;
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

=head2 index

 Title   : index
 Usage   : $obj->index($newval)
 Function: access the index type structure ($DB_BTREE, $DB_HASH, $DB_RECNO)
           that initialized this instance
 Example : 
 Returns : value of index (a hashref)
 Args    : 

=cut

sub index {
    my $self = shift;
    return $self->{'index'};
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

=head2 Array object helper functions : internal

=cut

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
	$self->warn("Element $index did not exist");
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
    my ($sth, $do_cursor);
    for ($flags) {
	(!defined || $_ == R_SETCURSOR) && do { # put or upd
	    if ($self->dup) { # just make a new one
		$pk = $self->_get_pk;
		$sth = $self->put_seq_sth;
		$do_cursor = sub {
		    push @$SEQIDX, $pk;
		    $$CURSOR = $#$SEQIDX if $flags;
		    $self->_reindex if $self->index->{type} eq 'BINARY';
		};
	    }
	    else {
		$self->FETCH($key);
		$pk = $self->_last_pk || $self->_get_pk;
		$sth = ($self->_last_pk ? 
			   $self->upd_seq_sth :
			   $self->put_seq_sth);
		$do_cursor = sub {
		    push @$SEQIDX, $pk if !$self->_last_pk;
		    $flags && do { # R_SETCURSOR
			if ( $pk = $$SEQIDX[-1] ) {
			    $$CURSOR = $#$SEQIDX;
			}
			else {
			    $$CURSOR = _find_idx($pk, $SEQIDX);
			};
			$self->_reindex if $self->index->{type} eq 'BINARY';
		    };
		};
	    }
	    last;
	};
	$_ == R_IAFTER && do {
	    $self->_wring_SEQIDX unless $$SEQIDX[$$CURSOR];
	    # duplicate protect
	    return 1 unless ($self->ref eq 'ARRAY') || $self->dup || !$self->EXISTS($key);
	    $self->throw("R_IAFTER flag meaningful only for RECNO type") unless
		$self->index->{type} eq 'RECNO';
	    $pk = $self->_get_pk;
	    $sth = $self->put_seq_sth;
	    $_[0] = $$CURSOR+1;
	    $do_cursor = sub {
		if ($$CURSOR == $#$SEQIDX) {
		    push @$SEQIDX, $pk;
		}
		else {
		    splice(@$SEQIDX,$$CURSOR,0,$pk);
		}
	    };
	    last;
	};
	$_ == R_IBEFORE && do {
	    $self->_wring_SEQIDX unless $$SEQIDX[$$CURSOR];
	    # duplicate protect
	    return 1 unless ($self->ref eq 'ARRAY') || $self->dup || !$self->EXISTS($key);
	    $self->throw("R_IBEFORE flag meaningful only for RECNO type") unless
		$self->index->{type} eq 'RECNO';
	    $pk = $self->_get_pk;
	    $sth = $self->put_seq_sth;
	    $_[0] = $$CURSOR;
	    $do_cursor = sub {
		if ($$CURSOR) {
		    splice(@$SEQIDX,$$CURSOR-1,0,$pk);
		}
		else {
		    unshift(@$SEQIDX, $pk);
		}
		$$CURSOR++; # preserve cursor
	    };
	    last;
	};
	$_ == R_CURSOR && do { # upd only
	    $self->_wring_SEQIDX unless $$SEQIDX[$$CURSOR];
	    # duplicate protect
	    return 1 unless ($self->ref eq 'ARRAY') || $self->dup || !$self->EXISTS($key);
	    $pk = $$SEQIDX[$$CURSOR];
	    $sth = $self->upd_seq_sth;
	    $do_cursor = sub {
		$self->_reindex if $self->index->{type} eq 'BINARY';
	    };
	    last;
	};
	$_ == R_NOOVERWRITE && do { # put only/add to the "end"
	    #will create a duplicate if $self->dup is set!
	    return 1 unless ($self->ref eq 'ARRAY') || $self->dup || !$self->EXISTS($key);
	    $pk = $self->_get_pk;
	    $sth = $self->put_seq_sth;
	    $do_cursor = sub {
		push @$SEQIDX, $pk;
		$self->_reindex if $self->index->{type} eq 'BINARY';
	    };
	    last;
	};
    }
    if ($self->ref eq 'ARRAY') {
	$sth->bind_param(1, $value, SQL_BLOB);
	$sth->bind_param(2, $pk);
    }
    else {
	$sth->bind_param(1, $key);
	$sth->bind_param(2, $value, SQL_BLOB);
	$sth->bind_param(3, $pk);
    }
    $status = !$sth->execute;
    $do_cursor->() if !$status;
    $self->{pending} = 1;
    $self->{_stale} = 0 if $self->index->{type} eq 'BINARY';
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
    $self->_reindex if ($self->index->{type} eq 'BINARY' and $self->_index_is_stale);
    my $SEQIDX = $self->SEQIDX;
    my $CURSOR = $self->CURSOR;
    my $status;
    if ($flags eq R_CURSOR) {
	_wring_SEQIDX($self->SEQIDX) unless $$SEQIDX[$$CURSOR];
	my $pk = $$SEQIDX[$$CURSOR];
#	my $col = ($self->ref eq 'ARRAY' ? 'id' : 'pk');
#	$status = $self->dbh->do("DELETE FROM hash WHERE $col = $pk");
	$self->del_seq_sth->execute($pk);
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
    $self->{_stale} = 1;
    $self->{pending} = 1;
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
    $self->commit;
    my $status;
    # to modify $key, set $_[0]
    # to modify $value, set $_[1]
    $self->_reindex if ($self->index->{type} eq 'BINARY' and $self->_index_is_stale);
    my $SEQIDX = $self->SEQIDX;
    my $CURSOR = $self->CURSOR;
    for ($flags) {
	$_ eq R_CURSOR && do {
	    last;
	};
	$_ eq R_FIRST && do {
	    $$CURSOR  = 0;
	    last;
	};
	$_ eq R_LAST && do {
	    $$CURSOR = $#$SEQIDX;
	    last;
	};
	$_ eq R_NEXT && do {
	    return 1 if ($$CURSOR >= $#$SEQIDX);
	    ($$CURSOR)++;
	    last;
	};
	$_ eq R_PREV && do {
	    return 1 if $$CURSOR == 0;
	    ($$CURSOR)--;
	    last;
	};
    }
    $self->_wring_SEQIDX() unless defined $$SEQIDX[$$CURSOR];
    # get by pk, set key and value.
    if (($flags == R_CURSOR ) && $self->ref eq 'HASH') {
	$status = $self->partial_match($key, $value);
	$_[0] = $key; $_[1] = $value;
	return $status;
    }
    else {
	$self->get_seq_sth->execute($$SEQIDX[$$CURSOR]);
	my $ret = $self->get_seq_sth->fetch;
	($_[0], $_[1]) = (($self->ref eq 'ARRAY' ? $$CURSOR : $$ret[0]), $$ret[1]);
    }
    return 0;
}

=head2 sync()

 Title   : sync
 Usage   : 
 Function: stub for BDB sync 
 Returns : 0
 Args    : 

=cut

sub sync { !shift->commit };

=head2 BDB API Emulation : internals

=head2 partial_match()

 Title   : partial_match
 Usage   : 
 Function: emulate the partial matching of DB_File::seq() with
           R_CURSOR flag
 Returns : 
 Args    : $key

=cut

sub partial_match {
    my $self = shift;
    my ($key, $value) = @_;
    my ($status,$ret, $pk);
    unless ($self->ref ne 'ARRAY') {
	$self->throw("Partial matches not meaningful for arrays");
    }
    my $SEQIDX = $self->SEQIDX;
    my $CURSOR = $self->CURSOR;
    $status = !$self->part_seq_sth->execute( $key );
    if (!$status) { # success
	if ($ret = $self->{part_seq_sth}->fetch) {
	    $_[0] = $ret->[0]; $_[1] = $ret->[1];
	    $pk = $ret->[2];
 	    unless (defined($$CURSOR = _find_idx($pk,$SEQIDX))) {
		$self->throw("Primary key value disappeared!");
	    }
	    return 0;
	}
    }
    return 1;
}

=head2 _index_is_stale()

 Title   : _index_is_stale
 Usage   : 
 Function: predicate indicating whether a _reindex has been
           performed since adding or updating the db
 Returns : 
 Args    : none

=cut

sub _index_is_stale {
    my $self = shift;
    return $self->{_stale};
}

=head2 _reindex()

 Title   : _reindex
 Usage   : 
 Function: reorder SEQIDX to reflect BTREE ordering,
           preserving cursor
 Returns : true on success
 Args    : none

=cut

sub _reindex {
    my $self = shift;
    my ($q, @order);
    my $SEQIDX = $self->SEQIDX;
    my $CURSOR = $self->CURSOR;
    $self->_wring_SEQIDX;
    $q = $self->dbh->selectall_arrayref("SELECT pk, id FROM hash ORDER BY id");
    unless ($q) {
	return 0;
    }
    @order = map { $$_[0] } @$q;
    $$CURSOR = _find_idx($$SEQIDX[$$CURSOR],\@order);
    $self->{SEQIDX} = \@order;
    $self->{_stale} = 0;
    return 1;
}

=head2 _find_idx()

 Title   : _find_idx
 Usage   : 
 Function: search of array for index corresponding
           to a given value
 Returns : scalar int (target array index)
 Args    : scalar int (target value), array ref (index array)

=cut


sub _find_idx {
    my ($pk, $seqidx) = @_;
    my $i;
    for (0..$#$seqidx) {
	$i = $_;
	next unless defined $$seqidx[$_];
	last if $pk == $$seqidx[$_];
    }
    return (defined $$seqidx[$i] and $pk == $$seqidx[$i] ? $i : undef);
}

=head2 _wring_SEQIDX()

 Title   : _wring_SEQIDX
 Usage   : 
 Function: remove undef'ed values from SEQIDX,
           preserving cursor
 Returns : 
 Args    : none

=cut

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
    $self->commit;
    unless ($self->dup) {
	$self->warn("DB not created in dup context; ignoring");
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
    $self->commit;
    unless ($self->dup) {
	$self->warn("DB not created in dup context; ignoring");
	return;
    }
    $self->get_sth->execute($key);
    my $ret = $self->get_sth->fetchall_arrayref;
    return 0 if grep(/^$value$/, map {$_->[0]} @$ret);
    return;
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
    my $ret;
    return unless $self->dbh;
    unless ($self->dup) {
	$self->warn("DB not created in dup context; ignoring");
	return;
    }
#    my $sth = $self->dbh->prepare("SELECT pk FROM hash WHERE id = '$key' AND obj = ?");
    $self->sel_dup_sth->bind_param(1, $key);
    $self->sel_dup_sth->bind_param(2, $value, SQL_BLOB);
    $self->sel_dup_sth->execute;
    $ret = $self->sel_dup_sth->fetchall_arrayref;
    unless ($ret) {
	return 1;
    }
#    $sth = $self->dbh->prepare("DELETE FROM hash WHERE id = '$key' and obj = ?");
    
    $self->del_dup_sth->bind_param(1, $key);
    $self->del_dup_sth->bind_param(2, $value, SQL_BLOB);
    if ($self->del_dup_sth->execute) {
	# update SEQIDX
	foreach (map { $$_[0] } @$ret) {
	    delete ${$self->SEQIDX}[_find_idx($_,$self->SEQIDX)];
	}
	$self->_wring_SEQIDX;
	$self->{pending} = 1;
	return 0; # success
    }
    else {
	return 1; # fail
    }
}

1;

package Bio::DB::SQLite_File::HASHINFO;
use strict;
use warnings;

# a hashinfo class stub
sub new {
    my $class = shift;
    my $self = bless({}, $class);
    $self->{type} = 'HASH';
    return $self;
}

1;

package Bio::DB::SQLite_File::BTREEINFO;
use strict;
use warnings;

# a btreeinfo class stub
sub new {
    my $class = shift;
    my $self = bless({}, $class);
    $self->{type} = 'BINARY';
    return $self;
}

1;

package Bio::DB::SQLite_File::RECNOINFO;
use strict;
use warnings;

# a recnoinfo class stub
sub new {
    my $class = shift;
    my $self = bless({}, $class);
    $self->{type} = 'RECNO';
    return $self;
}

1;

