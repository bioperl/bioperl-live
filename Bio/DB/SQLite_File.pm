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

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;

use base qw(Bio::Root::Root );

 BEGIN {
     unless (eval "require DBD::SQLite; 1") {
 	Bio::Root::Root->throw( "SQLite_File requires DBD::SQLite" );
     }
 }

@Bio::DB::SQLite_File::ISA = qw(Bio::DB::SQLiteHASH);

sub TIEHASH {
    my $class = shift;
    my ($file, $flags, $mode, @args) = @_;
    my $keep = pop @args;
    my $self = Bio::DB::SQLiteHASH->TIEHASH($file, $flags, $keep);
    bless ($self, $class);
}

1;

package Bio::DB::SQLiteHASH;
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

sub TIEHASH {
    my $class = shift;
    my ($file, $flags, $keep) = @_;
    my $self = {};
    bless($self, $class);
    my $infix;
    # you'll love this...
    for ($flags) {
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

    my $dbh = DBI->connect("dbi:SQLite:dbname=".$self->file,"","",
			   {RaiseError => 1, AutoCommit => 1});
    $self->dbh( $dbh );
    $self->dbh->do("CREATE TABLE IF NOT EXISTS hash $hash_tbl");
    return $self;
}

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
    $self->put_sth->execute($key, $value);
}

sub DELETE {
    my $self = shift;
    my $key = shift;
    return unless $self->dbh;
    my $oldval = $self->FETCH($key);
    $self->dbh->do("DELETE FROM hash WHERE id = '$key'") if $oldval;
    return $oldval;
}

sub CLEAR {
    my $self = shift;
    return unless $self->dbh;
    $self->dbh->do("DELETE FROM hash");
    return 1;
}

sub EXISTS {
    my $self = shift;
    my $key = shift;
    return unless $self->dbh;
    return $self->FETCH($key) ? 1 : 0;
}

sub FIRSTKEY {
    my $self = shift;
    return unless $self->dbh;
    my $ids = $self->dbh->selectall_arrayref("SELECT id FROM hash");
    return unless $ids;
    return $self->_keys($ids);
}

sub NEXTKEY {
    my $self = shift;
    my $lastkey = shift;
    return $self->_keys;
}

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

=head2 Bio::DB::SQLiteHASH 

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
