package Bio::DB::GFF::Adaptor::dbi::caching_handle;

use strict;
use DBI;
use vars '$AUTOLOAD';
use base qw(Bio::Root::Root);

=head1 NAME

Bio::DB::GFF::Adaptor::dbi::caching_handle -- Cache for database handles

=head1 SYNOPSIS

 use Bio::DB::GFF::Adaptor::dbi::caching_handle;
 $db  = Bio::DB::GFF::Adaptor::dbi::caching_handle->new('dbi:mysql:test');
 $sth = $db->prepare('select * from foo');
 @h   = $sth->fetch_rowarray;
 $sth->finish

=head1 DESCRIPTION

This module handles a pool of database handles.  It was motivated by
the MYSQL driver's {mysql_use_result} attribute, which dramatically
improves query speed and memory usage, but forbids additional query
statements from being evaluated while an existing one is in use.

This module is a plug-in replacement for vanilla DBI.  It
automatically activates the {mysql_use_result} attribute for the mysql
driver, but avoids problems with multiple active statement handlers by
creating new database handles as needed.

=head1 USAGE

The object constructor is
Bio::DB::GFF::Adaptor::dbi::caching_handle-E<gt>new().  This is called
like DBI-E<gt>connect() and takes the same arguments.  The returned object
looks and acts like a conventional database handle.

In addition to all the standard DBI handle methods, this package adds
the following:

=head2 dbi_quote

 Title   : dbi_quote
 Usage   : $string = $db->dbi_quote($sql,@args)
 Function: perform bind variable substitution
 Returns : query string
 Args    : the query string and bind arguments
 Status  : public

This method replaces the bind variable "?" in a SQL statement with
appropriately quoted bind arguments.  It is used internally to handle
drivers that don't support argument binding.

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

=head2 debug

 Title   : debug
 Usage   : $debug = $db->debug([$debug])
 Function: activate debugging messages
 Returns : current state of flag
 Args    : optional new setting of flag
 Status  : public

=cut

sub new {
  my $class    = shift;
  my @dbi_args = @_;
  my $self = bless {
		    dbh    => [],
		    args   => \@dbi_args,
		    debug => 0,
		   },$class;
  $self->dbh || $self->throw("Can't connect to database: " . DBI->errstr);
  $self;
}

sub AUTOLOAD {
  my($pack,$func_name) = $AUTOLOAD=~/(.+)::([^:]+)$/;
  return if $func_name eq 'DESTROY';
  my $self = shift or return DBI->$func_name(@_);
  $self->dbh->$func_name(@_);
}

sub debug {
  my $self = shift;
  my $d = $self->{debug};
  $self->{debug} = shift if @_;
  $d;
}

sub prepare {
  my $self  = shift;
  my $query = shift;

  # find a non-busy dbh
  my $dbh = $self->dbh || $self->throw("Can't connect to database: " . DBI->errstr);

  warn "Using prepare_cache\n" if $self->debug;
  my $sth = $dbh->prepare_cached($query, {}, 3) || $self->throw("Couldn't prepare query $query:\n ".DBI->errstr."\n");
  return $sth;
}

sub do_query {
  my $self = shift;
  my ($query,@args) = @_;
  warn $self->dbi_quote($query,@args),"\n" if $self->debug;
  my $sth = $self->prepare($query);
  $sth->execute(@args) || $self->throw("Couldn't execute query $query:\n ".DBI->errstr."\n");
  $sth;
}

sub dbh {
  my $self = shift;
  foreach (@{$self->{dbh}}) {
    return $_ if $_->inuse == 0;
  }
  # if we get here, we must create a new one
  warn "(Re)connecting to database\n" if $self->debug;
  my $dbh = DBI->connect(@{$self->{args}}) or return;

  $dbh->{PrintError} = 0;
  
  # for Oracle - to retrieve LOBs, need to define the length (Jul 15, 2002)
  $dbh->{LongReadLen} = 100*65535;
  $dbh->{LongTruncOk} = 0;
  $dbh->{mysql_auto_reconnect} = 1;

  my $wrapper = Bio::DB::GFF::Adaptor::dbi::faux_dbh->new($dbh);
  push @{$self->{dbh}},$wrapper;
  $wrapper;
}

# The clone method should only be called in child processes after a fork().
# It does two things: (1) it sets the "real" dbh's InactiveDestroy to 1,
# thereby preventing the database connection from being destroyed in
# the parent when the dbh's destructor is called; (2) it replaces the
# "real" dbh with the result of dbh->clone(), so that we now have an
# independent handle.
sub clone {
    my $self = shift;
    foreach (@{$self->{dbh}}) { $_->clone };
}

=head2 attribute

 Title   : attribute
 Usage   : $value = $db->attribute(AttributeName , [$newvalue])
 Function: get/set DBI::db handle attribute
 Returns : current state of the attribute
 Args    : name of the attribute and optional new setting of attribute
 Status  : public

  Under Bio::DB::GFF::Adaptor::dbi::caching_handle the DBI::db
  attributes that are usually set using hashref calls are unavailable.
  Use attribute() instead.  For example, instead of:

    $dbh->{AutoCommit} = 0;

  use

    $dbh->attribute(AutoCommit=>0);

=cut

sub attribute {
  my $self = shift;
  my $dbh = $self->dbh->{dbh};
  return $dbh->{$_[0]} = $_[1] if @_ == 2;
  return $dbh->{$_[0]}         if @_ == 1;
  return;
}

sub disconnect {
  my $self = shift;
  $_ && $_->disconnect foreach @{$self->{dbh}};
  $self->{dbh} = [];
}

sub dbi_quote {
  my $self = shift;
  my ($query,@args) = @_;
  my $dbh = $self->dbh;
  $query =~ s/\?/$dbh->quote(shift @args)/eg;
  $query;
}

package Bio::DB::GFF::Adaptor::dbi::faux_dbh;
use vars '$AUTOLOAD';

sub new {
  my $class = shift;
  my $dbh   = shift;
  bless {dbh=>$dbh},$class;
}

sub prepare {
  my $self = shift;
  my $sth = $self->{dbh}->prepare(@_) or return;
  $sth->{mysql_use_result} = 1 if $self->{dbh}->{Driver}{Name} eq 'mysql';
  $sth;
}

sub prepare_delayed {
  my $self = shift;
  my $sth = $self->{dbh}->prepare(@_) or return;
  $sth;
}

sub inuse {
    shift->{dbh}->{ActiveKids};
}

# The clone method should only be called in child processes after a fork().
# It does two things: (1) it sets the "real" dbh's InactiveDestroy to 1,
# thereby preventing the database connection from being destroyed in
# the parent when the dbh's destructor is called; (2) it replaces the
# "real" dbh with the result of dbh->clone(), so that we now have an
# independent handle.
sub clone {
    my $self = shift;
    $self->{dbh}{InactiveDestroy} = 1;
    $self->{dbh} = $self->{dbh}->clone;
}

sub DESTROY { }

sub AUTOLOAD {
  my($pack,$func_name) = $AUTOLOAD=~/(.+)::([^:]+)$/;
  return if $func_name eq 'DESTROY';
  my $self = shift;
  if( defined $self->{dbh} ) {
      $self->{dbh}->$func_name(@_);
  }
}

1;

__END__

=head1 BUGS

Report to the author.

=head1 SEE ALSO

L<DBI>, L<Bio::DB::GFF>, L<bioperl>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

