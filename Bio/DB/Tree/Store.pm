#
# BioPerl module for Bio::DB::Tree::Store.pm
#
# Copyright Hilmar Lapp and Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Tree::Store.pm - The abstract interface for a local (or remote) persistent storage database. 

=head1 SYNOPSIS

  use Bio::DB::Tree::Store;

  # Open the tree database
  my $db = Bio::DB::Tree::Store->new(-adaptor => 'DBI::SQLite',
                                     -dsn     => '/path/to/database.db');

  my $id = 1;
  my $tree = $db->fetch($id);


=head1 DESCRIPTION

This is the generic interface for a persistent storage database for tree and node objects.  It is essentially an interface, and to use it you will need to provide the adaptor for the actual database. 

You can find the supported adaptors in the DBI sub-directory, such as
the SQLite adaptor (DBI::SQLite).

=head1 Using the an actual database adaptor

Create the connection to the actual storage by calling
Bio::DB::Tree::Store-E<gt>new() with the -adaptor parameter giving the
adaptor (e.g., -adaptor=-E<gt>'DBI::SQLite'), and any additional
parameters needed by that adaptor to connect to its database. See the
documentation for the adaptor to determine what those parameters
are. In general, for any DBI-compatible parameter you can expect at
least -dsn for specifying the database, and -dbuser and -dbpass for
authentication.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Hilmar Lapp, Jason Stajich

Email lapp-at-bioperl.org
Jason Stajich - jason-at-bioperl.org


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::DB::Tree::Store.pm;

use strict;

use base qw(Bio::Root::Root);

=head2 insert_tree

 Title   : insert_tree
 Usage   : $store->insert_tree($tree)
 Function: Inserts (adds) the given tree into the store.
 Example :
 Returns : True upon success, and false otherwise.
 Args    : The tree to be inserted, as a Bio::Tree::TreeI compliant object.

=cut

sub insert_tree{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 fetch_node

 Title   : fetch_node
 Usage   :
 Function: Fetch a tree node from the store.
 Example :
 Returns : A Bio::Tree::NodeI compliant object
 Args    : The primary key of the node to fetch.

=cut

sub fetch_node{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 max_id

 Title   : max_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub max_id{
   my ($self,@args) = @_;


}

sub optimize {
    my $self = shift;
    $self->throw_not_implemented();
}

sub index_tables {
    my $self = shift;
    $self->throw_not_implemented();
}

sub enable_keys  { }  # nullop
sub disable_keys { }  # nullop

