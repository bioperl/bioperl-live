
package Bio::DB::Tree::Store::DBI::SQLite;

use strict;
use Scalar::Util qw(blessed);
use DBI qw(:sql_types);
use Cwd 'abs_path';
use File::Spec;

use Data::Dumper;

use Bio::DB::Tree::Tree;
use Bio::DB::Tree::Node;
use base qw(Bio::DB::Tree::Store);


###
# object initialization
#
sub init {
  my $self          = shift;
  my ($dsn,
      $is_temporary,
      $autoindex,
      $dump_dir,
      $user,
      $pass,
      $dbi_options,
      $writeable,
      $create,
     ) = $self->_rearrange([qw(DSN TEMPORARY AUTOINDEX DUMPDIR
			       USER PASS OPTIONS WRITE
			       CREATE)
			   ],@_);
  $dbi_options  ||= {};
  $writeable    = 1 if $is_temporary or $dump_dir;
  $dsn or $self->throw("Usage: ".__PACKAGE__."->init(-dsn => \$dbh || \$dsn)");

  my $dbh;
  if (ref $dsn) {
    $dbh = $dsn;
  } else {
    $dsn = "dbi:SQLite:$dsn" unless $dsn =~ /^dbi:/;
    $dbh = DBI->connect($dsn,$user,$pass,$dbi_options) or $self->throw($DBI::errstr);
    $dbh->do("PRAGMA synchronous = OFF;"); # makes writes much faster
    $dbh->do("PRAGMA temp_store = MEMORY;"); # less disk I/O; some speedup
    $dbh->do("PRAGMA cache_size = 20000;"); # less disk I/O; some speedup
  }
  $self->dbh($dbh);
  $self->{is_temp}   = $is_temporary;
  $self->{writeable} = $writeable;

#  $self->default_settings;
#  $self->autoindex($autoindex)                   if defined $autoindex;
#  $self->dumpdir($dump_dir)                      if $dump_dir;
  if ($self->is_temp) {
    $self->init_tmp_database();
  } elsif ($create) {
    $self->init_database('erase');
  }
}

sub table_definitions {
  my $self = shift;
  my $defs =
    {
     tree => <<END,
(
       tree_id integer primary key autoincrement,
       label text not null,
       is_rooted integer default 1,
       root_id integer,
       score real,
       annotations text
);
create unique index ui_tree_id on tree(tree_id);
create index i_root_node_id ON tree(root_id);
create index i_treelabel on tree(label);
create index i_rooted on tree(is_rooted);
END

     node => <<END,
(
       node_id integer primary key autoincrement,
       parent_id integer,
       label text,
       distance_to_parent real,
       annotations text,
       left_idx integer,
       right_idx integer
);
create index ui_nodeid ON node (node_id);
create index i_nodelabel ON node (label);
create index i_leftidx ON node (left_idx);
create index i_rightidx ON node (right_idx);

END
    };
  return $defs;
}

=head2 Node methods

=cut

=head2 insert_node

 Title   : insert_node
 Usage   : $dbkey = $store->insert_node($nodeh, $parent);

 Function: Inserts a single node into the underlying storage. The
           default assumption is that the node does not exist yet,
           however the database defines (and/or enforces) that.

 Example :
 Returns : If successful, the primary key of the node that has been created. 

 Args :    The node to be inserted as a hashref or a Bio::Tree::NodeI
           object, and the parent as a primary key or a
           Bio::DB::Tree::Node object.

           If the node ia a hashref, the key '-parent' may be used to
           provide the primary key to the parent node instead of a
           second argument, and the key '-tree' to provide the primary
           key to the tree of which this node is part.  Depending on
           the backing database, the parent and/or the tree may be
           optional (and may be provided later through an update).

=cut

sub insert_node {
    my $self = shift;
    my ($nodeh, $parent) = @_;

    # unify the parent from the different ways it can be provided
    if (ref($parent)) {
        $parent = $parent->node_id();
    } elsif (!defined($parent)) {
        $parent = $nodeh->{'-parent'} if exists($nodeh->{'-parent'});
    }
    
    # unify between node object or hashref
    my $data = $self->_node_to_hash($nodeh);

    # store in database
    my $sth = $self->{sths}->{'insertNode'};
    if (! $sth) {
        $sth = $self->_prepare(
            "INSERT INTO node (parent_id,label,distance_to_parent,annotations)"
            ." VALUES (?,?,?,?)");
        $self->{sths}->{'insertNode'} = $sth;
    }
    $sth->execute($parent, 
                  $data->{'-id'}, $data->{'-branch_length'}, 
                  $data->{'-flatAnnotations'});
    my $pk = $self->dbh->func('last_insert_rowid');
    $nodeh->node_id($pk) 
        if $pk && blessed($nodeh) && $nodeh->isa("Bio::DB::Tree::Node");
    # done
    return $pk;
}

=head2 update_node

 Title   : update_node
 Usage   : $store->update_node($nodeh, $parent);

 Function: Updates a single node in the underlying storage. The node
           must exist already.

 Example :
 Returns : True if success and false otherwise.

 Args :    The node data to which the database is to be updated as a
           hashref or a Bio::DB::Tree::Node object that has been
           obtained from this store, and the parent as a primary key
           to the parent node or also as a Bio::DB::Tree::Node
           object. The parent is optional if it is not being updated.

           If the node is a hashref, the primary key is expected as
           the value for the '-pk' key, and the key '-parent' may be
           used to provide the primary key to the parent node instead
           of a second argument, and if present the key '-tree' can
           provide a new value (as primary key) for the tree of which
           this node is part.

=cut

sub update_node{
    my $self = shift;
    my ($nodeh, $parent) = @_;

    # unify the parent from the different ways it can be provided
    if (ref($parent)) {
        $parent = $parent->node_id();
    } elsif (!defined($parent)) {
        $parent = $nodeh->{'-parent'} if exists($nodeh->{'-parent'});
    }

    # unify the hashref and object options to a single format
    my $data = $self->_node_to_hash($nodeh);

    # store in database
    my $sth = $self->{sths}->{'updateNode'};
    if (! $sth) {
        $sth = $self->_prepare(
            "UPDATE node SET "
            ."parent_id = IFNULL(?,parent_id), "
            ."label = IFNULL(?,label), "
            ."distance_to_parent = IFNULL(?,distance_to_parent), "
            ."annotations = IFNULL(?,annotations) "
            ."WHERE node_id = ?");
        $self->{sths}->{'updateNode'} = $sth;
    }
    my $rv = $sth->execute($parent, 
                           $data->{'-id'}, $data->{'-branch_length'}, 
                           $data->{'-flatAnnotations'},
                           $data->{'-pk'});
    # done
    return $rv;
}

=head2 fetch_node

 Title   : fetch_node
 Usage   :
 Function: Fetch a tree node from the store.
 Example :
 Returns : A Bio::Tree::NodeI compliant object
 Args    : The primary key of the node to fetch.

=cut

sub fetch_node {
    my $self = shift;
    my $pk   = shift;

    my $sth = $self->{sths}->{'fetchNode'};
    if (! $sth) {
        $sth = $self->dbh->prepare(
            "SELECT node_id FROM node WHERE node_id = ?");
        $self->{sths}->{'fetchNode'} = $sth;
    }
    $sth->execute($pk);
    my $row = $sth->fetchrow_arrayref;
    return undef unless $row;
    return Bio::DB::Tree::Node->new(-node_id => $row->[0],
                                    -store   => $self);
}

=head2 populate_node

 Title   : populate_node
 Usage   :
 Function: Populates a node object's state from the store.
 Example :
 Returns : True on success and false otherwise.
 Args    : The Bio::DB::Tree::Node object to be populated.

=cut

sub populate_node {
    my $self = shift;
    my $node = shift;

    my $sth = $self->{sths}->{'populateNode'};
    if (! $sth) {
        $sth = $self->dbh->prepare(
            "SELECT parent_id,label,distance_to_parent "
            ."FROM node WHERE node_id = ?");
        $self->{sths}->{'populateNode'} = $sth;
    }
    $sth->execute($node->node_id);
    my $row = $sth->fetchrow_arrayref;
    return undef unless $row;
    my $dirty = $node->_dirty; # save so we can reset to current value
    $node->parent_id($row->[0]);
    $node->id($row->[1]);
    $node->branch_length($row->[2]);
    $node->_dirty($dirty); # restore dirty flag to where it was
    return 1;
}

=head2 fetch_nodes_by_parent

 Title   : fetch_nodes_by_parent
 Usage   : @children = $store->fetch_nodes_by_parent($parent);
 Function: Fetches the child nodes of a node from the database. 
 Example :
 Returns : An array of Bio::Tree::NodeI compliant objects
 Args :    The parent, either as a node object obtained from this store,
           or as its primary key.

=cut

sub fetch_nodes_by_parent{
    my $self   = shift;
    my $parent = shift;
    my $sortby = shift; # we aren't doing anything with this yet

    # unify parent
    $parent = $parent->node_id if ref($parent);

    my $sth = $self->{sths}->{'nodesByParent'};
    if (! $sth) {
        $sth = $self->_prepare("SELECT node_id FROM node WHERE parent_id = ?");
        $self->{sths}->{'nodesByParent'} = $sth;
    }
    $sth->execute($parent);

    my @nodes = ();
    while (my $row = $sth->fetchrow_arrayref) {
        push @nodes, Bio::DB::Tree::Node->new(-node_id => $row->[0],
                                              -store => $self);
    }
    return @nodes;
}

=head2 Tree methods

=cut

=head2 insert_tree

 Title   : insert_tree
 Usage   : $dbkey = $store->insert_tree($tree)
 Function: Inserts (adds) the given tree into the store.
 Example :
 Returns : If successful, the primary key to the newly created tree.

 Args :    The tree to be inserted, as a Bio::Tree::TreeI compliant
           object or a hashref.

=cut

sub insert_tree {
    my $self = shift;
    my $treeh = shift;

    # unify between node object or hashref
    my $data = $self->_tree_to_hash($treeh);

    # store in database
    my $sth = $self->{sths}->{'insertTree'};
    if (! $sth) {
        $sth = $self->_prepare(
            "INSERT INTO tree (label,is_rooted,root_id,score,annotations) "
            ."VALUES (?,?,?,?,?)");
        $self->{sths}->{'insertTree'} = $sth;
    }
    $sth->execute($data->{'-id'}, 
                  $data->{'-is_rooted'}, $data->{'-root'}, $data->{'-score'},
                  $data->{'-flatAnnotations'});
    my $pk = $self->dbh->func('last_insert_rowid');
    $treeh->tree_id($pk) 
        if $pk && blessed($treeh) && $treeh->isa("Bio::DB::Tree::Tree");
    # done
    return $pk;
}

=head2 fetch_tree

 Title   : fetch_tree
 Usage   : $tree = $store->fetch_tree($key);
 Function: Fetches a tree from the database. The tree is identified by
           its primary key.
 Example :
 Returns : A Bio::DB::Tree::Tree object if found and undef otherwise
 Args    : The primary key of the tree in the underlying storage.

=cut

sub fetch_tree{
    my $self = shift;
    my $pk = shift;

    my $sth = $self->{sths}->{'fetchTree'};
    if (! $sth) {
        $sth = $self->_prepare(
            "SELECT tree_id,label,is_rooted,root_id,score "
            ."FROM tree WHERE tree_id = ?");
        $self->{sths}->{'fetchTree'} = $sth;
    }
    $sth->execute($pk);
    my $row = $sth->fetchrow_arrayref;
    return undef unless $row;
    my $root;
    if ($row->[3]) {
        $root = Bio::DB::Tree::Node->new(-node_id => $row->[3],
                                         -store => $self);
    }
    return Bio::DB::Tree::Tree->new(-tree_id   => $row->[0],
                                    -store     => $self,
                                    -root_node => $root,
                                    -id        => $row->[1],
                                    -is_rooted => $row->[2],
                                    -score     => $row->[4]);
}

=head2 update_tree

 Title   : update_tree
 Usage   : $store->update_tree($treeh);

 Function: Updates a tree record in the underlying storage. The tree
           must exist already. This will not update the nodes of the tree.

 Example :
 Returns : True if success and false otherwise.

 Args :    The tree properties to which the database is to be updated as a
           hashref or a Bio::DB::Tree::Tree object that has been
           obtained from this store.

           If the argument is a hashref, the primary key is expected as
           the value for the '-pk' key, and the key '-root' may be
           used to provide the primary key to the root node.

=cut

sub update_tree{
    my $self = shift;
    my $treeh = shift;

    # unify the hashref and object options to a single format
    my $data = $self->_tree_to_hash($treeh);

    # update in database
    my $sth = $self->{sths}->{'updateTree'};
    if (! $sth) {
        $sth = $self->_prepare(
            "UPDATE tree SET "
            ."label = IFNULL(?,label), "
            ."root_id = IFNULL(?,root_id), "
            ."is_rooted = IFNULL(?,is_rooted), "
            ."annotations = IFNULL(?,annotations), "
            ."score = IFNULL(?,score) "
            ."WHERE tree_id = ?");
        $self->{sths}->{'updateTree'} = $sth;
    }
    my $rv = $sth->execute($data->{'-id'}, 
                           $data->{'-root'}, $data->{'-is_rooted'}, 
                           $data->{'-flatAnnotations'},
                           $data->{'-score'},
                           $data->{'-pk'});
    # done
    return $rv;
}

=head2 populate_tree

 Title   : populate_tree
 Usage   :
 Function: Populates a tree object's state from the store.
 Example :
 Returns : True on success and false otherwise.
 Args    : The Bio::DB::Tree::Tree object to be populated.

=cut

sub populate_tree {
    my $self = shift;
    my $tree = shift;

    my $sth = $self->{sths}->{'populateTree'};
    if (! $sth) {
        $sth = $self->_prepare(
            "SELECT label,is_rooted,root_id,score "
            ."FROM tree WHERE tree_id = ?");
        $self->{sths}->{'populateTree'} = $sth;
    }
    $sth->execute($tree->tree_id);
    my $row = $sth->fetchrow_arrayref;
    return undef unless $row;
    my $root;
    if ($row->[2]) {
        $root = Bio::DB::Tree::Node->new(-node_id => $row->[2],
                                         -store => $self);
    }
    my $dirty = $tree->_dirty; # save so we can reset to current value
    $tree->root($root);
    $tree->id($row->[0]);
    $tree->rooted($row->[1]);
    $tree->score($row->[3]);
    $tree->_dirty($dirty); # restore flag to where it was
    return 1;
}


sub _set_parent_ids {
  my $self = shift @_;
  my $parent_id = shift @_;
  my @nodes = @_;
  my $sth= $self->_prepare(<<END);
UPDATE node SET parent_id = ? WHERE node_id = ?
END

  for my $node ( @nodes ) {
    $sth->execute($parent_id,$node);
  }
  $sth->finish;
}


sub _fetch_node_all_children {
  my $self = shift;
  my $id   = shift;
  my $sortby = shift;

  my $sth = $self->_prepare(<<END);
SELECT node_id FROM node WHERE left_idx <= ? AND right_idx >= ?
END
  $sth->execute($id);
  my @nodes;
  for my $nid ( @{$sth->fetchrow_arrayref}) {
    push @nodes, Bio::DB::Tree::Node->new(-node_id => $nid);
  }
  return @nodes;
}

sub _unset_node_parent {
  my $self = shift;
  my $id   = shift;
  my $sth = $self->_prepare(<<END);
UPDATE node set parent_id = NULL where node_id = ?
END
  $sth->execute($id);
  $sth->finish;
}

sub _delete_node {
  my $self = shift;
  my $id   = shift;

  my $sth = $self->_prepare(<<END);
DELETE FROM node self, node n WHERE n.left_idx >= s.left_idx AND n.right_idx <= s.right_idx AND n.node_id = ?
END
  $sth->execute($id);
  $sth->finish;
#  $self->dbh->commit;
}


sub _is_Leaf {
  my $self = shift;
  my $id = shift;
  my $sth = $self->_prepare(<<END);
SELECT COUNT(*) FROM node WHERE node_id = ? AND left_idx == right_idx
END
  my ($isleaf) = @{$sth->fetchrow_arrayref};
  return $isleaf;
}

=head2 DB management methods

=cut

sub max_id {
  my $self = shift;
  my $type = shift;
  my $tbl;
  if ( $type =~ /node/i ) {
    $tbl = $self->_node_table;
  } elsif ( $type =~ /tree/i ) {
    $tbl = $self->_tree_table;
  } else {
    $self->throw("cannot call max id without a type");
  }

  my $sth  = $self->_prepare("SELECT max(id) from $tbl");
  $sth->execute or $self->throw($sth->errstr);
  my ($id) = $sth->fetchrow_array;
  $id;
}

sub _init_database {
  my $self  = shift;
  my $erase = shift;

  my $dbh    = $self->dbh;
  my $tables = $self->table_definitions;
  for my $t (keys %$tables) {
    $dbh->do("DROP table IF EXISTS $t") if $erase;
    my $query = "CREATE TABLE IF NOT EXISTS $t $tables->{$t}";
    $self->_create_table($dbh,$query);
  }
  1;
}


sub _create_table {
    my $self         = shift;
    my ($dbh,$query) = @_;
    for my $q (split ';',$query) {
	chomp($q);
	next unless $q =~ /\S/;
	$dbh->do("$q;\n") or $self->throw("$q -- ".$dbh->errstr);
    }
}

sub init_tmp_database {
  my $self = shift;

  my $dbh    = $self->dbh;
  my $tables = $self->table_definitions;

  for my $t (keys %$tables) {
    next if $t eq 'meta';  # done earlier
    my $table = $self->_qualify($t);
    my $query = "CREATE TEMPORARY TABLE $table $tables->{$t}";
    $self->_create_table($dbh,$query);
  }
  1;
}

sub optimize {
  my $self = shift;
  $self->dbh->do("VACUUM TABLE $_") foreach $self->index_tables;
}

sub index_tables {
  my $self = shift;
  return qw(node tree)
}

sub _enable_keys  { }  # nullop
sub _disable_keys { }  # nullop

=head1 Private methods

Don't call from outside.

=head2 _node_to_hash

 Title   : _node_to_hash
 Usage   : $nodeHash = $self->_node_to_hash($nodeObj);
 Function: Unifies a hashref and node object representation to a
           hashref representation. This will also flatten out the
           key-value annotations.

 Example :
 Returns : A hashref with the standard keys set if the input object
           had a value for them, and with the flattened annotations
           under the -flatAnnotations key if the input object had
           annotations. If the input object was a hashref, the
           returned ref will have at least all keys that the input ref
           had.
 Args :    A hashref of node attribute names and their values, or a
           Bio::Tree::NodeI object.

=cut

sub _node_to_hash {
    my $self = shift;
    my $obj = shift;

    # unify the hashref and object options to a single format
    my $h;
    if (ref($obj) eq "HASH") {
        $h = {%$obj};
        # flatten out all key-value annotations (if we have any)
        my $flatAnnots = _flatten_keyvalues($obj->{'-annotations'});
        $h->{'-flatAnnotations'} = $flatAnnots if $flatAnnots;
    } elsif (blessed($obj) && $obj->isa("Bio::Tree::NodeI")) {
        $h = {};
        # flatten out all key-value annotations (if we have any)
        my $flatAnnots = _flatten_annotations($obj);
        $h->{'-flatAnnotations'} = $flatAnnots if $flatAnnots;
        $h->{'-id'} = $obj->id if $obj->id;
        $h->{'-branch_length'} = $obj->branch_length if $obj->branch_length();
        if ($obj->isa("Bio::DB::Tree::Node")) {
            $h->{'-pk'} = $obj->node_id;
            $h->{'-parent'} = $obj->parent_id if $obj->parent_id;
        }
    } else {
        $self->throw("don't know how to deal with ".$obj); 
    }
    return $h;
}

=head2 _tree_to_hash

 Title   : _tree_to_hash
 Usage   : $treeHash = $self->_tree_to_hash($treeObj);
 Function: Unifies a hashref and tree object representation to a
           hashref representation. This will also flatten out the
           key-value annotations.

 Example :
 Returns : A hashref with the standard keys set if the input object
           had a value for them, and with the flattened annotations
           under the -flatAnnotations key if the input object had
           annotations. If the input object was a hashref, the
           returned ref will have at least all keys that the input ref
           had.
 Args :    A hashref of tree attribute names and their values, or a
           Bio::Tree::TreeI object.

=cut

sub _tree_to_hash {
    my $self = shift;
    my $obj = shift;

    # unify the hashref and object options to a single format
    my $h;
    if (ref($obj) eq "HASH") {
        $h = {%$obj};
        # flatten out all key-value annotations (if we have any)
        my $flatAnnots = _flatten_keyvalues($obj->{'-annotations'});
        $h->{'-flatAnnotations'} = $flatAnnots if $flatAnnots;
    } elsif (blessed($obj) && $obj->isa("Bio::Tree::TreeI")) {
        $h = {};
        # flatten out all key-value annotations (if we have any)
        my $flatAnnots = _flatten_annotations($obj);
        $h->{'-flatAnnotations'} = $flatAnnots if $flatAnnots;
        $h->{'-id'} = $obj->id if $obj->id;
        $h->{'-is_rooted'} = $obj->is_rooted;
        $h->{'-score'} = $obj->score if $obj->score;
        if ($obj->isa("Bio::DB::Tree::Tree")) {
            $h->{'-pk'} = $obj->tree_id;
            $h->{'-root'} = $obj->root->node_id if $obj->root;
        }
    } else {
        $self->throw("don't know how to deal with ".$obj);         
    }
    return $h;
}

##
# flatten out into a string the tag/value annotations of a Bio::AnnotatableI
##
sub _flatten_annotations {
  my $annHolder = shift;
  my @tags = $annHolder->get_all_tags;
  my $flatAnnot;
  if (@tags) {
    $flatAnnot = 
      join(";",
	   map { sprintf("%s=%s",$_,
			 join(",",$annHolder->get_tag_values($_)))
	       } 
	   @tags);
  }
  return $flatAnnot;
}

##
# flatten out into a string the key/value pairs of a hash
##
sub _flatten_keyvalues {
  my $h = shift;
  my $flatAnnot;
  if ($h && keys(%$h)) {
    $flatAnnot = join(";",
		      map { 
			sprintf("%s=%s",$_,join(",",$h->{$_}))
		      } 
		      keys(%$h));
  }
  return $flatAnnot;
}

=head1 NAME

Bio::DB::Tree::Store::DBI::SQLite -

=head1 SYNOPSIS

  use Bio::DB::Tree::Store;

  # Open the tree database
  my $db = Bio::DB::Tree::Store->new(-adaptor => 'DBI::SQLite',
                                     -dsn     => '/path/to/database.db');

  my $id = 1;
  my $tree = $db->fetch($id);


=head1 DESCRIPTION

Bio::DB::Tree::Store::SQLite is the SQLite adaptor for
Bio::DB::Tree::Store. You will not create it directly, but
instead use Bio::DB::Tree::Store-E<gt>new() to do so.

See L<Bio::DB::Tree::Store> for complete usage instructions.

=head2 Using the SQLite adaptor

To establish a connection to the database, call
Bio::DB::Tree::Store-E<gt>new(-adaptor=E<gt>'DBI::SQLite',@more_args). The
additional arguments are as follows:


=head1 AUTHOR

Jason Stajich - jason-at-bioperl.org

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself. See the BioPerl license for
more details.


=cut

1;
