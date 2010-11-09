package Bio::DB::Tree::Store::DBI::SQLite;

use strict;
use DBI qw(:sql_types);
use Cwd 'abs_path';
use File::Spec;


sub table_definitions {
  my $self = shift;
  my $defs = 
      {
	  tree => <<END,
(
       tree_id integer primary autoincrement,
       label text not null,
       is_rooted integer default 1,
       root_id integer not null,
       annotations text
) ;

create unique index ui_tree_id on tree(tree_id);
create index i_tree_node_id ON tree(node_id);
create index i_treelabel on tree(label);
END

       node => <<END,
(
       node_id integer primary autoincrement,
       parent_id integer,
       label text,
       distance_to_parent real,
       annotations text,
       left_idx integer,
       right_idx integer,
);
create index ui_nodeid ON node (node_id);
create index i_nodelabel ON node (label);
create index i_leftidx ON node (left_idx);
create index i_rightidx ON node (right_idx);

END
};
  return $defs;
}

sub insert_node {
  my $self = shift;
  my $object = shift;
  my $index_flag = shift || 0;
  my $annotations = join(";",map { sprintf("%s=%s",$_,
					   join(",",$object->get_tag_values($_))) } 
			 $object->get_all_tags);
  my $sth = $self->_prepare(<<END);
INSERT INTO node (node_id,parent_id,label,distance_to_parent,annotations) VALUES (?,?,?)
END
  $sth->execute(undef,$object->parent_id,$object->id, $annotations);
  $sth->finish;
}

sub insert_tree {
  my $self = shift;
  my $object = shift;
  my $index_flag = shift || 0;
  my $annotations = join(";",map { sprintf("%s=%s",$_,
					   join(",",$object->get_tag_values($_))) } 
			 $object->get_all_tags);
  my $sth = $self->_prepare(<<END);
INSERT INTO node (tree_id,label,is_rooted,root_id,annotations) VALUES (?,?,?,?,?)
END
  $sth->execute(undef,$object->id,undef,$object->get_root_node,$annotations);
  $sth->finish;
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

sub _fetch_node {
  my $self = shift;
  my $id   = shift;

  my $sth = $self->_prepare(<<END);
SELECT (node_id, parent_id, label, distance_to_parent,annotations) FROM node where node_id = ?
END
  $sth->execute($id);
  my ($nid,$pid,$lbl,$distance,$annot) = @{$sth->fetchrow_arrayref};

  return Bio::DB::Tree::Node->new(-internal_id => $id,
				  -parent_id   => $pid,
				  -id          => $lbl,
				  -branch_lengths => $distance,
				  -annotations => $annot);
}

sub _fetch_node_children {
  my $self = shift;
  my $id   = shift;
  my $sortby = shift;

  my $sth = $self->_prepare(<<END);
SELECT node_id FROM node WHERE parent_id = ?
END
  $sth->execute($id);
  my @nodes;
  for my $nid ( @{$sth->fetchrow_arrayref}) {
    push @nodes, Bio::DB::Tree::Node->new(-node_id => $nid);
  }
  return @nodes;
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


sub _get_branch_length {
  my $self = shift;
  my $id   = shift;
  my $sth = $self->_prepare(<<END);
SELECT distance_to_parent FROM node WHERE node_id = ?
END
  my $d = @{$sth->fetchrow_array};
  $sth->finish;
  return $d;
}


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
    next if $t eq 'meta';      # don't get rid of meta data!
    my $table = $self->_qualify($t);
    $dbh->do("DROP table IF EXISTS $table") if $erase;
    my $query = "CREATE TABLE IF NOT EXISTS $table $tables->{$t}";
    $self->_create_table($dbh,$query);
  }
  1;
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
