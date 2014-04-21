#
# BioPerl module for Bio::Tree::TreeFunctionsI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code


=head1 NAME

Bio::Tree::TreeFunctionsI - Decorated Interface implementing basic Tree exploration methods

=head1 SYNOPSIS

  use Bio::TreeIO;
  my $in = Bio::TreeIO->new(-format => 'newick', -file => 'tree.tre');

  my $tree = $in->next_tree;

  my @nodes = $tree->find_node('id1');

  if( $tree->is_monophyletic(-nodes => \@nodes, -outgroup => $outnode) ){
   #...
  }

=head1 DESCRIPTION

This interface provides a set of implementated Tree functions which
only use the defined methods in the TreeI or NodeI interface.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich, Aaron Mackey, Justin Reese

Email jason-at-bioperl-dot-org
Email amackey-at-virginia.edu
Email jtr4v-at-virginia.edu

=head1 CONTRIBUTORS

Sendu Bala, bix@sendu.me.uk

Rerooting code was worked on by

  Daniel Barker d.barker-at-reading.ac.uk
  Ramiro Barrantes Ramiro.Barrantes-at-uvm.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::Tree::TreeFunctionsI;

use strict;
use base qw(Bio::Tree::TreeI);


=head2 find_node

 Title   : find_node
 Usage   : my @nodes = $self->find_node(-id => 'node1');
 Function: returns all nodes that match a specific field, by default this
           is id, but different branch_length, 
 Returns : List of nodes which matched search
 Args    : text string to search for
           OR
           -fieldname => $textstring

=cut

sub find_node {
   my ($self, $type, $field) = @_;
   if( ! defined $type ) { 
       $self->warn("Must request a either a string or field and string when searching");
   }

   # all this work for a '-' named field
   # is so that we could potentially 
   # expand to other constraints in 
   # different implementations
   # like 'find all nodes with boostrap < XX'

   if( ! defined $field ) { 
       # only 1 argument, default to searching by id
       $field = $type; 
       $type  = 'id';
   } else {   
       $type =~ s/^-//;
   }

   # could actually do this by testing $rootnode->can($type) but
   # it is possible that a tree is implemeted with different node types
   # - although it is unlikely that the root node would be richer than the
   # leaf nodes.  Can't handle NHX tags right now

   my @nodes = grep { $_->can($type) && defined $_->$type() &&
                      $_->$type() eq $field } $self->get_nodes();

   if ( wantarray) { 
       return @nodes;
   } else { 
       if( @nodes > 1 ) { 
           $self->warn("More than 1 node found but caller requested scalar, only returning first node");
       }
       return shift @nodes;
   }
}


=head2 remove_Node

 Title   : remove_Node
 Usage   : $tree->remove_Node($node)
 Function: Removes a node from the tree
 Returns : boolean represent status of success
 Args    : either Bio::Tree::NodeI or string of the node id

=cut

sub remove_Node {
   my ($self,$input) = @_;
   my $node = undef;
   unless( ref($input) ) {
       $node = $self->find_node($input);
   }  elsif( ! $input->isa('Bio::Tree::NodeI') ) {
       $self->warn("Did not provide either a valid Bio::Tree::NodeI object or id to remove_node");
       return 0;
   } else { 
       $node = $input;
   }
   if( ! $node->ancestor && 
       $self->get_root_node->internal_id != $node->internal_id) {
     $self->warn("Node (".$node->to_string . ") has no ancestor, can't remove!");
   } else { 
     $node->ancestor->remove_Descendent($node);
   }
}


=head2 get_lineage_nodes

 Title   : get_lineage_nodes
 Usage   : my @nodes = $tree->get_lineage_nodes($node);
 Function: Given a node or its ID, get its full lineage, i.e. all its ancestors,
           from the root to the most recent ancestor. Only use the node ID as
           input if the nodes have been added to the tree.
 Returns : list of nodes
 Args    : either Bio::Tree::NodeI (or string of the node id)

=cut

sub get_lineage_nodes {
    my ($self, $input) = @_;
    my $node;

    # Sanity checks
    if (ref $input) {
        if (not $input->isa('Bio::Tree::NodeI')) {
            $self->throw("Did not provide a valid Bio::Tree::NodeI object or ID string to get_lineage_nodes");
        }
        $node = $input;
    } else {
        $node = $self->find_node($input);
    }

    # When dealing with Bio::Taxon objects with databases, the root will always
    # be the database's root, ignoring this Tree's set root node; prefer the
    # Tree's idea of root.
    my $root = $self->get_root_node || '';

    my @lineage;
    while ($node) {
        $node = $node->ancestor || last;
        unshift(@lineage, $node);
        $node eq $root && last;
    }
    return @lineage;
}


=head2 get_lineage_string

 Title   : get_lineage_string
 Usage   : my $lineage = $tree->get_lineage_string($node);
 Function: Get the string representation of the full lineage of a node, e.g.
           for the Enterobacteriales node, return
           Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales.
           This method uses get_lineage_nodes internally and therefore inherits
           of all of its caveats.
 Returns : string
 Args    : * either Bio::Tree::NodeI (or string of the node id)
           * an optional separator (default: ';')

=cut

sub get_lineage_string {
    my ($self, $input, $sep) = @_;
    $sep ||= ';';
    my $node;
    unless (ref $input) {
        $node = $self->find_node($input);
    }
    elsif (! $input->isa('Bio::Tree::NodeI')) {
        $self->warn("Did not provide either a valid Bio::Tree::NodeI object or id to get_lineage_nodes");
        return;
    }
    else {
        $node = $input;
    }
    my @nodes = ($self->get_lineage_nodes($node), $node);
    for my $i (0 .. scalar @nodes - 1) {
        my $node_name = $nodes[$i]->node_name || '';
        if ($node_name =~ m/$sep/) {
           $self->warn("Separator '$sep' is not safe to use because the node ".
               "called '$node_name' contains it. Consider using another separator".
               " or sanitizing the node name.");
        }
        $nodes[$i] = $node_name;
    }
    return join $sep, @nodes;
}


=head2 splice

 Title   : splice
 Usage   : $tree->splice(-remove_id => \@ids);
 Function: Remove all the nodes from a tree that correspond to the supplied
           args, making all the descendents of a removed node the descendents
           of the removed node's ancestor.
           You can ask to explicitly remove certain nodes by using -remove_*,
           remove them conditionally by using -remove_* in combination with
           -keep_*, or remove everything except certain nodes by using only
           -keep_*.
 Returns : n/a
 Args    : just a list of Bio::Tree::NodeI objects to remove, OR
           -key => value pairs, where -key has the prefix 'remove' or 'keep',
           followed by an underscore, followed by a fieldname (like for the
           method find_node). Value should be a scalar or an array ref of
           scalars (again, like you might supply to find_node).

           So (-remove_id => [1, 2]) will remove all nodes from the tree that
           have an id() of '1' or '2', while
           (-remove_id => [1, 2], -keep_id => [2]) will remove all nodes with
           an id() of '1'.
           (-keep_id => [2]) will remove all nodes unless they have an id() of
           '2' (note, no -remove_*).

           -preserve_lengths => 1 : setting this argument will splice out
           intermediate nodes, preserving the original total length between
           the ancestor and the descendants of the spliced node. Undef 
           by default.

=cut

sub splice {
    my ($self, @args) = @_;
    $self->throw("Must supply some arguments") unless @args > 0;
    my $preserve_lengths = 0;
    my @nodes_to_remove;
    if (ref($args[0])) {
        $self->throw("When supplying just a list of Nodes, they must be Bio::Tree::NodeI objects") unless $args[0]->isa('Bio::Tree::NodeI');
        @nodes_to_remove = @args;
    }
    else {
        $self->throw("When supplying -key => value pairs, must be an even number of args") unless @args % 2 == 0;
        my %args = @args;
        my @keep_nodes;
        my @remove_nodes;
        my $remove_all = 1;
        while (my ($key, $value) = each %args) {
            my @values = ref($value) ? @{$value} : ($value);

            if ($key =~ s/remove_//) {
                $remove_all = 0;
                foreach my $value (@values) {
                    push(@remove_nodes, $self->find_node($key => $value));
                }
            }
            elsif ($key =~ s/keep_//) {
                foreach my $value (@values) {
                    push(@keep_nodes, $self->find_node($key => $value));
                }
            }
            elsif ($key =~ /preserve/) {
                $preserve_lengths = $value;
            }
        }

        if ($remove_all) {
            if (@keep_nodes == 0) {
                $self->warn("Requested to remove everything except certain nodes, but those nodes were not found; doing nothing instead");
                return;
            }

            @remove_nodes = $self->get_nodes;
        }
        if (@keep_nodes > 0) {
            my %keep_iids = map { $_->internal_id => 1 } @keep_nodes;
            foreach my $node (@remove_nodes) {
                push(@nodes_to_remove, $node) unless exists $keep_iids{$node->internal_id};
            }
        }
        else {
            @nodes_to_remove = @remove_nodes;
        }
    }
    # do the splicing
    #*** the algorithm here hasn't really been thought through and tested much,
    #    will probably need revising
    my %root_descs;
    my $reroot = 0;
    foreach my $node (@nodes_to_remove) {
        my @descs = $node->each_Descendent;
        my $ancestor = $node->ancestor;
        if (! $ancestor && ! $reroot) {
            # we're going to remove the tree root, so will have to re-root the
            # tree later
            $reroot = 1;
            %root_descs = map { $_->internal_id => $_ } @descs;
            $node->remove_all_Descendents;
            next;
        }
        if (exists $root_descs{$node->internal_id}) {
            # well, this one can't be the future root anymore
            delete $root_descs{$node->internal_id};
            # but maybe one of this one's descs will become the root
            foreach my $desc (@descs) {
                $root_descs{$desc->internal_id} = $desc;
            }
        }
        # make the ancestor of our descendents our own ancestor, and give us
        # no ancestor of our own to remove us from the tree
        foreach my $desc (@descs) {
            $desc->ancestor($ancestor);
            $desc->branch_length($desc->branch_length + $node->branch_length) if $preserve_lengths;
        }
        $node->ancestor(undef);
    }
    if ($reroot) {
        my @candidates = values %root_descs;
        $self->throw("After splicing, there was no tree root!") unless @candidates > 0;
        $self->throw("After splicing, the original root was removed but there are multiple candidates for the new root!") unless @candidates == 1;
        $self->set_root_node($candidates[0]); # not sure its valid to use the reroot() method
    }
}


=head2 get_lca

 Title   : get_lca
 Usage   : get_lca(-nodes => \@nodes ); OR
           get_lca(@nodes);
 Function: given two or more nodes, returns the lowest common ancestor (aka most
           recent common ancestor)
 Returns : node object or undef if there is no common ancestor
 Args    : -nodes => arrayref of nodes to test, OR
           just a list of nodes

=cut

sub get_lca {
    my ($self, @args) = @_;
    my ($nodes) = $self->_rearrange([qw(NODES)],@args);
    my @nodes;
    if (ref($nodes) eq 'ARRAY') {
        @nodes = @{$nodes};
    }
    else {
        @nodes = @args;
    }
    @nodes >= 2 or $self->throw("At least 2 nodes are required");
    # We must go root->leaf to get the correct answer to lca (in a world where
    # internal_id might not be uniquely assigned), but leaf->root is more
    # forgiving (eg. lineages may not all have the same root, or they may have
    # different numbers of 'minor' taxa inbeteen 'major' ones).
    #
    # I use root->leaf so that we can easily do multiple nodes at once - no
    # matter what taxa are below the lca, the lca and all its ancestors ought to
    # be identical.
    my @paths;
    foreach my $node (@nodes) {
        unless(ref($node) && $node->isa('Bio::Tree::NodeI')) {
            $self->throw("Cannot process get_lca() with a non-NodeI object ($node)\n");
        }
        my @path = ($self->get_lineage_nodes($node), $node);
        push(@paths, \@path);
    }
    return unless @paths >= 2;
    my $lca;
    LEVEL: while ($paths[0] > 0) {
        my %node_ids;
        my $node;
        foreach my $path (@paths) {
            $node = shift(@{$path}) || last LEVEL;
            my $node_id = $node->internal_id;
            unless (defined $node_id) {
                $self->warn("One of the lineages had a node with no internal_id, can't calculate the common ancestor");
                return;
            }
            $node_ids{$node_id}++;
        }
        if (keys %node_ids == 1) {
            $lca = $node;
        }
        else {
            # at this point in the lineage the nodes are different; the previous
            # loop had the lca
            last LEVEL;
        }
    }
    # If the tree that we are contains the lca (get_lca could have been called
    # on an empty tree, since it works with plain Nodes), prefer to return the
    # node object that belongs to us
    if ($lca && $self->number_nodes > 0) {
        my $own_lca = $self->find_node(-internal_id => $lca->internal_id);
        $lca = $own_lca if $own_lca;
    }
    return $lca;
}


=head2 merge_lineage

 Title   : merge_lineage
 Usage   : merge_lineage($node)
 Function: Merge a lineage of nodes with this tree.
 Returns : true for success, false (and a warning) otherwise
 Args    : Bio::Tree::TreeI with only one leaf, OR
           Bio::Tree::NodeI which has an ancestor

 For example, if we are the tree $tree:

 +---B
 |
 A
 |
 +---C

 and we want to merge the lineage $other_tree:

 A---C---D

 After calling $tree->merge_lineage($other_tree), $tree looks like:

 +---B
 |
 A
 |
 +---C---D

=cut

sub merge_lineage {
    my ($self, $thing) = @_;
    $self->throw("Must supply an object reference") unless ref($thing);

    my $lineage_leaf;
    if ($thing->isa('Bio::Tree::TreeI')) {
        my @leaves = $thing->get_leaf_nodes;
        $self->throw("The supplied Tree can only have one leaf") unless @leaves == 1;
        $lineage_leaf = shift(@leaves);
    }
    elsif ($thing->isa('Bio::Tree::NodeI')) {
        $self->throw("The supplied Node must have an ancestor") unless $thing->ancestor;
        $lineage_leaf = $thing;
    }

    # Find the lowest node in the supplied lineage that is in the tree
    # That will be our lca and we can merge at the node below
    my @lineage = ($lineage_leaf, reverse($self->get_lineage_nodes($lineage_leaf)));
    my $merged = 0;
    my $node;
    my $i = 0;
    while ($i <= $#lineage) {
        $node = $self->find_node(-internal_id => $lineage[$i]->internal_id);
        if (defined $node) {
            $merged = 1;
            last;
        }
        $i++;
    }
    if (not $merged) {
        $self->warn("Could not merge the lineage of ".$lineage_leaf->id." with the rest of the tree");
    }

    # Merge descendents, recursively
    while ($i > 0) {
        $node->add_Descendent($lineage[$i-1]);
        $node = $self->find_node(-internal_id => $lineage[$i-1]->internal_id);
        $i--;
    }

    return $merged;
}


=head2 contract_linear_paths

 Title   : contract_linear_paths
 Usage   : contract_linear_paths()
 Function: Splices out all nodes in the tree that have an ancestor and only one
           descendent.
 Returns : n/a
 Args    : none for normal behaviour, true to dis-regard the ancestor requirment
           and re-root the tree as necessary

 For example, if we are the tree $tree:

             +---E
             |
 A---B---C---D
             |
             +---F

 After calling $tree->contract_linear_paths(), $tree looks like:

     +---E
     |
 A---D
     |
     +---F

 Instead, $tree->contract_linear_paths(1) would have given:

 +---E
 |
 D
 |
 +---F

=cut

sub contract_linear_paths {
    my $self = shift;
    my $reroot = shift;
    my @remove;
    foreach my $node ($self->get_nodes) {
        if ($node->ancestor && $node->each_Descendent == 1) {
            push(@remove, $node);
        }
    }
    $self->splice(@remove) if @remove;
    if ($reroot) {
        my $root = $self->get_root_node;
        my @descs = $root->each_Descendent;
        if (@descs == 1) {
            my $new_root = shift(@descs);
            $self->set_root_node($new_root);
            $new_root->ancestor(undef);
        }
    }
}


=head2 is_binary

  Example    : is_binary(); is_binary($node);
  Description: Finds if the tree or subtree defined by
               the internal node is a true binary tree
               without polytomies
  Returns    : boolean
  Exceptions : 
  Args       : Internal node Bio::Tree::NodeI, optional


=cut

sub is_binary {
    my $self = shift;
    my $node = shift || $self->get_root_node;

    my $binary = 1;
    my @descs = $node->each_Descendent;
    $binary = 0 unless @descs == 2 or @descs == 0;
    #print "$binary, ", scalar @descs, "\n";

    # recurse
    foreach my $desc (@descs) {
        $binary += $self->is_binary($desc) -1;
    }
    $binary = 0 if $binary < 0;
    return $binary;
}


=head2 force_binary

 Title   : force_binary
 Usage   : force_binary()
 Function: Forces the tree into a binary tree, splitting branches arbitrarily
           and creating extra nodes as necessary, such that all nodes have
           exactly two or zero descendants.
 Returns : n/a
 Args    : none

 For example, if we are the tree $tree:

 +---G
 |
 +---F
 |
 +---E
 |
 A
 |
 +---D
 |
 +---C
 |
 +---B

 (A has 6 descendants B-G)

 After calling $tree->force_binary(), $tree looks like:

         +---X
         |
     +---X
     |   |
     |   +---X
     |
 +---X
 |   |
 |   |   +---G
 |   |   |
 |   +---X
 |       |
 |       +---F
 A
 |       +---E
 |       |
 |   +---X
 |   |   |
 |   |   +---D
 |   |
 +---X
     |
     |   +---C
     |   |
     +---X
         |
         +---B

 (Where X are artificially created nodes with ids 'artificial_n', where n is
 an integer making the id unique within the tree)

=cut

sub force_binary {
    my $self = shift;
    my $node = shift || $self->get_root_node;

    my @descs = $node->each_Descendent;
    if (@descs > 2) {
        # Removed overly verbose warning - cjfields 3-12-11
        
        # Many nodes have no identifying names, a simple warning is probably
        # enough.

        $self->warn("Node has more than two descendants\nWill do an arbitrary balanced split");
        my @working = @descs;
        # create an even set of artifical nodes on which to later hang the descs
        my $half = @working / 2;
        $half++ if $half > int($half);
        $half = int($half);
        my @artificials;
        while ($half > 1) {
            my @this_level;
            foreach my $top_node (@artificials || $node) {
                for (1..2) {
                    my $art = $top_node->new(-id => "artificial_".++$self->{_art_num});
                    $top_node->add_Descendent($art);
                    push(@this_level, $art);
                }
            }
            @artificials = @this_level;
            $half--;
        }
        # attach two descs to each artifical leaf
        foreach my $art (@artificials) {
            for (1..2) {
                my $desc = shift(@working) || $node->new(-id => "artificial_".++$self->{_art_num});
                $desc->ancestor($art);
            }
        }
    }
    elsif (@descs == 1) {
        # ensure that all nodes have 2 descs
        $node->add_Descendent($node->new(-id => "artificial_".++$self->{_art_num}));
    }
    # recurse
    foreach my $desc (@descs) {
        $self->force_binary($desc);
    }
}


=head2 simplify_to_leaves_string

 Title   : simplify_to_leaves_string
 Usage   : my $leaves_string = $tree->simplify_to_leaves_string()
 Function: Creates a simple textual representation of the relationship between
           leaves in self. It forces the tree to be binary, so the result may
           not strictly correspond to the tree (if the tree wasn't binary), but
           will be as close as possible. The tree object is not altered. Only
           leaf node ids are output, in a newick-like format.
 Returns : string
 Args    : none

=cut

sub simplify_to_leaves_string {
    my $self = shift;

    # Before contracting and forcing binary we need to clone self, but Clone.pm
    # clone() seg faults and fails to make the clone, whilst Storable dclone
    # needs $self->{_root_cleanup_methods} deleted (code ref) and seg faults at
    # end of script. Let's make our own clone...
    my $tree = $self->_clone;

    $tree->contract_linear_paths(1);
    $tree->force_binary;
    foreach my $node ($tree->get_nodes) {
        my $id = $node->id;
        $id = ($node->is_Leaf && $id !~ /^artificial/) ? $id : '';
        $node->id($id);
    }

    my %paired;
    my @data = $self->_simplify_helper($tree->get_root_node, \%paired);

    return join(',', @data);
}


# alias
sub _clone { shift->clone(@_) }


# safe node clone that doesn't seg fault, but deliberately loses ancestors and
# descendents
sub _clone_node {
    my ($self, $node) = @_;
    my $clone = $node->new;

    while (my ($key, $val) = each %{$node}) {
        if ($key eq '_desc' || $key eq '_ancestor') {
            next;
        }
        ${$clone}{$key} = $val;
    }

    return $clone;
}


# tree string generator for simplify_to_leaves_string, based on
# Bio::TreeIO::newick::_write_tree_Helper
sub _simplify_helper {
    my ($self, $node, $paired) = @_;
    return () if (!defined $node);

    my @data = ();
    foreach my $node ($node->each_Descendent()) {
        push(@data, $self->_simplify_helper($node, $paired));
    }

    my $id = $node->id_output || '';
    if (@data) {
        unless (exists ${$paired}{"@data"} || @data == 1)  {
            $data[0] = "(" . $data[0];
            $data[-1] .= ")";
            ${$paired}{"@data"} = 1;
        }
    }
    elsif ($id) {
        push(@data, $id);
    }

    return @data;
}


=head2 distance

 Title   : distance
 Usage   : distance(-nodes => \@nodes )
 Function: returns the distance between two given nodes
 Returns : numerical distance
 Args    : -nodes => arrayref of nodes to test
           or ($node1, $node2)

=cut

sub distance {
    my ($self,@args) = @_;
    my ($nodes) = $self->_rearrange([qw(NODES)],@args);
    if( ! defined $nodes ) {
        $self->warn("Must supply two nodes or -nodes parameter to distance() method");
        return;
    }
    elsif (ref($nodes) eq 'ARRAY') {
        1;
    }
    elsif ( @args == 2) { # assume these are nodes...
            $nodes = \@args;
    }
    else {
        $self->warn("Must supply two nodes or -nodes parameter to distance() method");
        return;
    }
    $self->throw("Must provide 2 nodes") unless @{$nodes} == 2;

    my $lca = $self->get_lca(@{$nodes});
    unless($lca) { 
        $self->warn("could not find the lca of supplied nodes; can't find distance either");
        return;
    }

    my $cumul_dist = 0;
    my $warned = 0;
    foreach my $current_node (@{$nodes}) {
        while (1) {
            last if $current_node eq $lca;
            if ($current_node->branch_length) {
                $cumul_dist += $current_node->branch_length;
            }
            elsif (! $warned) {
                $self->warn("At least some nodes do not have a branch length, the distance returned could be wrong");
                $warned = 1;
            }

            $current_node = $current_node->ancestor || last;
        }
    }

    return $cumul_dist;
}


=head2 is_monophyletic

 Title   : is_monophyletic
 Usage   : if( $tree->is_monophyletic(-nodes => \@nodes, 
                                      -outgroup => $outgroup)
 Function: Will do a test of monophyly for the nodes specified
           in comparison to a chosen outgroup
 Returns : boolean
 Args    : -nodes    => arrayref of nodes to test
           -outgroup => outgroup to serve as a reference

=cut

sub is_monophyletic{
   my ($self,@args) = @_;
   my ($nodes,$outgroup) = $self->_rearrange([qw(NODES OUTGROUP)],@args);

   if( ! defined $nodes || ! defined $outgroup ) {
       $self->warn("Must supply -nodes and -outgroup parameters to the method
is_monophyletic");
       return;
   }
   if( ref($nodes) !~ /ARRAY/i ) {
       $self->warn("Must provide a valid array reference for -nodes");
   }

   my $clade_root = $self->get_lca(@{$nodes});
   unless( defined $clade_root ) { 
       $self->warn("could not find clade root via lca");
       return;
   }

   my $og_ancestor = $outgroup->ancestor;
   while( defined ($og_ancestor ) ) {
       if( $og_ancestor->internal_id == $clade_root->internal_id ) {
           # monophyly is violated
           return 0;
       }
       $og_ancestor = $og_ancestor->ancestor;
   }
   return 1;
}


=head2 is_paraphyletic

 Title   : is_paraphyletic
 Usage   : if( $tree->is_paraphyletic(-nodes =>\@nodes,
                                      -outgroup => $node) ){ }
 Function: Tests whether or not a given set of nodes are paraphyletic
           (representing the full clade) given an outgroup
 Returns : [-1,0,1] , -1 if the group is not monophyletic
                       0 if the group is not paraphyletic
                       1 if the group is paraphyletic
 Args    : -nodes => Array of Bio::Tree::NodeI objects which are in the tree
           -outgroup => a Bio::Tree::NodeI to compare the nodes to


=cut

sub is_paraphyletic{
   my ($self,@args) = @_;
   my ($nodes,$outgroup) = $self->_rearrange([qw(NODES OUTGROUP)],@args);

   if( ! defined $nodes || ! defined $outgroup ) {
       $self->warn("Must suply -nodes and -outgroup parameters to the method is_paraphyletic");
       return;
   }
   if( ref($nodes) !~ /ARRAY/i ) { 
       $self->warn("Must provide a valid array reference for -nodes");
       return;
   }

   # Algorithm
   # Find the lca
   # Find all the nodes beneath the lca
   # Test to see that none are missing from the nodes list
   my %nodehash;
   foreach my $n ( @$nodes ) {
       $nodehash{$n->internal_id} = $n;
   }

   my $clade_root = $self->get_lca(-nodes => $nodes );
   unless( defined $clade_root ) { 
       $self->warn("could not find clade root via lca");
       return;
   }

   my $og_ancestor = $outgroup->ancestor;

   # Is this necessary/correct for paraphyly test?
   while( defined ($og_ancestor ) ) {
       if( $og_ancestor->internal_id == $clade_root->internal_id ) {
           # monophyly is violated, could be paraphyletic
           return -1;
       }
       $og_ancestor = $og_ancestor->ancestor;
   }
   my $tree = Bio::Tree::Tree->new(-root     => $clade_root,
                                   -nodelete => 1);

   foreach my $n ( $tree->get_nodes() ) { 
       next unless $n->is_Leaf();
       # if any leaf node is not in the list
       # then it is part of the clade and so the list
       # must be paraphyletic
       return 1 unless (  $nodehash{$n->internal_id} );
   }
   return 0;
}


=head2 reroot

 Title   : reroot
 Usage   : $tree->reroot($node);
 Function: Reroots a tree making a new node the root
 Returns : 1 on success, 0 on failure
 Args    : Bio::Tree::NodeI that is in the tree, but is not the current root

=cut

sub reroot {
    my ($self,$new_root) = @_;
    unless (defined $new_root && $new_root->isa("Bio::Tree::NodeI")) {
        $self->warn("Must provide a valid Bio::Tree::NodeI when rerooting");
        return 0;
    }

    my $old_root = $self->get_root_node;
    if( $new_root == $old_root ) {
        $self->warn("Node requested for reroot is already the root node!");
        return 0;
    }
    my $anc = $new_root->ancestor;
    unless( $anc ) {
        # this is already the root
        $self->warn("Node requested for reroot is already the root node!");
        return 0;
    }
    my $tmp_node = $new_root->create_node_on_branch(-position=>0,-force=>1);
    # reverse the ancestor & children pointers
    my $former_anc = $tmp_node->ancestor;
    my @path_from_oldroot = ($self->get_lineage_nodes($tmp_node), $tmp_node);
    for (my $i = 0; $i < $#path_from_oldroot; $i++) {
        my $current = $path_from_oldroot[$i];
        my $next = $path_from_oldroot[$i + 1];
        $current->remove_Descendent($next);
        $current->branch_length($next->branch_length);
        $current->bootstrap($next->bootstrap) if defined $next->bootstrap;
        $next->remove_tag('B');
        $next->add_Descendent($current);
    }

    $new_root->add_Descendent($former_anc);
    $tmp_node->remove_Descendent($former_anc);
    
    $tmp_node = undef;
    $new_root->branch_length(undef);

    $old_root = undef;
    $self->set_root_node($new_root);

    return 1;
}


=head2 reroot_at_midpoint

 Title   : reroot_at_midpoint
 Usage   : $tree->reroot_at_midpoint($node, $new_root_id);
 Function: Reroots a tree on a new node created halfway between the 
           argument and its ancestor
 Returns : the new midpoint Bio::Tree::NodeIon success, 0 on failure
 Args    : non-root Bio::Tree::NodeI currently in $tree
           scalar string, id for new node (optional)

=cut

sub reroot_at_midpoint {
    my $self = shift;
    my $node = shift;
    my $id = shift;

    unless (defined $node && $node->isa("Bio::Tree::NodeI")) {
        $self->warn("Must provide a valid Bio::Tree::NodeI when rerooting");
        return 0;
    }

    my $midpt = $node->create_node_on_branch(-FRACTION=>0.5);
    if (defined $id) {
        $self->warn("ID argument is not a scalar") if (ref $id);
        $midpt->id($id) if defined($id) && !ref($id);
    }
    $self->reroot($midpt);
    return $midpt;
}


=head2 findnode_by_id

 Title   : findnode_by_id
 Usage   : my $node = $tree->findnode_by_id($id);
 Function: Get a node by its id (which should be 
           unique for the tree)
 Returns : L<Bio::Tree::NodeI>
 Args    : node id


=cut

sub findnode_by_id {
    my $tree = shift;
    $tree->deprecated("use of findnode_by_id() is deprecated; ".
                      "use find_node() instead");
    my $id = shift;
    my $rootnode = $tree->get_root_node;
    if ( ($rootnode->id) and ($rootnode->id eq $id) ) {
        return $rootnode;
    }
    # process all the children
    foreach my $node ( $rootnode->get_Descendents ) {
        if ( ($node->id) and ($node->id eq $id ) ) {
            return $node;
        }
    }
}


=head2 move_id_to_bootstrap

 Title   : move_id_to_bootstrap
 Usage   : $tree->move_id_to_bootstrap
 Function: Move internal IDs to bootstrap slot
 Returns : undef
 Args    : undef

=cut

sub move_id_to_bootstrap{
   my ($tree) = shift;
   for my $node ( grep { ! $_->is_Leaf } $tree->get_nodes ) {
       $node->bootstrap(defined $node->id ? $node->id : '');
       $node->id('');
   }
}


=head2 add_trait

 Title   : add_trait
 Usage   : my $key = $tree->add_trait($trait_file, 3);
 Function: Add traits to the leaf nodes of a Bio::Tree:Tree from a file.
           The trait file is a tab-delimited text file and needs to have a
           header line giving names to traits. The first column contains the
           leaf node ids. Subsequent columns contain different trait value sets.
           Single or double quotes are removed from the trait values. Traits
           are added to leaf nodes as a tag named $key using the add_tag_value()
           method. This means that you can retrieve the trait values using the
           get_tag_values() method (see the documentation for Bio::Tree::Node).
 Returns : Trait name (a scalar) on success, undef on failure (for example, if
           the column index requested was too large).
 Args    : * Name of trait file (scalar string).
           * Index of trait file column (scalar int). Note that numbering starts
             at 0. Default: 1 (second column).
           * Ignore missing values. Typically, if a leaf node has no value in
             the trait file, an exception is thrown. If you set this option to
             1, then no trait will be given to the node (no exception thrown).

=cut

sub _read_trait_file {
    my ($self, $file, $column) = @_;
    $column ||= 1;

    my $trait_name;
    my $trait_values;
    open my $TRAIT, '<', $file or $self->throw("Could not read file '$file': $!");

    my $first_line = 1;
    while (<$TRAIT>) {
        chomp;
        s/['"]//g;
        my @line = split /\t/;
        if ($first_line) {
            $first_line = 0;
            $trait_name = $line[$column];
            next;
        }

        my $id = $line[0];
        last if (not defined $id) or ($id eq '');

        # Skip empty trait values
        my $value = $line[$column];
        next if (not defined $value) or ($value eq '');

        $trait_values->{$id} = $value;
    }

    close $TRAIT;
    return $trait_name, $trait_values;
}

sub add_trait {
    my ($self, $file, $column, $ignore) = @_;
    $ignore = 0 if not defined $ignore;

    my ($trait_name, $trait_values) = $self->_read_trait_file($file, $column);

    if (defined $trait_name) {

        for my $node ($self->get_leaf_nodes) {

            # strip quotes from the node id
            $node->id($1) if $node->id =~ /^['"]+(.*)['"]+$/;

            if ( not exists $trait_values->{$node->id} ) {
                if ($ignore) {
                    next;
                } else {
                    $self->throw("No trait for node [".$node->id."/".$node->internal_id."]");
                }
            }

            $node->add_tag_value($trait_name, $trait_values->{ $node->id } );

        }
    }
    return $trait_name;
}


1;
