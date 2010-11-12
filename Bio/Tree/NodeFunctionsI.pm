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

Bio::Tree::NodeFunctionsI - Decorated Interface implementing basic Node methods

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

  http://bugzilla.open-bio.org/

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

package Bio::Tree::NodeFunctionsI;
use strict;

use base qw(Bio::Tree::NodeI);


=head2 id_output

 Title   : id_output
 Usage   : my $id = $node->id_output;
 Function: Return an id suitable for output in format like newick
           so that if it contains spaces or ():; characters it is properly 
           quoted
 Returns : $id string if $node->id has a value
 Args    : none


=cut

sub id_output{
    my $node = shift;
    my $id = $node->id;
    return unless( defined $id && length($id ) );
    # single quotes must become double quotes
    # $id =~ s/'/''/g;
    if( $id =~ /[\(\);:,\s]/ ) {
	$id = '"'.$id.'"';
    }
    return $id;
}

sub depth_to_root {
   my ($self) = @_;
   
   my $depth = 1;
   my $node = $self;
   while( defined $node->ancestor ) { 
       $depth += 1;
       $node = $node->ancestor;
   }
   return $depth;
}

sub height_to_root {
    my ($self) = @_;

   my $height = 0;
   my $node = $self;
   while( defined $node->ancestor ) { 
       $height += $node->branch_length;
       $node = $node->ancestor;
   }
   return $height;
}
sub distance_to_root { shift->height_to_root(@_) }


sub max_distance_to_leaf{
    my ($self) = @_;

    return 0 if( $self->is_leaf );
    
    my $max = 0;
    foreach my $subnode ( $self->children ) { 
	# Collect a 'safe' version of the branch length which returns a 'default' 
	# 1 if the branch length is not defined.
	my $s = $subnode->max_distance_to_leaf + $subnode->branch_length_or_one;
	if( $s > $max ) { $max = $s; }
    }
    return $max;
}
sub max_height_to_leaf { shift->max_distance_to_leaf(@_) }
sub height { shift->max_distance_to_leaf(@_) }


sub max_depth_to_leaf{
    my ($self) = @_;

    return 1 if( $self->is_leaf );
    
    my $max = 0;
    foreach my $subnode ( $self->children ) { 
	my $s = $subnode->max_depth_to_leaf;
	if( $s > $max ) { $max = $s; }
    }
    return $max;
}
sub depth { shift->max_depth_to_leaf(@_) }


sub total_branch_length {
    my $self = shift;
    my $sum = 0;
    for ( $self->nodes ) {
    $sum += $_->branch_length || 0;
    }
    return $sum;
}
sub subtree_length { shift->total_branch_length(@_) }
sub subtree_size { shift->total_branch_length(@_) }


sub leaves {
    my $self = shift;

    my @all_nodes = $self->nodes();
    return grep {$_->is_leaf} @all_nodes;
}
sub get_all_leaves { shift->leaves(@_) }


sub nodes {
   my ($self) = @_;
   my @nodes = ();
   foreach my $node ( $self->children() ) {
       push @nodes, ($node,$node->nodes());
   }
   return @nodes;
}
sub get_all_Descendents { shift->nodes(@_) }
sub get_Descendents { 
    my $self = shift;
    $self->deprecated(-message     => 'get_Descendents() is deprecated... use nodes() instead');
    $self->nodes(@_);
}
sub get_all_nodes { shift->nodes(@_) }


sub _ordered_nodes {
    my $self = shift;
    my $order = shift;
   
    $order ||= 'depth';

    if ($order =~ m/^b|(breadth)$/oi) {
       my @children = ($self);
       for (@children) {
        push @children, $_->children;
       }
       return @children;
   }

   if ($order =~ m/^d|(depth)$/oi) {
       # this is depth-first search I believe
       my @children = ($self,$self->nodes);
       return @children;
   }
}

sub nodes_breadth_first {
    my $self = shift;
    return $self->_ordered_nodes('b');
}

sub nodes_depth_first {
    my $self = shift;
    return $self->_ordered_nodes('d');
}


sub lineage {
    my $self = shift;

    my @lineage = ();
    my $node = $self->parent;
    while ($node) {
        unshift(@lineage, $node);
        $node = $node->parent;
    }
    return @lineage;
}


sub change_child_to_parent {
    my $self = shift;
    my $child_to_become_parent = shift;

    if ($self->parent) {
	# Climb up the tree, flipping the direction of all nodes up to the root.
	$self->parent->change_child_to_parent($self);
    }

    # Store the branch length between us and our child.
    my $bl = $child_to_become_parent->branch_length;

    # Remove $child from me, add me to $child
    $self->remove_child($child_to_become_parent);
    $child_to_become_parent->add_child($self);
    # Set my branch length appropriately
    $self->branch_length($bl);
}

sub child_count {
    my $self = shift;
    return scalar($self->children);
}

sub node_count {
   my ($self) = @_;
   return scalar($self->nodes);
}
sub descendent_count { shift->node_count(@_) }
sub number_nodes { shift->node_count(@_) }

sub leaf_count {
   my ($self) = @_;
   return scalar($self->leaves);
}

sub reverse_children {
    my $self = shift;

    my @children = $self->children;
    @children = reverse @children;
    $self->set_child_order(@children);
}

sub flip_subtree {
    my $self = shift;

    foreach my $node ($self->nodes) {
	$node->reverse_children;
    }
    $self->reverse_children;
}

sub split_branch_with_node {
    # Splits the branch above this node with the given node.
    my $self = shift;
    my $node_to_insert = shift;
    my $fraction_above_node = shift;

    $fraction_above_node = 0.5 unless (defined $fraction_above_node);

    if (!defined $self->parent) {
	$self->warn("Cannot split branch above the root node -- Not inserting node!");
	return undef;
    }

    my $parent = $self->parent;
    my $bl = $self->branch_length;

    $parent->add_child($node_to_insert);
    $parent->remove_child($self);
    $node_to_insert->add_child($self);

    $node_to_insert->branch_length($bl * (1-$fraction_above_node));
    $self->branch_length($bl * $fraction_above_node);
    
    return $node_to_insert;
}


sub remove_children{
   my ($self) = @_;

   foreach my $node ($self->children) {
       $self->remove_child($node);
   }   
}
sub remove_Descendents { shift->remove_children(@_) }

sub remove {
    my $self = shift;
    if ($self->parent) {
	$self->parent->remove_child($self);
    }
}

sub splice {
    my $self = shift;

    if ($self->parent) {
	my $parent = $self->parent;
	my $bl = $self->branch_length || 0;
	foreach my $child ($self->children) {
	    # Remove the child from me, add it to my parent.
	    $self->remove_child($child);
	    $parent->add_child($child);

	    # Add my branch length to my child -- hence this method is called "splice"
	    my $child_bl = $child->branch_length || 0;
	    $child->branch_length($child_bl + $bl);
	}
	$parent->remove_child($self);
    }
}

sub force_binary {
    my $self = shift;

    foreach my $child ($self->children) {
	# Recurse through tree.
	$child->force_binary;
    }

    my @children = $self->children;
    while (scalar(@children) > 2) {
	my $first_child = shift @children;
	
	# Make a new node to hold a binary split.
	my $new_node = new $self;
	$first_child->parent->add_child($new_node);
	foreach my $cur_child (@children) {
	    # Remove the "extra" children from the direct parent
	    $cur_child->remove;
	    # and add them to $new_node
	    $new_node->add_child($cur_child);
	}
    }
}

sub is_binary {
    my $self = shift;
    return ($self->is_leaf || scalar($self->children) == 2);
}

sub is_subtree_binary {
    my $self = shift;

    foreach my $node ($self->nodes) {
	return 0 unless ($node->is_binary);
    }
    return 1;
}


sub contract_linear_paths {
    my $self = shift;

    foreach my $child ($self->children) {
	$child->contract_linear_paths;
    }
    
    if ($self->parent && $self->child_count == 1) {
	$self->splice;
    }
}
sub remove_elbow_nodes { shift->contract_linear_paths(@_) }


sub branch_length_or_one {
    my $self = shift;
    
    my $bl = $self->branch_length;
    return 1 unless (defined $bl);
    return $bl;
}


=head2 ascii

 Title   : ascii
 Usage   : print STDERR $tree->ascii."\n";
 Function: Return a representation of this tree as ASCII text. Shamelessly
           copied / "ported" from PyCogent. http://pycogent.sourceforge.net/
 Returns : a multi-line String, suitable for printing to the console
 Args    : show_internal - 0 to hide internal nodes, 1 to show them (default 1)
           compact - 1 to use a 'compact' mode with exactly 1 line per node. (default 0)
           ignore_branchlengths - 0 to ignore branch lengths, 1 to use them. (default 0)
=cut

sub ascii {
    my $self = shift;
    my $show_internal = shift;
    my $compact = shift;
    my $ignore_branchlengths = shift;

    my $max_bl = $self->max_distance_to_leaf;
    if ($max_bl == 0) {
	$max_bl = $self->max_depth_to_leaf;
    }
    my $width = 80;

    my $char_per_bl = int($width / $max_bl);

    my ($lines_arrayref,$mid) = $self->_ascii_art($self,'-',$char_per_bl,$show_internal,$compact,$ignore_branchlengths);
    my @lines = @{$lines_arrayref};
    return join("\n",@lines)."\n";
}

# Private helper method for $tree->ascii_art
sub _ascii_art {
    my $self = shift;
    my $node = shift;
    my $char1 = shift;
    my $char_per_bl = shift;
    my $show_internal = shift;
    my $compact = shift;
    my $ignore_branchlengths = shift;

    $char1 = '-' unless (defined $char1);
    $show_internal = 1 unless (defined $show_internal);
    $compact = 0 unless (defined $compact);
    $ignore_branchlengths = 1 unless (defined $ignore_branchlengths);

    my $len = 10;
    if (!$ignore_branchlengths) {
	my $bl = $node->branch_length;
	$bl = 1 unless (defined $bl && $bl =~ /^\-?\d+(\.\d+)?$/);
	$len = int($bl * $char_per_bl);
	$len = 1 if ($len < 1);
    }
    my $pad = ' ' x $len;
    my $pa = ' ' x ($len-1);
    my $name_str = $node->id;
    $name_str = '' unless (defined $name_str);

    if (!$node->is_leaf) {
	my @mids;
	my @results;
	my @children = $node->children;
	for (my $i=0; $i < scalar(@children); $i++) { 
	    my $child = $children[$i];
	    my $char2;
	    if ($i == 0) {
		# First child.
		$char2 = '/';
	    } elsif ($i == scalar(@children) -1) {
		# Last child.
		$char2 = '\\';
	    } else {
		# Middle child.
		$char2 = '-';
	    }
	    my ($clines_arrayref, $mid) = $self->_ascii_art($child,$char2,$char_per_bl,$show_internal,$compact,$ignore_branchlengths);
	    push @mids,($mid+scalar(@results));
	    my @clines = @$clines_arrayref;
	    foreach my $line (@clines) {
		push @results, $line;
	    }
	    if (!$compact) {
		push @results,'';
	    }
	}
	if (!$compact) {
	    pop @results;
	}
	my $lo = $mids[0];
	my $hi = $mids[scalar(@mids)-1];
	my $end = scalar(@results);

	my @prefixes = ();
	push @prefixes, ($pad) x ($lo+1);
	push @prefixes, ($pa.'|') x ($hi-$lo-1);
	push @prefixes, ($pad) x ($end-$hi);
	
	my $mid = int(($lo + $hi) / 2);
	$prefixes[$mid] = $char1 . '-'x($len-2) . substr($prefixes[$mid],length($prefixes[$mid])-1,1);
	my @new_results;
	for (my $i=0; $i < scalar(@prefixes); $i++) {
	    my $p = $prefixes[$i] || '';
	    my $r = $results[$i] || '';
	    push @new_results, ($p . $r);
	}
	@results = @new_results;
	if ($show_internal) {
	    my $stem = $results[$mid];
	    my $str = substr($stem,0,1);
	    $str .= $name_str;
	    $str .= substr($stem,length($name_str)+1);
	    $results[$mid] = $str;
	}
	return (\@results,$mid);
    } else {
	my @results = ($char1 . ('-'x$len) . $name_str);
	if ($ignore_branchlengths) {
	    @results = ($char1 . '-' . $name_str);
	}
	return (\@results,0);
    }
}

sub find_by_id {
    my $self = shift;
    my $value = shift;

    my @nodes = grep { $_->id eq $value } $self->nodes;
    
    if ( wantarray) { 
	return @nodes;
    } else { 
	if( @nodes > 1 ) { 
	    $self->warn("More than 1 node found but caller requested scalar, only returning first node");
       }
	return shift @nodes;
    }
}
sub find { shift->find_by_id(@_) }




1;
