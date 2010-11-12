#
# BioPerl module for Bio::Tree::Node
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::Node - A Simple Tree Node

=head1 SYNOPSIS

    use Bio::Tree::Node;
    my $nodeA = Bio::Tree::Node->new();
    my $nodeL = Bio::Tree::Node->new();
    my $nodeR = Bio::Tree::Node->new();

    my $node = Bio::Tree::Node->new();
    $node->add_Descendent($nodeL);
    $node->add_Descendent($nodeR);

    print "node is not a leaf \n" if( $node->is_leaf);

=head1 DESCRIPTION

Makes a Tree Node suitable for building a Tree.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 CONTRIBUTORS

Aaron Mackey, amackey-at-virginia-dot-edu
Sendu Bala,   bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Tree::Node;
use vars qw($CREATIONORDER);
use strict;

use base qw(Bio::Root::Root Bio::Tree::NodeI Bio::Tree::NodeFunctionsI Bio::Tree::TagValueHolder);

BEGIN {
    $CREATIONORDER = 1;
}

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tree::Node->new();
 Function: Builds a new Bio::Tree::Node object
 Returns : Bio::Tree::Node
 Args    : -descendents   => arrayref of descendents (they will be
                             updated s.t. their ancestor point is this
                             node)
           -branch_length => branch length [integer] (optional)
           -bootstrap     => value   bootstrap value (string)
           -description   => description of node
           -id            => human readable id for node

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($children, $branchlen,$id,
      $bootstrap, $desc,$d) = $self->_rearrange([qw(
						    DESCENDENTS
						    BRANCH_LENGTH
						    ID
						    BOOTSTRAP
						    DESC
						    DESCRIPTION
						 )],
					     @args);
  $self->_register_for_cleanup(\&node_cleanup);
  $self->{'_desc'} = {}; # for descendents
  if( defined $d && defined $desc ) { 
      $self->warn("can only accept -desc or -description, not both, accepting -description");
      $desc = $d;
  } elsif( defined $d && ! defined $desc ) {
      $desc = $d;
  }
  defined $desc && $self->description($desc);
  defined $bootstrap && $self->bootstrap($bootstrap);
  defined $id && $self->id($id);
  defined $branchlen && $self->branch_length($branchlen);
  if( defined $children ) {
      if( ref($children) !~ /ARRAY/i ) {
	  $self->warn("Must specify a valid ARRAY reference to initialize a Node's Descendents");
      }
      foreach my $c ( @$children ) { 	
	  $self->add_Descendent($c);
      }
  }
  $self->_creation_id($CREATIONORDER++);
  return $self;
}

=head2 add_child

 Title   : add_child
 Usage   : $node->add_child($node);
 Function: Adds a child to a node
 Returns : nothing.
 Args    : Bio::Node::NodeI
           boolean flag, true if you want to ignore the fact that you are
           adding a second node with the same unique id (typically memory 
           location reference in this implementation).  default is false and 
           will throw an error if you try and overwrite an existing node.

=cut

sub add_child{
   my ($self,$node,$ignoreoverwrite) = @_;
   return -1 if( ! defined $node );
   
   if( ! ref($node) ||
       ref($node) =~ /HASH/ ||
       ! $node->isa('Bio::Tree::NodeI') ) {
       $self->throw("Trying to add a Descendent who is not a Bio::Tree::NodeI");
       return -1;
   }
   
   $self->{_adding_descendent} = 1;
   # avoid infinite recurse
   $node->ancestor($self) unless $node->{_setting_ancestor}; 
   $self->{_adding_descendent} = 0;
   
   if( $self->{'_desc'}->{$node->internal_id} && ! $ignoreoverwrite ) {
       $self->throw("Going to overwrite a node which is $node that is already stored here, set the ignore overwrite flag (parameter 2) to true to ignore this in the future");
   }
   $self->{'_desc'}->{$node->internal_id} = $node; # is this safely unique - we've tested before at any rate??
   
}

=head2 remove_child

 Title   : remove_child
 Usage   : $node->remove_child($node_foo);
 Function: Removes a specific node from being a child of this node
 Returns : nothing
 Args    : An array of Bio::Node::NodeI objects

=cut

sub remove_child{
   my ($self,@nodes) = @_;
   my $c= 0;
   foreach my $n ( @nodes ) { 
       if( $self->{'_desc'}->{$n->internal_id} ) {
        $self->{_removing_descendent} = 1;
        $n->ancestor(undef);
        $self->{_removing_descendent} = 0;
	   # should be redundant
	   $self->{'_desc'}->{$n->internal_id}->ancestor(undef);
	   delete $self->{'_desc'}->{$n->internal_id};
	   $c++;
       } else { 
	   if( $self->verbose ) {
	       $self->debug(sprintf("no node %s (%s) listed as a descendent in this node %s (%s)\n",$n->id, $n,$self->id,$self));
	       $self->debug("Descendents are " . join(',', keys %{$self->{'_desc'}})."\n");
	   }
       }
   }
   $c;
}


=head2 children

 Title   : children
 Usage   : my @nodes = $node->children;
 Function: all the descendents for this Node (but not their descendents
					      i.e. not a recursive fetchall)
 Returns : Array of Bio::Tree::NodeI objects
 Args    : none.

=cut

sub children {
   my ($self) = @_;

   return map { $_->[0] }
   sort { $a->[1] <=> $b->[1] } 
   map { [$_, $_->internal_id ] }
   grep {defined $_}
   values %{$self->{'_desc'}};

   # Ignore everything below -- we don't do sorting anymore!
   
   # order can be based on branch length (and sub branchlength)   
   my $sortby;
   $sortby ||= 'none';
   if (ref $sortby eq 'CODE') {
       my @values = sort { $sortby->($a,$b) } values %{$self->{'_desc'}};
       return @values;
   } elsif ($sortby eq 'height') {
       return map { $_->[0] }
       sort { $a->[1] <=> $b->[1] || 
		  $a->[2] <=> $b->[2] } 
       map { [$_, $_->height, $_->internal_id ] } 
       values %{$self->{'_desc'}};
   } elsif( $sortby eq 'alpha' ) {
       my @set;
       for my $v ( values %{$self->{'_desc'}} ) {
	   unless( $v->is_Leaf ) {
	       my @lst = ( sort { $a cmp $b } map { $_->id } 
                          grep { $_->is_Leaf } 
			   $v->get_all_Descendents($sortby));
	       push @set, [$v, $lst[0], $v->internal_id];
	   } else {
	       push @set, [$v, $v->id, $v->internal_id];
	   }
       } 
       return map { $_->[0] }
       sort {$a->[1] cmp $b->[1] || $a->[2] <=> $b->[2] } @set;       
   } elsif( $sortby eq 'revalpha' ) {
       my @set;
       for my $v ( values %{$self->{'_desc'}} ) {
	   if( ! defined $v->id && 
	       ! $v->is_Leaf ) {
	       my ($l) = ( sort { $b cmp $a } map { $_->id } 
			   grep { $_->is_Leaf } 
			   $v->get_all_Descendents($sortby));
	       push @set, [$v, $l, $v->internal_id];
	   } else { 
	       push @set, [$v, $v->id, $v->internal_id];
	   }
       } 
       return map { $_->[0] }
       sort {$b->[1] cmp $a->[1] || $b->[2] <=> $a->[2] } @set;
   } else { # creation
       return map { $_->[0] }
       sort { $a->[1] <=> $b->[1] } 
       map { [$_, $_->internal_id ] }
       grep {defined $_}
       values %{$self->{'_desc'}};	   
   }
}


=head2 set_child_order

 Title   : set_child_order
 Usage   : $node->set_child_order($child_a,$child_b,$child_c);
 Function: Explicitly sets the sort order of the children beneath this node
 Returns : nothing
 Args    : List of Bio::Tree::Node objects, which *must* contain *all* of the direct
           descendents of the current node

=cut
sub set_child_order {
    my $self = shift;
    my @nodes_in_new_order = @_;

    # Collect the sorted internal_ids of all the children nodes.
    my @internal_ids = map {$_->_creation_id} $self->children;
    @internal_ids = sort {$a <=> $b}  @internal_ids;

    # Explicitly re-apply the sorted creation IDs to the new order of children
    foreach my $child (@nodes_in_new_order) {
	$child->_creation_id(shift @internal_ids);
    }
}


=head2 is_leaf

 Title   : is_leaf
 Usage   : if( $node->is_leaf )
 Function: Get leaf status
 Returns : boolean
 Args    : none

=cut

sub is_leaf {
    my ($self) = @_;
    my $isleaf = ! (defined $self->{'_desc'} &&
		 (keys %{$self->{'_desc'}} > 0) );
    return $isleaf;
}

=head2 branch_length

 Title   : branch_length
 Usage   : $obj->branch_length()
 Function: Get/Set the branch length
 Returns : value of branch_length
 Args    : newvalue (optional)

=cut

sub branch_length{
    my $self = shift;
    if( @_ ) {
	my $bl = shift;
	if( defined $bl &&
	    $bl =~ s/\[(\d+)\]// ) {
	    $self->bootstrap($1);
	}
	$self->{'_branch_length'} = $bl;
    }
    my $current_bl = $self->{'_branch_length'};
    #$current_bl = 1 unless (defined $current_bl && $current_bl =~ /^\-?\d+(\.\d+)?$/);
    return $current_bl;
}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: The human readable identifier for the node 
 Returns : value of human readable id
 Args    : newvalue (optional)

"A name can be any string of printable characters except blanks,
colons, semicolons, parentheses, and square brackets. Because you may
want to include a blank in a name, it is assumed that an underscore
character ("_") stands for a blank; any of these in a name will be
converted to a blank when it is read in."  

from L<http://evolution.genetics.washington.edu/phylip/newicktree.html>

Also note that these objects now support spaces, ();: because we can
automatically quote the strings if they contain these characters.  The
L<id_output> method does this for you so use the id() method to get
the raw string while L<id_output> to get the pre-escaped string.

=cut

sub id {
    my ($self, $value) = @_;
    if (defined $value) {
        #$self->warn("Illegal characters ();:  and space in the id [$value], converting to _ ")
	# if $value =~ /\(\);:/ and $self->verbose >= 0;
        #$value =~ s/[\(\);:\s]/_/g;
        $self->{'_id'} = $value;
    }
    my $val = $self->{'_id'};
    $val = '' if (!defined $val);
    return $val;
}

sub internal_id {
   return $_[0]->_creation_id;
}

sub bootstrap { 
    my $self = shift;
    if( @_ ) {
	if( $self->has_tag('B') ) {
	    $self->remove_tag('B');
	}
	$self->add_tag_value('B',shift);
    }
    return ($self->get_tag_values('B'))[0];
}

sub description { 
    my $self = shift;
    if( @_ ) {
	if( $self->has_tag('description') ) {
	    $self->remove_tag('description');
	}
	$self->add_tag_value('description',shift);
    }
    return ($self->get_tag_values('description'))[0];
}


sub parent {
    my $self = shift;
    if (@_) {
        my $new_ancestor = shift;
        
        # we can set ancestor to undef
        if ($new_ancestor) {
            $self->throw("This is [$new_ancestor], not a Bio::Tree::NodeI")
		unless $new_ancestor->isa('Bio::Tree::NodeI');
        }
        
        my $old_ancestor = $self->{'_ancestor'} || '';
        if (!$old_ancestor || 
	    ($old_ancestor && ( !$new_ancestor || 
			       $new_ancestor ne $old_ancestor)) ) {
            if( $old_ancestor && ! $old_ancestor->{_removing_descendent}) {
		$old_ancestor->remove_Descendent($self);
	    }
            if ($new_ancestor && 
		! $new_ancestor->{_adding_descendent} ) { # avoid infinite recurse
                $self->{_setting_ancestor} = 1;
                $new_ancestor->add_Descendent($self, 1);
                $self->{_setting_ancestor} = 0;
            }
        }
        $self->{'_ancestor'} = $new_ancestor;
    }
    
    return $self->{'_ancestor'};
}


sub _creation_id {
    my $self = shift @_;
    $self->{'_creation_id'} = shift @_ if( @_);
    return $self->{'_creation_id'} || 0;
}

sub node_cleanup {
    my $self = shift;
    return unless defined $self;
    
    #*** below is wrong, cleanup doesn't actually occur. Will replace with:
    # $self->remove_all_Descendents; once further fixes in place..
    #if( defined $self->{'_desc'} &&
    #    ref($self->{'_desc'}) =~ /HASH/i ) {
    #    while( my ($nodeid,$node) = each %{ $self->{'_desc'} } ) {
    #        $node->ancestor(undef); # insure no circular references
    #        $node = undef;
    #    }
    #}
    $self->remove_children;
    
    #$self->{'_desc'} = {};
    1;
}

sub find_by_tag_value {
    my $self = shift;
    my $tag = shift;
    my $value = shift;

    my @nodes;
    foreach my $node ($self->nodes) {
	next unless ($node->can('get_tag_values'));
	my @values = $node->get_tag_values($tag);
	if (grep {$_ =~ $value} @values) {
	    push @nodes, $node ;
	}
    }

    if ( wantarray) { 
	return @nodes;
    } else { 
	if( @nodes > 1 ) { 
	    $self->warn("More than 1 node found but caller requested scalar, only returning first node");
       }
	return shift @nodes;
    }
}


1;
