# $Id$
#
# BioPerl module for Bio::Tree::TreeFunctionsI
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::TreeFunctionsI - Decorated Interface implementing basic Tree exploration methods

=head1 SYNOPSIS

use Bio::TreeIO;
my $in = new Bio::TreeIO(-format => 'newick', -file => 'tree.tre');

my $tree = $in->next_tree;

my @nodes = $tree->find_nodes('id1');

if( $tree->is_monophyletic(-clade => @nodes, -outgroup => $outnode) ){

}

=head1 DESCRIPTION

Describe the interface here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich, Aaron Mackey

Email jason@bioperl.org
Email amackey@virginia.edu

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tree::TreeFunctionsI;
use vars qw(@ISA);
use strict;
use  Bio::Tree::TreeI;

@ISA = qw(Bio::Tree::TreeI);


=head2 is_monophyletic

 Title   : is_monophyletic
 Usage   : if( $tree->is_monophyletic(-nodes => \@nodes, 
				      -outgroup => $outgroup)
 Function: Will do a test of monophyly for the nodes specified
           in comparison to a chosen outgroup
 Returns : boolean
 Args    : -nodes => arrayref of nodes to test
           -outgroup => outgroup to serve as a reference


=cut

sub is_monophyletic{
   my ($self,@args) = @_;
   my ($nodes,$outgroup) = $self->_rearrange([qw(NODES OUTGROUP)],@args);
   
   if( ! defined $nodes || ! defined $outgroup ) {
       $self->warn("Must suply -nodes and -outgroup parameters to the method is_monophyletic");
       return undef;
   }
   if( ref($nodes) !~ /ARRAY/i ) { 
       $self->warn("Must provide a valid array reference for -nodes");
   }
   
   # algorithm:
   # For each node, walk up to its ancestor, whenever
   # nodes have the name ancestor, collapse the ancestor
   # list, when the ancestors collapse to 1.
   #
   # See if this collapsed node is an ancestor of 
   # outgroup, if it is that violates monophyly.
   #
   # We use the internal_id method here to be sure
   # we are talking about the same object.
   
   # first just start by filling the hash
   my (%ancestors);
   for my $n ( @$nodes ) {
       $ancestors{$n->internal_id} = $n;
   }  
       
   while( (keys %ancestors) > 1 ) {
       foreach my $n ( values %ancestors ) {
	   my $a = $n->ancestor;
	   if( $a ) {
	       $ancestors{$a->internal_id} = $a;
	       my $id = $n->internal_id;
	       delete $ancestors{$id};
	   }
       } 
   }
   
   my ($clade_root) = values %ancestors;
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
   my ($self,$type,$field) = @_;
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
       $field= $type; 
       $type = 'id';
   } else {   
       $type =~ s/^-//;
   }

   # could actually do this by testing $rootnode->can($type) but
   # it is possible that a tree is implemeted with different node types
   # - although it is unlikely that the root node would be richer than the 
   # leaf nodes.  Can't handle NHX tags right now

   unless( $type eq 'id' || $type eq 'name' ||
	   $type eq 'bootstrap' || $type eq 'description' ||
	   $type eq 'internal_id') {
       $self->warn("unknown search type $type - will try anyways");
   } 
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

1;
