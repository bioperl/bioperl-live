#
# BioPerl module for Bio::Taxonomy::Tree
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Dan Kortschak but pilfered extensively from Bio::Tree::Tree by Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Taxonomy::Tree - An Organism Level Implementation of TreeI interface.

=head1 SYNOPSIS

    # NB: This module is deprecated. Use Bio::Taxon in combination with
    # Bio::Tree::Tree instead

    # like from a TreeIO
    my $treeio = Bio::TreeIO->new(-format => 'newick', -file => 'treefile.dnd');
    my $tree = $treeio->next_tree;
    my @nodes = $tree->get_nodes;
    my $root = $tree->get_root_node;
    my @leaves = $tree->get_leaves;


=head1 DESCRIPTION

This object holds handles to Taxonomic Nodes which make up a tree.

=head1 EXAMPLES

  use Bio::Species;
  use Bio::Taxonomy::Tree;

  my $human=Bio::Species->new();
  my $chimp=Bio::Species->new();
  my $bonobo=Bio::Species->new();

  $human->classification(qw( sapiens Homo Hominidae
                             Catarrhini Primates Eutheria
                             Mammalia Euteleostomi Vertebrata 
                             Craniata Chordata
                             Metazoa Eukaryota ));
  $chimp->classification(qw( troglodytes Pan Hominidae
                             Catarrhini Primates Eutheria
                             Mammalia Euteleostomi Vertebrata 
                             Craniata Chordata
                             Metazoa Eukaryota ));
  $bonobo->classification(qw( paniscus Pan Hominidae
                              Catarrhini Primates Eutheria
                              Mammalia Euteleostomi Vertebrata 
                              Craniata Chordata
                              Metazoa Eukaryota ));

  # ranks passed to $taxonomy match ranks of species
  my @ranks = ('superkingdom','kingdom','phylum','subphylum',
               'no rank 1','no rank 2','class','no rank 3','order',
               'suborder','family','genus','species');

  my $taxonomy=Bio::Taxonomy->new(-ranks => \@ranks,
                                 -method => 'trust',
                                 -order => -1);


  my $tree1=Bio::Taxonomy::Tree->new();
  my $tree2=Bio::Taxonomy::Tree->new();

  $tree1->make_species_branch($human,$taxonomy);
  $tree2->make_species_branch($chimp,$taxonomy);

  my ($homo_sapiens)=$tree1->get_leaves;

  $tree1->splice($tree2);

  $tree1->add_species($bonobo,$taxonomy);

  my @taxa;
  foreach my $leaf ($tree1->get_leaves) {
     push @taxa,$leaf->taxon;
  }
  print join(", ",@taxa)."\n";

  @taxa=();
  $tree1->remove_branch($homo_sapiens);
  foreach my $leaf ($tree1->get_leaves) {
     push @taxa,$leaf->taxon;
  }
  print join(", ",@taxa)."\n";

=head1 FEEDBACK

See AUTHOR

=head1 AUTHOR - Dan Kortschak

Email kortschak@rsbs.anu.edu.au

=head1 CONTRIBUTORS

Mainly Jason Stajich

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Code begins...


package Bio::Taxonomy::Tree;
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Taxonomy::Taxon;

# Import rank information from Bio::Taxonomy.pm
use vars qw(@RANK %RANK);

use base qw(Bio::Root::Root Bio::Tree::TreeI Bio::Tree::TreeFunctionsI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Taxonomy::Tree->new();
 Function: Builds a new Bio::Taxonomy::Tree object 
 Returns : Bio::Taxonomy::Tree
 Args    : 


=cut

sub new {
  my($class,@args) = @_;
  
  my $self = $class->SUPER::new(@args);
  $self->warn("Bio::Taxonomy::Tree is deprecated. Use Bio::Taxon in combination with Bio::Tree::Tree instead.");
  
  $self->{'_rootnode'} = undef;
  $self->{'_maxbranchlen'} = 0;

  my ($root)= $self->_rearrange([qw(ROOT)], @args);
  if( $root ) { $self->set_root_node($root); }
  return $self;
}


=head2 get_nodes

 Title   : get_nodes
 Usage   : my @nodes = $tree->get_nodes()
 Function: Return list of Bio::Taxonomy::Taxon objects
 Returns : array of Bio::Taxonomy::Taxon objects
 Args    : (named values) hash with one value 
           order => 'b|breadth' first order or 'd|depth' first order

=cut

sub get_nodes{
   my ($self, @args) = @_;
   
   my ($order, $sortby) = $self->_rearrange([qw(ORDER SORTBY)],@args);
   $order ||= 'depth';
   $sortby ||= 'height';

   if ($order =~ m/^b|(breadth)$/oi) {
      my $node = $self->get_root_node;
      my @children = ($node);
      for (@children) {
	   push @children, $_->each_Descendent($sortby);
      }
      return @children;
   }

   if ($order =~ m/^d|(depth)$/oi) {
       # this is depth-first search I believe
       my $node = $self->get_root_node;
       my @children = ($node,$node->get_Descendents($sortby));
       return @children;
   }
}

=head2 get_root_node

 Title   : get_root_node
 Usage   : my $node = $tree->get_root_node();
 Function: Get the Top Node in the tree, in this implementation
           Trees only have one top node.
 Returns : Bio::Taxonomy::Taxon object
 Args    : none

=cut


sub get_root_node{
   my ($self) = @_;
   return $self->{'_rootnode'};
}

=head2 set_root_node

 Title   : set_root_node
 Usage   : $tree->set_root_node($node)
 Function: Set the Root Node for the Tree
 Returns : Bio::Taxonomy::Taxon
 Args    : Bio::Taxonomy::Taxon

=cut


sub set_root_node{
   my ($self,$value) = @_;
   if( defined $value ) { 
      if( ! $value->isa('Bio::Taxonomy::Taxon') ) { 
	   $self->warn("Trying to set the root node to $value which is not a Bio::Taxonomy::Taxon");
	   return $self->get_root_node;
      }
      $self->{'_rootnode'} = $value;
   }
   return $self->get_root_node;
}


=head2 get_leaves

 Title   : get_leaves
 Usage   : my @nodes = $tree->get_leaves()
 Function: Return list of Bio::Taxonomy::Taxon objects
 Returns : array of Bio::Taxonomy::Taxon objects
 Args    : 

=cut


sub get_leaves{
   my ($self) = @_;
   
   my $node = $self->get_root_node;
   my @leaves;
   my @children = ($node);
   for (@children) {
      push @children, $_->each_Descendent();
   }
   for (@children) {
      push @leaves, $_ if $_->is_Leaf;
   }
   return @leaves;
}

=head2 make_species_branch

 Title   : make_species_branch
 Usage   : @nodes = $tree->make_species_branch($species,$taxonomy)
 Function: Return list of Bio::Taxonomy::Taxon objects based on a Bio::Species object
 Returns : array of Bio::Taxonomy::Taxon objects
 Args    : Bio::Species and Bio::Taxonomy objects

=cut

# I'm not happy that make_species_branch and make_branch are seperate routines
# should be able to just make_branch and have it sort things out

sub make_species_branch{
   my ($self,$species,$taxonomy) = @_;
   
   if (! $species->isa('Bio::Species') ) {
      $self->throw("Trying to classify $species which is not a Bio::Species object");
   }
   if (! $taxonomy->isa('Bio::Taxonomy') ) {
      $self->throw("Trying to classify with $taxonomy which is not a Bio::Taxonomy object");
   }

   # this is done to make sure we aren't duplicating a path (let God sort them out)
   if (defined $self->get_root_node) {
      $self->get_root_node->remove_all_Descendents;
   }

   my @nodes;

   # nb taxa in [i][0] and ranks in [i][1]
   my @taxa=$taxonomy->classify($species);

   for (my $i = 0; $i < @taxa; $i++) {
      $nodes[$i]=Bio::Taxonomy::Taxon->new(-taxon => $taxa[$i][0],
                                           -rank  => $taxa[$i][1]);
   }

   for (my $i = 0; $i < @taxa-1; $i++) {
      $nodes[$i]->add_Descendent($nodes[$i+1]);
   }

   $self->set_root_node($nodes[0]);

   return @nodes;
}


=head2 make_branch

 Title   : make_branch
 Usage   : $tree->make_branch($node)
 Function: Make a linear Bio::Taxonomy::Tree object from a leafish node
 Returns :
 Args    : Bio::Taxonomy::Taxon object

=cut


sub make_branch{
   my ($self,$node) = @_;

   # this is done to make sure we aren't duplicating a path (let God sort them out)
   # note that if you are using a linked set of node which include node 
   # already in the tree, this will break
   $self->get_root_node->remove_all_Descendents;
   
   while (defined $node->ancestor) {
      $self->set_root_node($node);
      $node=$node->ancestor;
   }
}


=head2 splice

 Title   : splice
 Usage   : @nodes = $tree->splice($tree)
 Function: Return a of Bio::Taxonomy::Tree object that is a fusion of two
 Returns : array of Bio::Taxonomy::Taxon added to tree
 Args    : Bio::Taxonomy::Tree object

=cut


sub splice{
   my ($self,$tree) = @_;

   my @nodes;

   my @newleaves = $tree->get_leaves;
   foreach my $leaf (@newleaves) {
      push @nodes,$self->add_branch($leaf);
   }

   return @nodes;
}

=head2 add_species

 Title   : add_species
 Usage   : @nodes = $tree->add_species($species,$taxonomy)
 Function: Return a of Bio::Taxonomy::Tree object with a new species added
 Returns : array of Bio::Taxonomy::Taxon added to tree
 Args    : Bio::Species object

=cut


sub add_species{
   my ($self,$species,$taxonomy) = @_;

   my $branch=Bio::Taxonomy::Tree->new;
   my @nodes=$branch->make_species_branch($species,$taxonomy);

   my ($newleaf)=$branch->get_leaves;
  
   return $self->add_branch($newleaf);
}

=head2 add_branch

 Title   : add_branch
 Usage   : $tree->add_branch($node,boolean)
 Function: Return a of Bio::Taxonomy::Tree object with a new branch added
 Returns : array of Bio::Taxonomy::Taxon objects of the resulting tree
 Args    : Bio::Taxonomy::Taxon object
           boolean flag to force overwrite of descendent
             (see Bio::Node->add_Descendent)

=cut


sub add_branch {
   my ($self,$node,$force) = @_;

   my $best_node_level=0;
   my ($best_node,@nodes,$common);

   my @leaves=$self->get_leaves;
   foreach my $leaf (@leaves) {
      $common=$node->recent_common_ancestor($leaf); # the root of the part to add
      if (defined $common && ($common->distance_to_root > $best_node_level)) {
         $best_node_level = $common->distance_to_root;
         $best_node = $common;
      }
   }

   return unless defined $best_node;

   push @nodes,($self->get_root_node,$self->get_root_node->get_Descendents);
   foreach my $node (@nodes) {
      if ((defined $best_node->id && $best_node->id == $node->id) ||
         ($best_node->rank eq $node->rank && $best_node->taxon eq $node->taxon) &&
         ($best_node->rank ne 'no rank')) {
         foreach my $descendent ($common->each_Descendent) {
            $node->add_Descendent($descendent,$force);
         }
      }

      $self->set_root_node($node) if $node->distance_to_root==0;
   }

   return ($common->get_Descendents);
}

=head2 remove_branch

 Title   : remove_branch
 Usage   : $tree->remove_branch($node)
 Function: remove a branch up to the next multifurcation
 Returns :
 Args    : Bio::Taxonomy::Taxon object

=cut


sub remove_branch{
   my ($self,$node) = @_;

   # we can define a branch at any point along it
   
   while (defined $node->ancestor) {
      last if $node->ancestor->each_Descendent > 1;
      $node=$node->ancestor;
   }
   $node->remove_all_Descendents; # I'm not sure if this is necessary,
                                  # but I don't see that remove_Descendent
                                  # has the side effect of deleting
                                  # descendent nodes of the deletee
   $node->ancestor->remove_Descendent($node);
}



1;
