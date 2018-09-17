#
# BioPerl module for Bio::Taxonomy::Taxon
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Dan Kortschak but pilfered extensively from 
# the Bio::Tree::Node code of Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Taxonomy::Taxon - Generic Taxonomic Entity object

=head1 SYNOPSIS

    # NB: This module is deprecated. Use Bio::Taxon instead.

    use Bio::Taxonomy::Taxon;
    my $taxonA = Bio::Taxonomy::Taxon->new();
    my $taxonL = Bio::Taxonomy::Taxon->new();
    my $taxonR = Bio::Taxonomy::Taxon->new();

    my $taxon = Bio::Taxonomy::Taxon->new();
    $taxon->add_Descendents($taxonL);
    $taxon->add_Descendents($taxonR);

    my $species = $taxon->species;

=head1 DESCRIPTION

Makes a taxonomic unit suitable for use in a taxonomic tree

=head1 AUTHOR

Dan Kortschak email B<kortschak@rsbs.anu.edu.au>

=head1 CONTRIBUTORS

Sendu Bala: bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# code begins...

package Bio::Taxonomy::Taxon;
use vars qw($CREATIONORDER);
use strict;

use Bio::Species;

use base qw(Bio::Root::Root Bio::Tree::NodeI);

BEGIN { 
    $CREATIONORDER = 0;
}

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Taxonomy::Taxon->new();
 Function: Builds a new Bio::Taxonomy::Taxon object
 Returns : Bio::Taxonomy::Taxon
 Args    : -descendents   => array pointer to descendents (optional)
     	   -branch_length => branch length [integer] (optional)
     	   -taxon     => taxon
           -id     => unique taxon id for node (from NCBI's list preferably)
           -rank  => the taxonomic level of the node (also from NCBI)

=cut

#' for emacs

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  $self->warn("Bio::Taxonomy::Taxon is deprecated. Use Bio::Taxon instead.");
  
  my ($children,$branchlen,$id,$taxon,$rank,$desc) = 

                      $self->_rearrange([qw(DESCENDENTS
                                            BRANCH_LENGTH
                                            ID
                                            TAXON
                                            RANK
                                            DESC)], @args);
  
  $self->{'_desc'} = {};
  defined $desc && $self->description($desc);
  defined $taxon && $self->taxon($taxon);
  defined $id && $self->id($id);
  defined $branchlen && $self->branch_length($branchlen);
  defined $rank && $self->rank($rank);

  if( defined $children ) {
      if( ref($children) !~ /ARRAY/i ) {
	  $self->warn("Must specify a valid ARRAY reference to initialize a Taxon's Descendents");
      }
      foreach my $c ( @$children ) { 	
 	  $self->add_Descendent($c);
      }
  }
  $self->_creation_id($CREATIONORDER++);
  return $self;
}

=head2 add_Descendent

 Title   : add_Descendent
 Usage   : $taxon->add_Descendent($taxon);
 Function: Adds a descendent to a taxon
 Returns : number of current descendents for this taxon
 Args    : Bio::Taxonomy::Taxon
           boolean flag, true if you want to ignore the fact that you are
           adding a second node with the same unique id (typically memory 
           location reference in this implementation).  default is false and 
           will throw an error if you try and overwrite an existing node.

=cut

sub add_Descendent{

   my ($self,$node,$ignoreoverwrite) = @_;

   return -1 if( ! defined $node ) ;
   if( ! $node->isa('Bio::Taxonomy::Taxon') ) {
       $self->warn("Trying to add a Descendent who is not a Bio::Taxonomy::Taxon");
       return -1;
   }
   # do we care about order?
   $node->{'_ancestor'} = $self;
   if( $self->{'_desc'}->{$node->internal_id} && ! $ignoreoverwrite ) {
       $self->throw("Going to overwrite a taxon which is $node that is already stored here, set the ignore overwrite flag (parameter 2) to true to ignore this in the future");
   }
   
   $self->{'_desc'}->{$node->internal_id} = $node; # is this safely unique - we've tested before at any rate??
   
   $self->invalidate_height();
   
   return scalar keys %{$self->{'_desc'}};
}

=head2 each_Descendent

 Title   : each_Descendent($sortby)
 Usage   : my @taxa = $taxon->each_Descendent;
 Function: all the descendents for this taxon (but not their descendents
					      i.e. not a recursive fetchall)
 Returns : Array of Bio::Taxonomy::Taxon objects
 Args    : $sortby [optional] "height", "creation" or coderef to be used
           to sort the order of children taxa.

=cut

sub each_Descendent{
   my ($self, $sortby) = @_;

   # order can be based on branch length (and sub branchlength)

   $sortby ||= 'height';

   if (ref $sortby eq 'CODE') {
       my @values = sort $sortby values %{$self->{'_desc'}};
       return @values;
   } else  {
       if ($sortby eq 'height') {
	   return map { $_->[0] }
		  sort { $a->[1] <=> $b->[1] || 
			 $a->[2] <=> $b->[2] } 
	       map { [$_, $_->height, $_->internal_id ] } 
	   values %{$self->{'_desc'}};
       } else {
	   return map { $_->[0] }
	          sort { $a->[1] <=> $b->[1] } 
	          map { [$_, $_->height ] }
	          values %{$self->{'_desc'}};	   
       }
   }
}

=head2 remove_Descendent

 Title   : remove_Descendent
 Usage   : $taxon->remove_Descedent($taxon_foo);
 Function: Removes a specific taxon from being a Descendent of this taxon
 Returns : nothing
 Args    : An array of Bio::taxonomy::Taxon objects which have be previously
           passed to the add_Descendent call of this object.

=cut

sub remove_Descendent{
   my ($self,@nodes) = @_;
   foreach my $n ( @nodes ) { 
       if( $self->{'_desc'}->{$n->internal_id} ) {
	   $n->{'_ancestor'} = undef;
	   $self->{'_desc'}->{$n->internal_id}->{'_ancestor'} = undef;
	   delete $self->{'_desc'}->{$n->internal_id};
	   
       } else { 
	   $self->debug(sprintf("no taxon %s (%s) listed as a descendent in this taxon %s (%s)\n",$n->id, $n,$self->id,$self));
	   $self->debug("Descendents are " . join(',', keys %{$self->{'_desc'}})."\n");
       }
   }
   1;
}

=head2 remove_all_Descendents

 Title   : remove_all_Descendents
 Usage   : $taxon->remove_All_Descendents()
 Function: Cleanup the taxon's reference to descendents and reset
           their ancestor pointers to undef, if you don't have a reference
           to these objects after this call they will be cleanedup - so
           a get_nodes from the Tree object would be a safe thing to do first
 Returns : nothing
 Args    : none

=cut

sub remove_all_Descendents{
   my ($self) = @_;
   # this won't cleanup the taxa themselves if you also have
   # a copy/pointer of them (I think)...
   while( my ($node,$val) = each %{ $self->{'_desc'} } ) {
       $val->{'_ancestor'} = undef;
   }
   $self->{'_desc'} = {};
   1;
}

=head2 get_Descendents

 Title   : get_Descendents
 Usage   : my @taxa = $taxon->get_Descendents;
 Function: Recursively fetch all the taxa and their descendents
           *NOTE* This is different from each_Descendent
 Returns : Array or Bio::Taxonomy::Taxon objects
 Args    : none

=cut

# implemented in the interface 

=head2 ancestor

 Title   : ancestor
 Usage   : $taxon->ancestor($newval)
 Function: Set the Ancestor
 Returns : value of ancestor
 Args    : newvalue (optional)

=cut

sub ancestor {
   my ($self, $value) = @_;
   if (defined $value) {
       $self->{'_ancestor'} = $value;
   }
   return $self->{'_ancestor'};
}

=head2 branch_length

 Title   : branch_length
 Usage   : $obj->branch_length($newval)
 Function:
 Example :
 Returns : value of branch_length
 Args    : newvalue (optional)

=cut

sub branch_length {
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'branch_length'} = $value;
    }
    return $self->{'branch_length'};
}

=head2 description

 Title   : description
 Usage   : $obj->description($newval)
 Function:
 Returns : value of description
 Args    : newvalue (optional)

=cut

sub description {
   my ($self,$value) = @_;
   if( defined $value  ) {
       $self->{'_description'} = $value;
   }
   return $self->{'_description'};
}

=head2 rank

 Title   : rank
 Usage   : $obj->rank($newval)
 Function: Set the taxonomic rank
 Returns : taxonomic rank of taxon
 Args    : newvalue (optional)

=cut

sub rank {
   my ($self,$value) = @_;
   if (defined $value) {
      $self->{'_rank'} = $value;
   }
   return $self->{'_rank'};
}

=head2 taxon

 Title   : taxon
 Usage   : $obj->taxon($newtaxon)
 Function: Set the name of the taxon
 Example :
 Returns : name of taxon
 Args    : newtaxon (optional)

=cut

# because internal taxa have names too...
sub taxon {
   my ($self,$value) = @_;
   if( defined $value  ) {
       $self->{'_taxon'} = $value;
   }
   return $self->{'_taxon'};
}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function:
 Example :
 Returns : value of id
 Args    : newvalue (optional)

=cut

sub id {
   my ($self,$value) = @_;
   if( defined $value ) {
       $self->{'_id'} = $value;
   }
   return $self->{'_id'};
}

sub DESTROY {
    my ($self) = @_;
    # try to insure that everything is cleaned up
    $self->SUPER::DESTROY();
    if( defined $self->{'_desc'} &&
	ref($self->{'_desc'}) =~ /ARRAY/i ) {
	while( my ($nodeid,$node) = each %{ $self->{'_desc'} } ) {
	    $node->{'_ancestor'} = undef; # ensure no circular references
	    $node->DESTROY();
	    $node = undef;
	}
	$self->{'_desc'} = {};
    }
}

=head2 internal_id

 Title   : internal_id
 Usage   : my $internalid = $taxon->internal_id
 Function: Returns the internal unique id for this taxon
           (a monotonically increasing number for this in-memory implementation
            but could be a database determined unique id in other 
	    implementations)
 Returns : unique id
 Args    : none

=cut

sub internal_id {
   return $_[0]->_creation_id;
}

=head2 _creation_id

 Title   : _creation_id
 Usage   : $obj->_creation_id($newval)
 Function: a private method signifying the internal creation order
 Returns : value of _creation_id
 Args    : newvalue (optional)


=cut

sub _creation_id {
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'_creation_id'} = $value;
    }
    return $self->{'_creation_id'} || 0;
}

# The following methods are implemented by NodeI decorated interface

=head2 is_Leaf

 Title   : is_Leaf
 Usage   : if( $node->is_Leaf )
 Function: Get Leaf status
 Returns : boolean
 Args    : none

=cut

sub is_Leaf {
    my ($self) = @_;
    my $rc = 0;
    $rc = 1 if( ! defined $self->{'_desc'} ||	
		keys %{$self->{'_desc'}} == 0);
    return $rc;
}

=head2 to_string

 Title   : to_string
 Usage   : my $str = $taxon->to_string()
 Function: For debugging, provide a taxon as a string
 Returns : string
 Args    : none

=cut

=head2 height

 Title   : height
 Usage   : my $len = $taxon->height
 Function: Returns the height of the tree starting at this
           taxon.  Height is the maximum branchlength.
 Returns : The longest length (weighting branches with branch_length) to a leaf
 Args    : none

=cut

sub height { 
    my ($self) = @_;

    return $self->{'_height'} if( defined $self->{'_height'} );
    
    if( $self->is_Leaf ) { 
      if( !defined $self->branch_length ) { 
	      $self->debug(sprintf("Trying to calculate height of a taxon when a taxon (%s) has an undefined branch_length",$self->id || '?' ));
	      return 0;
      }
      return $self->branch_length;
   }
   my $max = 0;
   foreach my $subnode ( $self->each_Descendent ) { 
       my $s = $subnode->height;
       if( $s > $max ) { $max = $s; }
   }
   return ($self->{'_height'} = $max + ($self->branch_length || 1));
}

=head2 invalidate_height

 Title   : invalidate_height
 Usage   : private helper method
 Function: Invalidate our cached value of the taxon's height in the tree
 Returns : nothing
 Args    : none

=cut

sub invalidate_height { 
    my ($self) = @_;
    
    $self->{'_height'} = undef;
    if( $self->ancestor ) {
	    $self->ancestor->invalidate_height;
    }
}

=head2 classify

 Title   : classify
 Usage   : @obj->classify()
 Function: a method to return the classification of a species
 Returns : name of taxon and ancestor's taxon recursively
 Args    : boolean to specify whether we want all taxa not just ranked 
           levels

=cut

sub classify {
   my ($self,$allnodes) = @_;

   my @classification=($self->taxon);
   my $node=$self;

   while (defined $node->ancestor) {
      push @classification, $node->ancestor->taxon if $allnodes==1;
      $node=$node->ancestor;
   }

   return (@classification);
}

=head2 has_rank

 Title   : has_rank
 Usage   : $obj->has_rank($rank)
 Function: a method to query ancestors' rank
 Returns : boolean
 Args    : $rank

=cut

sub has_rank {
   my ($self,$rank) = @_;

   return $self if $self->rank eq $rank;

   while (defined $self->ancestor) {
      return $self if $self->ancestor->rank eq $rank;
      $self=$self->ancestor;
   }

   return;
}

=head2 has_taxon

 Title   : has_taxon
 Usage   : $obj->has_taxon($taxon)
 Function: a method to query ancestors' taxa
 Returns : boolean
 Args    : Bio::Taxonomy::Taxon object

=cut

sub has_taxon {
   my ($self,$taxon) = @_;

   return $self if 
      ((defined $self->id && $self->id == $taxon->id) ||
      ($self->taxon eq $taxon->taxon && $self->rank eq $taxon->rank));

   while (defined $self->ancestor) {
      return $self if 
         ((defined $self->id && $self->id == $taxon->id) ||
         ($self->taxon eq $taxon->taxon && $self->rank eq $taxon->rank) &&
         ($self->taxon ne 'no rank'));
      $self=$self->ancestor;
   }

   return;
}

=head2 distance_to_root

 Title   : distance_to_root
 Usage   : $obj->distance_to_root
 Function: a method to query ancestors' taxa
 Returns : number of links to root
 Args    :

=cut

sub distance_to_root {
   my ($self,$taxon) = @_;

   my $count=0;

   while (defined $self->ancestor) {
      $count++;
      $self=$self->ancestor;
   }

   return $count;
}

=head2 recent_common_ancestor

 Title   : recent_common_ancestor
 Usage   : $obj->recent_common_ancestor($taxon)
 Function: a method to query find common ancestors
 Returns : Bio::Taxonomy::Taxon of query or undef if no ancestor of rank
 Args    : Bio::Taxonomy::Taxon

=cut

sub recent_common_ancestor {
   my ($self,$node) = @_;

   while (defined $node->ancestor) {
      my $common=$self->has_taxon($node);
      return $common if defined $common;
      $node=$node->ancestor;
   }

   return;
}

=head2 species

 Title   : species
 Usage   : $obj=$taxon->species;
 Function: Returns a Bio::Species object reflecting the taxon's tree position
 Returns : a Bio::Species object
 Args    : none

=cut

sub species {
   my ($self) = @_;
   my $species;

   if ($self->has_rank('subspecies') && $self->ancestor->rank eq 'species') {
      $species = Bio::Species->new(-classification => $self->ancestor->classify);
      $species->genus($self->ancestor->ancestor->taxon);
      $species->species($self->ancestor->taxon);
      $species->sub_species($self->taxon);
   } elsif ($self->has_rank('species')) {
      $species = Bio::Species->new(-classification => $self->classify);
      $species->genus($self->ancestor->taxon);
      $species->species($self->taxon);
   } else {
      $self->throw("Trying to create a species from a taxonomic entity without species rank. Use classify instead of species.\n");
   }
   return $species;
}

1;
