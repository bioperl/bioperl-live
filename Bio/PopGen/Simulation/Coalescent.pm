#
# BioPerl module for Bio::PopGen::Simulation::Coalescent
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

Bio::PopGen::Simulation::Coalescent - A Coalescent simulation factory

=head1 SYNOPSIS

    use Bio::PopGen::Simulation::Coalescent;
    my @taxonnames = qw(SpeciesA SpeciesB SpeciesC SpeciesD);
    my $sim1 = Bio::PopGen::Simulation::Coalescent->new(-samples => \@taxonnames);

    my $tree = $sim1->next_tree;

    # add 20 mutations randomly to the tree
    $sim1->add_Mutations($tree,20);

    # or for anonymous samples

    my $sim2 = Bio::PopGen::Simulation::Coalescent->new( -sample_size => 6,
							 -maxcount => 50);
    my $tree2 = $sim2->next_tree;
    # add 20 mutations randomly to the tree
    $sim2->add_Mutations($tree2,20);

=head1 DESCRIPTION

Builds a random tree every time next_tree is called or up to -maxcount
times with branch lengths and provides the ability to randomly add
mutations onto the tree with a probabilty proportional to the branch
lengths.

This algorithm is based on the make_tree algorithm from Richard Hudson 1990.

Hudson, R. R. 1990. Gene genealogies and the coalescent
       process. Pp. 1-44 in D. Futuyma and J.  Antonovics, eds. Oxford
       surveys in evolutionary biology. Vol. 7. Oxford University
       Press, New York.

This module was previously named Bio::Tree::RandomTree

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

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich, Matthew Hahn

Email jason-at-bioperl-dot-org
Email matthew-dot-hahn-at-duke-dot-edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::PopGen::Simulation::Coalescent;
use vars qw($PRECISION_DIGITS);
use strict;

$PRECISION_DIGITS = 3; # Precision for the branchlength

use Bio::Tree::AlleleNode;
use Bio::PopGen::Genotype;
use Bio::Tree::Tree;

use base qw(Bio::Root::Root Bio::Factory::TreeFactoryI);


=head2 new

 Title   : new
 Usage   : my $obj = Bio::PopGen::Simulation::Coalescent->new();
 Function: Builds a new Bio::PopGen::Simulation::Coalescent object 
 Returns : an instance of Bio::PopGen::Simulation::Coalescent
 Args    : -samples => arrayref of sample names
           OR
           -sample_size=> number of samples (samps will get a systematic name)
           -maxcount   => [optional] maximum number of trees to provide

=cut

sub new{
   my ($class,@args) = @_;
   my $self = $class->SUPER::new(@args);
   
   $self->{'_treecounter'} = 0;
   $self->{'_maxcount'} = 0;
   my ($maxcount, $samps,$samplesize ) = $self->_rearrange([qw(MAXCOUNT
							       SAMPLES
							       SAMPLE_SIZE)],
							   @args);
   my @samples;
   
   if( ! defined $samps ) { 
       if( ! defined $samplesize || $samplesize <= 0 ) { 
	   $self->throw("Must specify a valid samplesize if parameter -SAMPLE is not specified (sampsize is $samplesize)");
       }
       foreach ( 1..$samplesize ) { push @samples, "Samp$_"; }      
   } else { 
       if( ref($samps) !~ /ARRAY/i ) { 
	   $self->throw("Must specify a valid ARRAY reference to the parameter -SAMPLES, did you forget a leading '\\'?");
       }
       @samples = @$samps;
   }
   
   $self->samples(\@samples);
   $self->sample_size(scalar @samples);
   defined $maxcount && $self->maxcount($maxcount);   
   return $self;
}

=head2 next_tree

 Title   : next_tree
 Usage   : my $tree = $factory->next_tree
 Function: Returns a random tree based on the initialized number of nodes
           NOTE: if maxcount is not specified on initialization or
                 set to a valid integer, subsequent calls to next_tree will 
                 continue to return random trees and never return undef
 Returns : Bio::Tree::TreeI object
 Args    : none

=cut

sub next_tree{
   my ($self) = @_;
   # If maxcount is set to something non-zero then next tree will
   # continue to return valid trees until maxcount is reached
   # otherwise will always return trees 
   return if( $self->maxcount &&
	      $self->{'_treecounter'}++ >= $self->maxcount );
   my $size = $self->sample_size;
   
   my $in;
   my @tree = ();
   my @list = ();
   
   for($in=0;$in < 2*$size -1; $in++ ) { 
       push @tree, { 'nodenum' => "Node$in" };
   }
   # in C we would have 2 arrays
   # an array of nodes (tree)
   # and array of pointers to these nodes (list)
   # and we just shuffle the list items to do the 
   # tree topology generation
   # instead in perl, we will have a list of hashes (nodes) called @tree
   # and a list of integers representing the indexes in tree called @list

   for($in=0;$in < $size;$in++)  {
       $tree[$in]->{'time'} = 0;
       $tree[$in]->{'desc1'} = undef;
       $tree[$in]->{'desc2'} = undef;
       push @list, $in;
   }

   my $t=0;
   # generate times for the nodes
   for($in = $size; $in > 1; $in-- ) {
	$t+= -2.0 * log(1 - $self->random(1)) / ( $in * ($in-1) );    
	$tree[2 * $size - $in]->{'time'} =$t;
    }
   # topology generation
   for ($in = $size; $in > 1; $in-- ) {
       my $pick = int $self->random($in);    
       my $nodeindex = $list[$pick];       
       my $swap = 2 * $size - $in;       
       $tree[$swap]->{'desc1'} = $nodeindex;	
       $list[$pick] = $list[$in-1];       
       $pick = int rand($in - 1);    
       $nodeindex = $list[$pick];
       $tree[$swap]->{'desc2'} = $nodeindex;	
       $list[$pick] = $swap;
   }
   # Let's convert the hashes into nodes

   my @nodes = ();   
   foreach my $n ( @tree ) { 
       push @nodes, 
	   Bio::Tree::AlleleNode->new(-id => $n->{'nodenum'},
				     -branch_length => $n->{'time'});
   }
   my $ct = 0;
   foreach my $node ( @nodes ) { 
       my $n = $tree[$ct++];
       if( defined $n->{'desc1'} ) {
	   $node->add_Descendent($nodes[$n->{'desc1'}]);
       }
       if( defined $n->{'desc2'} ) { 
	   $node->add_Descendent($nodes[$n->{'desc2'}]);
       }
   }   
   my $T = Bio::Tree::Tree->new(-root => pop @nodes );
   return $T;
}

=head2 add_Mutations

 Title   : add_Mutations
 Usage   : $factory->add_Mutations($tree, $mutcount);
 Function: Adds mutations to a tree via a random process weighted by 
           branch length (it is a poisson distribution 
			  as part of a coalescent process) 
 Returns : none
 Args    : $tree - Bio::Tree::TreeI 
           $nummut - number of mutations
           $precision - optional # of digits for precision


=cut

sub add_Mutations{
   my ($self,$tree, $nummut,$precision) = @_;
   $precision ||= $PRECISION_DIGITS;
   $precision = 10**$precision;

   my @branches;
   my @lens;
   my $branchlen = 0;
   my $last = 0;
   my @nodes = $tree->get_nodes();
   my $i = 0;

   # Jason's somewhat simplistics way of doing a poission
   # distribution for a fixed number of mutations
   # build an array and put the node number in a slot
   # representing the branch to put a mutation on
   # but weight the number of slots per branch by the 
   # length of the branch ( ancestor's time - node time)
   
   foreach my $node ( @nodes ) {
       if( $node->ancestor ) { 
	   my $len = int ( ($node->ancestor->branch_length - 
			    $node->branch_length) * $precision);
	   if ( $len > 0 ) {
	       for( my $j =0;$j < $len;$j++) {
		   push @branches, $i;
	       }
	       $last += $len;
	   }
	   $branchlen += $len;
       }
       if( ! $node->isa('Bio::Tree::AlleleNode') ) {
	   bless $node, 'Bio::Tree::AlleleNode'; # rebless it to the right node
       } 
       # This let's us reset the stored genotypes so we can keep reusing the 
       # same tree topology, but throw down mutations multiple times
       $node->reset_Genotypes;
       $i++;
   }
   # sanity check
   $self->throw("branch len is $branchlen arraylen is $last")
        unless ( $branchlen == $last );
   my @mutations;
   for( my $j = 0; $j < $nummut; $j++)  {
       my $index = int(rand($branchlen));
       my $branch = $branches[$index];

       # We're using an infinite sites model so every new
       # mutation is a new site
       my $g = Bio::PopGen::Genotype->new(-marker_name  => "Mutation$j",
					 -alleles => [1]);
       $nodes[$branch]->add_Genotype($g);
       push @mutations, "Mutation$j";
       # Let's add this mutation to all the children (push it down
       # the branches to the tips)
       foreach my $child ( $nodes[$branch]->get_all_Descendents ) {
	   $child->add_Genotype($g);
       }
   }
   # Insure that everyone who doesn't have the mutation
   # has the ancestral state, which is '0'
   foreach my $node ( @nodes ) {
       foreach my $m ( @mutations ) {
	   if( ! $node->has_Marker($m) ) {
	       my $emptyg = Bio::PopGen::Genotype->new(-marker_name => $m,
						      -alleles     => [0]);
	       $node->add_Genotype($emptyg);
	   }
       }
   }
}

=head2 maxcount

 Title   : maxcount
 Usage   : $obj->maxcount($newval)
 Function: 
 Returns : Maxcount value
 Args    : newvalue (optional)


=cut

sub maxcount{
   my ($self,$value) = @_;
   if( defined $value) {
       if( $value =~ /^(\d+)/ ) { 
	   $self->{'maxcount'} = $1;
       } else { 
	   $self->warn("Must specify a valid Positive integer to maxcount");
	   $self->{'maxcount'} = 0;
       }
  }
   return $self->{'_maxcount'};
}

=head2 samples

 Title   : samples
 Usage   : $obj->samples($newval)
 Function: 
 Example : 
 Returns : value of samples
 Args    : newvalue (optional)


=cut

sub samples{
   my ($self,$value) = @_;
   if( defined $value) {
       if( ref($value) !~ /ARRAY/i ) { 
	   $self->warn("Must specify a valid array ref to the method 'samples'");
	   $value = [];
       } 
      $self->{'samples'} = $value;
    }
    return $self->{'samples'};

}

=head2 sample_size

 Title   : sample_size
 Usage   : $obj->sample_size($newval)
 Function: 
 Example : 
 Returns : value of sample_size
 Args    : newvalue (optional)


=cut

sub sample_size{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'sample_size'} = $value;
    }
    return $self->{'sample_size'};

}

=head2 random

 Title   : random
 Usage   : my $rfloat = $node->random($size)
 Function: Generates a random number between 0 and $size
           This is abstracted so that someone can override and provide their
           own special RNG.  This is expected to be a uniform RNG.
 Returns : Floating point random
 Args    : $maximum size for random number (defaults to 1)


=cut

sub random{
   my ($self,$max) = @_;
   return rand($max);
}

1;
