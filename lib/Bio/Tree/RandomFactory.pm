#
# BioPerl module for Bio::Tree::RandomFactory
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::RandomFactory - TreeFactory for generating Random Trees

=head1 SYNOPSIS

  use Bio::Tree::RandomFactory
  my @taxonnames;
  my $factory = Bio::Tree::RandomFactory->new( -taxa => \@taxonnames,
  					      -maxcount => 10);

  # or for anonymous samples

  my $factory = Bio::Tree::RandomFactory->new( -num_taxa => 6,
					      -maxcount => 50);


  my $tree = $factory->next_tree;

=head1 DESCRIPTION

Builds a random tree every time next_tree is called or up to -maxcount times.

This module was originally written for Coalescent simulations see
L<Bio::PopGen::Simulation::Coalescent>.  I've left the next_tree
method intact although it is not generating random trees in the
phylogenetic sense.  I would be happy for someone to provide
alternative implementations which can be used here.  As written it
will generate random topologies but the branch lengths are built from
assumptions in the coalescent and are not appropriate for phylogenetic
analyses.

This algorithm is based on the make_tree algorithm from Richard Hudson 1990.

Hudson, R. R. 1990. Gene genealogies and the coalescent
       process. Pp. 1-44 in D. Futuyma and J.  Antonovics, eds. Oxford
       surveys in evolutionary biology. Vol. 7. Oxford University
       Press, New York

Sanderson, M ... 

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

=head1 AUTHOR - Jason Stajich

Email jason-AT-bioperl.org

=head1 CONTRIBUTORS

Matthew Hahn, E<lt>matthew.hahn@duke.eduE<gt>
Mike Sanderson 

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tree::RandomFactory;
use vars qw($PRECISION_DIGITS $DefaultNodeType %Defaults);
use strict;

$PRECISION_DIGITS = 3; # Precision for the branchlength
$DefaultNodeType = 'Bio::Tree::Node';
%Defaults = ('YuleRate'          => 1.0, # as set by Sanderson in Rates
	     'Speciation'        => 1.0, #
	     'DefaultTreeMethod' => 'yule',
	     );

use Bio::Tools::RandomDistFunctions;
use Bio::Tree::Tree;

use base qw(Bio::Root::Root Bio::Factory::TreeFactoryI);

=head2 new

 Title   : new
 Usage   : my $factory = Bio::Tree::RandomFactory->new(-samples => \@samples,
						      -maxcount=> $N);
 Function: Initializes a Bio::Tree::RandomFactory object
 Returns : Bio::Tree::RandomFactory
 Args    : -nodetype => Type of Nodes to create [default Bio::Tree::Node]
           -maxcount => [optional] Maximum num trees to create
           -randtype => Type of random trees so far support
               - yule/backward_yule/BY [default]
               - forward_yule/FY
               - birthdeath_forward/BDF
               - birthdeath_backwards/BDB


          ONE of the following must be specified
           -taxa     => $arrayref of taxa names
           -num_taxa => integer indicating number of taxa in the tree

=cut

sub new{
   my ($class,@args) = @_;
   my $self = $class->SUPER::new(@args);
   
   $self->{'_treecounter'} = 0;
   $self->{'_maxcount'} = 0;
   my ($nodetype,$randtype,
       $maxcount, $samps,$samplesize,
       $taxa, $num_taxa) = $self->_rearrange([qw(NODETYPE
						 RANDTYPE
						 MAXCOUNT
						 SAMPLES
						 SAMPLE_SIZE
						 TAXA
						 NUM_TAXA)],
					     @args);
   my @taxa;
   $nodetype ||= $DefaultNodeType;
   $self->nodetype($nodetype);
   $taxa = $samps if defined $samps && ! defined $taxa;
   $num_taxa = $samplesize if $samplesize && ! $num_taxa;
   if( ! defined $taxa ) { 
       if( ! defined $num_taxa || $num_taxa <= 0 ) { 
	   $self->throw("Must specify a valid num_taxa if parameter -TAXA is not specified");
       }
       foreach ( 1..$num_taxa ) { push @taxa, "Taxon$_"; }      
   } else { 
       if( ref($taxa) !~ /ARRAY/i ) { 
	   $self->throw("Must specify a valid ARRAY reference to the parameter -TAXA, did you forget a leading '\\'? for $taxa");
       }
       @taxa = @$taxa;
   }
   
   $self->taxa(\@taxa);
   defined $maxcount && $self->maxcount($maxcount);   
   $self->{'_count'} = 0;
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
   my ($self,%options) = @_;
   return if $self->maxcount && 
       $self->{'_count'}++ >= $self->maxcount;
   my $rand_type = $options{'randtype'} || $self->random_tree_method;
   my $nodetype = $self->nodetype;
   my $treearray;

   if( $rand_type =~ /(birthdeath_forward|birth|BDF)/i ) {

   } elsif ( $rand_type =~ /(birthdeath_backward|BDB)/i ) {
       $treearray = $self->rand_birthdeath_backwards_tree;       
   } elsif( $rand_type =~ /(BY|backwards_yule)/i || 
	    $rand_type =~ /^yule/i ) {
       my $speciation = $options{'speciation'}; # can be undef
       $treearray = $self->rand_yule_c_tree($speciation);       
   } else { 
       $self->warn("unrecognized random type $rand_type");
   }
   
   my @nodes = ();   
   foreach my $n ( @$treearray ) { 
       for my $k ( qw(desc1 desc2) ) {
	   next unless defined $n->{$k};
	   push @{$n->{'descendents'}}, $nodes[$n->{$k}];
       }
       push @nodes, 
       $nodetype->new(-id            => $n->{'nodenum'},
		      -branch_length => $n->{'time'},
		      -descendents   => $n->{'descendents'},
		      );
   }
   my $T = Bio::Tree::Tree->new(-root => pop @nodes );
   return $T;
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
	   $self->{'_maxcount'} = $1;
       } else { 
	   $self->warn("Must specify a valid Positive integer to maxcount");
	   $self->{'_maxcount'} = 0;
       }
  }
   return $self->{'_maxcount'};
}


=head2 reset_tree_count

 Title   : reset_tree_count
 Usage   : $factory->reset_tree_count;
 Function: Reset the tree counter
 Returns : none
 Args    : none


=cut

sub reset_count{
    shift->{'_count'} = 0;
}

=head2 taxa

 Title   : taxa
 Usage   : $obj->taxa($newval)
 Function: Set the leaf node names
 Returns : value of taxa
 Args    : Arrayref of Taxon names


=cut

sub taxa {
    my ($self,$value) = @_;
    if( defined $value) {
	if( ref($value) !~ /ARRAY/i ) { 
	    $self->warn("Must specify a valid array ref to the method 'taxa'");
	    $value = [];
	} 
	$self->{'_taxa'} = $value;
	$self->{'_num_taxa'} = scalar @$value;
    }
    return $self->{'_taxa'};

}

=head2 num_taxa

 Title   : num_taxa
 Usage   : $obj->num_taxa($newval)
 Function: Get the number of Taxa
 Returns : value of num_taxa
 Args    : none


=cut

sub num_taxa {
    my ($self) = @_;
    return  $self->{'_num_taxa'};
}

# alias old methods
*num_samples = \&num_taxa;
*samples = \&taxa;

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


=head2 random_tree_method

 Title   : random_tree_method
 Usage   : $obj->random_tree_method($newval)
 Function: 
 Example : 
 Returns : value of random_tree_method (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub random_tree_method{
    my $self = shift;

    return $self->{'random_tree_method'} = shift if @_;
    return $self->{'random_tree_method'} || $Defaults{'DefaultTreeMethod'};
}

=head2 nodetype

 Title   : nodetype
 Usage   : $obj->nodetype($newval)
 Function: 
 Example : 
 Returns : value of nodetype (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub nodetype{
   my ($self,$value) = @_;
   if( defined $value) {
       eval "require $value";
       if( $@ ) { $self->throw("$@: Unrecognized Node type for ".ref($self). 
			       "'$value'");}
       
       my $a = bless {},$value;
       unless( $a->isa('Bio::Tree::NodeI')  ) {
	   $self->throw("Must provide a valid Bio::Tree::NodeI or child class to SeqFactory Not $value");
       }
      $self->{'nodetype'} = $value;
    }
    return $self->{'nodetype'};
}


# The assignment of times are based on Mike Sanderson's r8s code
# The topology assignment code is based on Richard Hudson's
# make_trees


sub rand_yule_c_tree {
    my ($self,$speciation) = @_;
    $speciation ||= $Defaults{'Speciation'};
    my $n_taxa = $self->num_taxa;
    my $taxa = $self->taxa || [];
    my $nodetype = $self->nodetype;
  
    my $randfuncs = Bio::Tools::RandomDistFunctions->new();
    my $rate = $Defaults{'YuleRate'};
    my (@tree,@list,@times,$i,$in);
    my $max = 2 * $n_taxa - 1;
    for($in=0;$in < $max; $in++ ) { 
	push @tree, { 'nodenum' => "Node$in" };
    }
    # setup leaf nodes
    for($in=0;$in < $n_taxa;$in++)  {
	$tree[$in]->{'time'} = 0;
	$tree[$in]->{'desc1'} = undef;
	$tree[$in]->{'desc2'} = undef;
	if( my $r = $taxa->[$in] ) { 
	    $tree[$in]->{'nodenum'} = $r;
	}
	push @list, $in;
    }
    
    for( $i = 0; $i < $n_taxa - 1; $i++ ) {
	# draw random interval times
	push @times, $randfuncs->rand_birth_distribution($speciation);
    }
    # sort smallest to largest
    @times = sort {$a <=> $b} @times;
     # topology generation
    for ($in = $n_taxa; $in > 1; $in-- ) {
	my $time = shift @times;
	
	my $pick = int $self->random($in);    
	my $nodeindex = $list[$pick];
	$tree[$list[$pick]]->{'time'} = $time;
	my $swap = 2 * $n_taxa - $in;
	$tree[$swap]->{'desc1'} = $nodeindex;	
	$list[$pick] = $list[$in-1];       
	
	$pick = int rand($in - 1);    
	$nodeindex = $list[$pick];
	$tree[$list[$pick]]->{'time'} = $time;
	$tree[$swap]->{'desc2'} = $nodeindex;	
	$list[$pick] = $swap;	
    }
    $tree[-1]->{'time'} = shift @times;
    return \@tree;
}



sub rand_birthdeath_backwards_tree {
    my ($self) = @_;
    my $n_taxa = $self->num_taxa;
    my $taxa = $self->taxa || [];
  
    my $randfuncs = Bio::Tools::RandomDistFunctions->new();
    my $rate = $Defaults{'YuleRate'};
    my (@tree,@list,@times,$i,$in);
    my $max = 2 * $n_taxa - 1;
    for($in=0;$in < $max; $in++ ) { 
	push @tree, { 'nodenum' => "Node$in" };
    }
    # setup leaf nodes
    for($in=0;$in < $n_taxa;$in++)  {
	$tree[$in]->{'time'} = 0;
	$tree[$in]->{'desc1'} = undef;
	$tree[$in]->{'desc2'} = undef;
	if( my $r = $taxa->[$in] ) { 
	    # deal with pre-labeled nodes
	    $tree[$in]->{'nodenum'} = $r;
	}
	push @list, $in;
    }
    my ($time) = (0);

     # topology generation
    for ($in = $n_taxa; $in > 1; $in-- ) {
	my $pick = int $self->random($in);    
	my $nodeindex = $list[$pick];
	my $swap = 2 * $n_taxa - $in;
	$time += $randfuncs->rand_geometric_distribution($n_taxa * $rate);;
	$tree[$list[$pick]]->{'time'} = $time;
	$tree[$swap]->{'desc1'} = $nodeindex;	
	$list[$pick] = $list[$in-1];       
	
	$pick = int rand($in - 1);    
	$nodeindex = $list[$pick];
	$tree[$list[$pick]]->{'time'} = $time;
	$tree[$swap]->{'desc2'} = $nodeindex;	
	$list[$pick] = $swap;	
    }
    my $root = $tree[-1];
    $time += $randfuncs->rand_geometric_distribution($n_taxa * $rate);;
    $root->{'time'} = $time;

    # Normalize times by the root node...
    for my $node ( @tree ) {
	$node->{'time'} /= $root->{'time'};
    }
    return \@tree;
}


# The assignment of times are based on Mike Sanderson's r8s code
# The topology assignment code is based on Richard Hudson's
# make_trees

sub rand_birth_death_tree {
# Still need to finish
#     my ($self,$spec_rate,$extinct_rate,$char_rate) = @_;
#     my $n_taxa =  $self->num_taxa;
#     my $dt = 0.1 / $n_taxa;
#     my @tree;
#     my $max = 3 * $n_taxa - 1;
#     # setup leaf nodes
    
#     for($in=0;$in < $size;$in++)  {
# 	push @tree, { 'nodenum' => $taxa->[$in] || "Node$in",
# 		      'time'    => 0,
# 		      'desc1'   => undef,
# 		      'desc2'   => undef, 
# 		  };
#     }
#     my $time = $dt;
#     my $idx = 0;
#     while( $n_taxa > 1 ) { 	
# 	if ( event($dt * $spec_rate, $n_taxa) ) {
# 	    my $pick = int $self->random($n_taxa);
# 	    my $pick2 = int $self->random($n_taxa);
# 	    while( $pick2 == $pick ) {
# 		$pick2 = int $self->random($n_taxa);
# 	    }
	    # to finish....
	    
# 	    $tree[$swap]->{'desc1'} = $nodeindex;		    
# 	}
#     }

	    

# 	$list[$pick] = $list[$in-1];       
	
# 	$pick = int rand($in - 1);    
# 	$nodeindex = $list[$pick];
# 	$tree[$swap]->{'desc2'} = $nodeindex;	
# 	$list[$pick] = $swap;	
# 	$tree[$swap]->{'time'} = $times[$ix++];
#     }
}


1;
