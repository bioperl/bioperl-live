#
# Module for Bio::PhyloNetwork::TreeFactoryMulti
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Gabriel Cardona <gabriel(dot)cardona(at)uib(dot)es>
#
# Copyright Gabriel Cardona
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::PhyloNetwork::TreeFactoryMulti - Module to sequentially generate
Phylogenetic Trees

=head1 SYNOPSIS

 use strict;
 use warnings;

 use Bio::PhyloNetwork;
 use Bio::PhyloNetwork::TreeFactory;

 # Will generate sequentially all the 15 binary phylogetic
 # trees with 4 leaves

 my $factory=Bio::PhyloNetwork::TreeFactory->new(-numleaves=>4);

 my @nets;

 while (my $net=$factory->next_network()) {
   push @nets,$net;
   print "".(scalar @nets).": ".$net->eNewick()."\n";
 }

=head1 DESCRIPTION

Sequentially builds a (binary) phylogenetic tree each time
next_network is called.

=head1 AUTHOR

Gabriel Cardona, gabriel(dot)cardona(at)uib(dot)es

=head1 SEE ALSO

L<Bio::PhyloNetwork>

=head1 APPENDIX

The rest of the documentation details each of the object methods.

=cut

package Bio::PhyloNetwork::TreeFactoryMulti;

use strict;
use warnings;

use base qw(Bio::Root::Root);

use Bio::PhyloNetwork;

=head2 new

 Title   : new
 Usage   : my $factory = new Bio::PhyloNetwork::TreeFactory();
 Function: Creates a new Bio::PhyloNetwork::TreeFactory
 Returns : Bio::PhyloNetwork::RandomFactory
 Args    : -numleaves => integer
            OR
           -leaves => reference to an array (of leaves names)

Returns a Bio::PhyloNetwork::TreeFactory object. Such an object will
sequentially create binary phylogenetic trees
each time next_network is called.

If the parameter -leaves=E<gt>\@leaves is given, then the set of leaves of
these networks will be @leaves. If it is given the parameter
-numleaves=E<gt>$numleaves, then the set of leaves will be "l1"..."l$numleaves".

=cut

sub new {
  my ($pkg,@args)=@_;

  my $self=$pkg->SUPER::new(@args);

  my ($leavesR,$numleaves,$numhybrids)=
    $self->_rearrange([qw(LEAVES
			  NUMLEAVES
			  NUMHYBRIDS)],@args);

  my @leaves;
  if ((! defined $leavesR) && (defined $numleaves)) {
    @leaves=map {"l$_"} (1..$numleaves);
    $leavesR=\@leaves;
  }
  if (! defined $leavesR) {
    $self->throw("No leaves set neither numleaves given");
  }
  @leaves=@$leavesR;
  $self->{leaves}=$leavesR;

  $numleaves=@leaves;
  $self->{numleaves}=$numleaves;
  if ($numleaves > 2) {
    my @leavesparent=@leaves;
    my $newleaf=pop @leavesparent;
    $self->{newleaf}=$newleaf;
    $self->{parent}=
      new($pkg,-leaves=>\@leavesparent);
    my $oldnet=$self->{parent}->next_network();
    $self->{oldnet}=$oldnet;
    my @candidates=$oldnet->nodes();
    $self->{candidates}=\@candidates;
    my @candidatesbis=$oldnet->internal_nodes();
    $self->{candidatesbis}=\@candidatesbis;
    $self->{processbis}=0;
  }
  $self->{index}=0;

  bless($self,$pkg);
}

=head2 next_network

 Title   : next_network
 Usage   : my $net=$factory->next_network()
 Function: returns a tree
 Returns : Bio::PhyloNetwork
 Args    : none

=cut

sub next_network {
  my ($self)=@_;

  my $n=$self->{numleaves};
  if ($self->{numleaves} == 2) {
    if ($self->{index} == 0) {
      my $graph=Graph::Directed->new();
      $graph->add_edges("t0",$self->{leaves}->[0],"t0",$self->{leaves}->[1]);
      my $net=Bio::PhyloNetwork->new(-graph=>$graph);
      $self->{index}++;
      return $net;
    }
    else {
      return 0;
    }
  }
  else {
    if (($self->{index} == (scalar @{$self->{candidatesbis}})) && $self->{processbis}
       ) {
      my $oldnet=$self->{parent}->next_network();
	if (! $oldnet) {
	return 0;
      }
      $self->{oldnet}=$oldnet;
      my @candidates=$oldnet->nodes();
      $self->{candidates}=\@candidates;
      my @candidatesbis=$oldnet->internal_nodes();
      $self->{candidatesbis}=\@candidatesbis;
      $self->{processbis}=0;
      $self->{index}=0;
      my $n1=scalar @candidates;
      my $n2=scalar @candidatesbis;
      print $oldnet->eNewick()."(".$self->{numleaves}.")($n1,$n2):\n";
    }
    if (($self->{index} == (scalar @{$self->{candidates}})) && $self->{processbis}==0
       ) {
      $self->{processbis}=1;
      $self->{index}=0;
      print "--\n";
    }
    if ($self->{processbis}==0) {
      my $graph=$self->{oldnet}->{graph}->copy();
      my $u=$self->{candidates}->[$self->{index}];
      foreach my $w ($graph->predecessors($u)) {
	$graph->delete_edge($w,$u);
	$graph->add_edge($w,"t$n");
      }
      $graph->add_edge("t$n",$u);
      $graph->add_edge("t$n",$self->{newleaf});
      my $net=Bio::PhyloNetwork->new(-graph=>$graph);
      $self->{index}++;
      return $net;
    } else {
      my $graph=$self->{oldnet}->{graph}->copy();
      my $u=$self->{candidatesbis}->[$self->{index}];
#      print "<<$u\n";
      $graph->add_edge($u,$self->{newleaf});
      my $net=Bio::PhyloNetwork->new(-graph=>$graph);
      $self->{index}++;
      return $net;
    }
  }
}

1;
