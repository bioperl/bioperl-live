#
# Module for Bio::PhyloNetwork::RandomFactory
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

Bio::PhyloNetwork::RandomFactory - Module to generate random
Phylogenetic Networks

=head1 SYNOPSIS

 use strict;
 use warnings;

 use Bio::PhyloNetwork;
 use Bio::PhyloNetwork::RandomFactory;

 # Will generate at random all the 66 binary tree-child phylogenetic
 # networks with 3 leaves

 my $factory=Bio::PhyloNetwork::RandomFactory->new(-numleaves=>3,-norepeat=>1);

 my @nets;

 for (my $i=0; $i<66; $i++) {
   my $net=$factory->next_network();
   push @nets,$net;
   print "".(scalar @nets).": ".$net->eNewick()."\n";
 }

=head1 DESCRIPTION

Builds a random (binary tree-child) phylogenetic network each time
next_network is called.

=head1 AUTHOR

Gabriel Cardona, gabriel(dot)cardona(at)uib(dot)es

=head1 SEE ALSO

L<Bio::PhyloNetwork>

=head1 APPENDIX

The rest of the documentation details each of the object methods.

=cut

package Bio::PhyloNetwork::RandomFactory;

use strict;
use warnings;

use base qw(Bio::Root::Root);

use Bio::PhyloNetwork;
use Bio::Tree::RandomFactory;

=head2 new

 Title   : new
 Usage   : my $factory = new Bio::PhyloNetwork::RandomFactory();
 Function: Creates a new Bio::PhyloNetwork::RandomFactory
 Returns : Bio::PhyloNetwork::RandomFactory
 Args    : -numleaves => integer
            OR
           -leaves => reference to an array (of leaves names)
           -numhybrids => integer [optional]
           -norepeat => boolean [optional]

Returns a Bio::PhyloNetwork::RandomFactory object. Such an object will create
random binary tree-child phylogenetic networks each time next_network
is called.

If the parameter -leaves=E<gt>\@leaves is given, then the set of leaves of
these networks will be @leaves. If it is given the parameter
-numleaves=E<gt>$numleaves, then the set of leaves will be "l1"..."l$numleaves".

If the parameter -numhybrids=E<gt>$numhybrids is given, then the generated
networks will have exactly $numhybrids hybrid nodes. Note that, necessarily,
$numhybrids E<lt> $numleaves. Otherwise, the number of hybrid nodes will be chosen
at random for each call of next_network.

If the parameter -norepeat=E<gt>1 is given, then successive calls of next_network
will give non-isomorphic networks.

=cut

sub new {
  my ($pkg,@args)=@_;

  my $self=$pkg->SUPER::new(@args);

  my ($leavesR,$numleaves,$numhybrids,$norepeat)=
    $self->_rearrange([qw(LEAVES
			  NUMLEAVES
			  NUMHYBRIDS
			  NOREPEAT)],@args);
  my @leaves;
  if ((! defined $leavesR) && (defined $numleaves)) {
    @leaves=map {"l$_"} (1..$numleaves);
    $leavesR=\@leaves;
  }
  if (! defined $leavesR) {
    $self->throw("No leaves set neither numleaves given");
  }
  $norepeat ||= 0;

  $self->{leaves}=\@leaves;
  $self->{numleaves}=$numleaves;
  $self->{numhybrids}=$numhybrids if defined $numhybrids;
  $self->{norepeat}=$norepeat;
  $self->{found}=[];
  $self->{tree_factory}=Bio::Tree::RandomFactory->new(-taxa => \@leaves);
  bless($self,$pkg);
}

=head2 next_network

 Title   : next_network
 Usage   : my $net=$factory->next_network()
 Function: returns a random network
 Returns : Bio::PhyloNetwork
 Args    : none

=cut

sub next_network {
  my ($self)=@_;

  my $numleaves=$self->{numleaves};
  my @found=@{$self->{found}};
  my $numhybrids;
 START:
  if (! defined $self->{numhybrids}) {
    $numhybrids=int(rand($numleaves));
  }
  else {
    $numhybrids=$self->{numhybrids};
  }
  my $tf=$self->{tree_factory};
  my $tree=$tf->next_tree;
  my $net=Bio::PhyloNetwork->new(-tree=>$tree);
  for (my $i=1; $i<=$numhybrids; $i++) {
    $net=random_attack($net,$i);
  }
  if ($self->{norepeat}) {
    foreach my $ant (@found) {
      goto START if $net->is_mu_isomorphic($ant);
    }
    push @found,$net;
    $self->{found}=\@found;
  }
  return $net;
}

sub random_attack {
  my ($net,$lbl)=@_;

  my $graph=$net->{graph};
  my ($u1,$v1,$u2,$v2);
  do {
    my $e1=$graph->random_edge;
    my $e2=$graph->random_edge;
    $u1=$e1->[0];
    $v1=$e1->[1];
    $u2=$e2->[0];
    $v2=$e2->[1];
  } while (! $net->is_attackable($u1,$v1,$u2,$v2,$lbl));
  $net->do_attack($u1,$v1,$u2,$v2,$lbl);
  return $net;
}

1;
