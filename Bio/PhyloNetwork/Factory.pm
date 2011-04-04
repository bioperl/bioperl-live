#
# Module for Bio::PhyloNetwork::Factory
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

Bio::PhyloNetwork::Factory - Module to sequentially generate
Phylogenetic Networks

=head1 SYNOPSIS

 use strict;
 use warnings;

 use Bio::PhyloNetwork;
 use Bio::PhyloNetwork::Factory;

 # Will generate sequentially all the 4059 binary tree-child phylogenetic
 # networks with 4 leaves

 my $factory=Bio::PhyloNetwork::Factory->new(-numleaves=>4);

 my @nets;

 while (my $net=$factory->next_network()) {
   push @nets,$net;
   print "".(scalar @nets).": ".$net->eNewick()."\n";
 }

=head1 DESCRIPTION

Sequentially builds a (binary tree-child) phylogenetic network each time
next_network is called.

=head1 AUTHOR

Gabriel Cardona, gabriel(dot)cardona(at)uib(dot)es

=head1 SEE ALSO

L<Bio::PhyloNetwork>

=head1 APPENDIX

The rest of the documentation details each of the object methods.

=cut

package Bio::PhyloNetwork::Factory;

use strict;
use warnings;

use base qw(Bio::Root::Root);

use Bio::PhyloNetwork;
use Bio::PhyloNetwork::TreeFactory;

=head2 new

 Title   : new
 Usage   : my $factory = new Bio::PhyloNetwork::Factory();
 Function: Creates a new Bio::PhyloNetwork::Factory
 Returns : Bio::PhyloNetwork::RandomFactory
 Args    : -numleaves => integer
            OR
           -leaves => reference to an array (of leaves names)
           -numhybrids => integer [default = numleaves -1]
           -recurse => boolean [optional]

Returns a Bio::PhyloNetwork::Factory object. Such an object will
sequentially create binary tree-child phylogenetic networks
each time next_network is called.

If the parameter -leaves=E<gt>\@leaves is given, then the set of leaves of
these networks will be @leaves. If it is given the parameter
-numleaves=E<gt>$numleaves, then the set of leaves will be "l1"..."l$numleaves".

If the parameter -numhybrids=E<gt>$numhybrids is given, then the generated
networks will have at most $numhybrids hybrid nodes. Note that, necessarily,
$numhybrids E<lt> $numleaves.

If the parameter -recurse=E<gt>1 is given, then all networks with number of hybrid
nodes less or equal to $numhybrids will be given; otherwise only those with
exactly $numhybrids hybrid nodes.

=cut

sub new {
  my ($pkg,@args)=@_;

  my $self=$pkg->SUPER::new(@args);

  my ($leavesR,$numleaves,$numhybrids,$recurse)=
    $self->_rearrange([qw(LEAVES
			  NUMLEAVES
			  NUMHYBRIDS
			  RECURSE)],@args);
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

  $recurse ||= 0;
  if (! defined $numhybrids) {
    $numhybrids=$numleaves-1;
    $recurse=1;
  }
  $self->{recurse}=$recurse;
  $self->{numhybrids}=$numhybrids;
  if ($numhybrids ==0) {
    return Bio::PhyloNetwork::TreeFactory->new(-leaves=>\@leaves);
  }
  my $parent;
  if ($numhybrids > 1) {
    $parent=new($pkg,'-leaves'=>\@leaves,
		'-numhybrids'=>($numhybrids-1),
		'-recurse'=>($recurse));
  }
  else {
    $parent=Bio::PhyloNetwork::TreeFactory->new(-leaves=>\@leaves);
  }
  $self->{parent}=$parent;
  my $oldnet=$parent->next_network();
  $self->{oldnet}=$oldnet;
  $self->update();
  $self->{found}=[];
  bless($self,$pkg);
}

sub update {
  my ($self)=@_;

  my @candidates=$self->{oldnet}->edges();
  $self->{candidates}=\@candidates;
  $self->{numcandidates}=(scalar @candidates);
  $self->{index1}=-$self->{recurse};
  $self->{index2}=0;
}

=head2 next_network

 Title   : next_network
 Usage   : my $net=$factory->next_network()
 Function: returns a network
 Returns : Bio::PhyloNetwork
 Args    : none

=cut

sub next_network {
  my ($self)=@_;
  my $numleaves=$self->{numleaves};
  my $numhybrids=$self->{numhybrids};
  START:
  if ($self->{index1}==-1) {
    $self->{index1}++;
    return $self->{oldnet};
  }
  if ($self->{index1} >= $self->{numcandidates}) {
    $self->{index2}++;
    $self->{index1}=0;
  }
  if ($self->{index2} >= $self->{numcandidates}) {
    my $oldnet=$self->{parent}->next_network();
    if (! $oldnet) {
      return 0;
    }
    $self->{oldnet}=$oldnet;
    $self->update();
    goto START;
  }
  if ((scalar $self->{oldnet}->hybrid_nodes())< $self->{numhybrids}-1) {
    $self->{candidates}=[];
    $self->{numcandidates}=0;
    goto START;
  }
  my $u1=$self->{candidates}->[$self->{index1}]->[0];
  my $v1=$self->{candidates}->[$self->{index1}]->[1];
  my $u2=$self->{candidates}->[$self->{index2}]->[0];
  my $v2=$self->{candidates}->[$self->{index2}]->[1];
  my $lbl=$self->{numhybrids};
  if ($self->{oldnet}->is_attackable($u1,$v1,$u2,$v2)) {
    my $net=Bio::PhyloNetwork->new(-graph=>$self->{oldnet}->graph);
    $net->do_attack($u1,$v1,$u2,$v2,$lbl);
    $self->{index1}++;
    my @found=@{$self->{found}};
    foreach my $netant (@found) {
      if ($net->is_mu_isomorphic($netant) ) {
	goto START;
      }
    }
    push @found,$net;
    $self->{found}=\@found;
    return $net;
  }
  else {
    $self->{index1}++;
    goto START;
  }
}

1;
