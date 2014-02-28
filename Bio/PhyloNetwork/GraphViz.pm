#
# Module for Bio::PhyloNetwork::GraphViz
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

Bio::PhyloNetwork::GraphViz - Interface between PhyloNetwork and GraphViz

=head1 SYNOPSIS

  use Bio::PhyloNetwork;
  use Bio::PhyloNetwork::GraphViz;

  my $net=Bio::PhyloNetwork->new(
      -eNewick=>'((H1,(H1,(H2,l))),H2)t0; (some long label)H1; ("quoted label")H2;'
  );

  my $gv=Bio::PhyloNetwork::GraphViz->new(-net=>$net,-short_labels=>1);

  foreach my $u ($net->nodes()) {
      print "$u:".$gv->nodePN_to_nodeGV->{$u}."\n";
  }

  print $gv->as_text;

  open my $PS, '>', "net.ps" or die "Could not write file 'net.ps': $!\n";
  print $PS $gv->as_ps;
  close $PS;

=head1 DESCRIPTION

This is a module to create GraphViz objects representing phylogenetic networks.

=head1 AUTHOR

Gabriel Cardona, gabriel(dot)cardona(at)uib(dot)es

=head1 SEE ALSO

L<Bio::PhyloNetwork>, L<GraphViz>

=head1 APPENDIX

The rest of the documentation details each of the object methods.

=cut

package Bio::PhyloNetwork::GraphViz;

use strict;
use warnings;

use base qw(Bio::Root::Root GraphViz);

use Bio::PhyloNetwork;

=head2 new

 Title   : new
 Usage   : my $graphv = new Bio::PhyloNetwork::GraphViz();
 Function: Creates a new Bio::PhyloNetwork::GraphViz object
 Returns : Bio::PhyloNetwork::GraphViz
 Args    : -net => Bio::PhyloNetwork object
           -short_labels => boolean (optional)

Returns a Bio::PhyloNetwork::GraphViz object, which is an extension of
a GraphViz object. The GraphViz object is a representation of the
phylogenetic network given. The extra information the created object
holds is a hash table with keys the nodes of the PhyloNetwork object
and values the nodes of the GraphViz object. If the optional argument
-short_labels=E<gt>1 is given, the labels of the nodes in GraphViz are
shortened to a maximum of 3 letters.

=cut

sub new {
  my ($pkg,@args)=@_;

  my $self=$pkg->SUPER::new(@args);

  my ($net,$short_labels)=
    $self->_rearrange([qw(NET
			  SHORT_LABELS)],@args);
  if (! defined $short_labels) {
    $short_labels=0;
  }
  my $gv=GraphViz->new();
  my $nodePN_to_nodeGV={};
  my @nodes=$net->nodes();
  foreach my $node (@nodes) {
#    my $namenode=generate_name($node);
#    $names->{$node}=$namenode;
    ###
    my $labelnodeint=$net->{labels}->{$node};
    ###
    my $labelnode=($short_labels ? find_short_label($labelnodeint) : find_label($labelnodeint));
    my $nodeGV=
      $gv->add_node(#$namenode,
		    label=>$labelnode,
		    shape=>($net->is_tree_node($node) ? 'circle' : 'box'));
    $nodePN_to_nodeGV->{$node}=$nodeGV;
  }
  my @edges=$net->edges();
  foreach my $edge (@edges) {
    my $node1=$edge->[0];
#    my $namenode1=generate_name($node1);
    my $node2=$edge->[1];
#    my $namenode2=generate_name($node2);
    $gv->add_edge($nodePN_to_nodeGV->{$node1},$nodePN_to_nodeGV->{$node2});
  }
  $self=$gv;
  $self->{nodePN_to_nodeGV}=$nodePN_to_nodeGV;
  bless($self,$pkg);
}

#sub generate_name {
#  my ($var)=@_;
#  if ($var =~ /\D/) {
#    print "$var contains a number.\b";
#    return $var;
#  }
#  return "N$var";
#}

sub find_short_label {
  my ($str)=@_;
  return substr(find_label($str),0,3);
}

sub find_label {
  my ($str)=@_;
  $str =~ tr/A-Za-z0-9//cd;
  return $str;
}

=head2 nodePN_to_nodeGV

 Title   : nodePN_to_nodeGV
 Usage   : my $hashR=$graphv->nodePN_to_nodeGV()
 Function: returns (a reference to) a hash holding the translation between
           nodes of the Bio::PhyloNetwork object and nodes of the GraphViz
           object
 Returns : reference to hash
 Args    : none

=cut

sub nodePN_to_nodeGV {
  my ($self)=@_;
  return $self->{nodePN_to_nodeGV};
}

1;
