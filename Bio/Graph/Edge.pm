# $Id$
#
# BioPerl module for Bio::Graph::Edge
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Graph::Edge - encapsulation of an interaction between 2 Bio::Seq objects

=head1 SYNOPSIS

  ## get an interaction between two nodes ##

  my $edge  = $gr->edge( $gr->nodes_by_id('P12345'),
                         $gr->nodes_by_id('P23456'));
  my $id    = $edge->object_id();
  my $wt    = $edge->weight();
  my @nodes = $edge->nodes();

=head1 DESCRIPTION

This class contains information about a bimolecular interaction.
At present it just contains data about its component node, a weight
(if set) and an identifier. Subclasses could hold more specific 
information.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Richard Adams

Email richard.adams@ed.ac.uk

=cut

use strict;
package Bio::Graph::Edge;
use base qw(Bio::Root::Root Bio::IdentifiableI);


=head2 new

 Name       : new
 Purpose    : constructor for an edge object
 Usage      : my $edge = Bio::Graph::Edge->new(nodes => [$node1,$node2]
                                               id => $id);
              $graph->add_edge($edge);
 Returns    : a new Bio::Graph::Edge object 
 Arguments  : hash nodes            => array reference of 2 nodes
                   id               => edge id
                   weight(optional) => weight score.

=cut

sub new {
      ##array based, not hash based ##..., therefore does not use 
      #Bio::Root::Root->new().

	my ($caller, @args) = @_;
	my $class  = ref ($caller) || $caller;
	my $self    = [];
	bless ($self, $class);

	my ($weight, $id, $nodes) = $self->_rearrange([qw( WEIGHT ID NODES)], @args);
	$self->[0] = $nodes->[0];
	$self->[1] = $nodes->[1];
	$self->[2] = defined($weight)?$weight:undef; 
	$self->[3] = defined($id)?$id:undef; 
	return $self;

}

=head2 weight

 Name      : weight
 Purpose   : get/setter for weight score
 Usage     : my $weight = $edge->weight();
 Returns   : anumber
 Arguments : void/ a number

=cut

sub weight {
	my $self = shift;
	if (@_) {$self->[2] = shift;}
	return defined($self->[2])?$self->[2]:undef;
}

=head2 object_id

 Name      : object_id
 Purpose   : get/setter for object_id
 Usage     : my $id = $edge->object_id();
 Returns   : a string identifier
 Arguments : void/ an identifier 

=cut

sub object_id {
	my $self            = shift;
	if (@_) {
		my $v  = shift;
		if (ref ($v)) {
			$self->throw ("Edge ID must be a text value, not a [".
							ref($v). "].");
			} 
		$self->[3] = shift;
	}
	return defined($self->[3])?$self->[3]:undef;
}

=head2  nodes

 Name      : nodes
 Purpose   : get/setter for nodes
 Usage     : my @nodes = $edge->nodes();
 Returns   : a 2 element list of nodes /void
 Arguments : void/ a 2 element list of nodes. 

=cut

sub nodes {
	my ($self, @args) = @_;
	if (@args >= 2 ) {
		$self->[0] =  $args[0];
		$self->[1] =  $args[1];
		}
	return ($self->[0], $self->[1]);
}

1;
