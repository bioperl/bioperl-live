# $Id: GraphAdaptor.pm 10525 2006-09-26 22:03:22Z sendu $
#
# BioPerl Graph adaptor for Bio::Ontology::SimpleGOEngine
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Nat Goodman <natg at shore.net>
#
# (c) Nathan Goodman natg@shore.net 2005
# (c) ISB, Institute for Systems Biology 2005
#
# You may distribute this module under the same terms as perl itself.
# Refer to the Perl Artistic License (see the license accompanying this
# software package, or see http://www.perl.com/language/misc/Artistic.html)
# for the terms under which you may use, modify, and redistribute this module.
#
# THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
# MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Ontology::SimpleGOEngine::GraphAdaptor - Graph adaptor for
Bio::Ontology::SimpleGOEngine

=head1 SYNOPSIS

  use Bio::Ontology::SimpleGOEngine::GraphAdaptor;

  my $graph = Bio::Ontology::SimpleGOEngine::GraphAdaptor;

=head1 DESCRIPTION

This is an adaptor to simplify use of versions of the standard CPAN Graph module
(old is versions 0.2x; new is 0.5x and beyond) within
Bio::Ontology::SimpleGOEngine. Prior versions of this module supported Graph
version older than 0.5, however we are removing support for these older version
post BioPerl 1.6.901. If you absolutely require an old version of Graph, please
use an older version of BioPerl.

This module implements only those Graph methods used by SimpleGOEngine. It is
far from a complete compatibility layer! It also implements workarounds for
certain performance problems in the current versions of Graph v0.5x.

This class provides implementations for the required graph methods using the new
version of Graph. In most cases, these are simple pass-throughs

The methods implemented here or in the subclasses are listed below.
In all cases, we implemented the Graph v0.5x interface.  Consult the
Graph v0.5x man page for details.

  add_vertex
  has_vertex
  add_edge
  has_edge
  vertices
  edges
  edges_at
  predecessors
  successors
  set_vertex_attribute
  get_vertex_attribute
  set_edge_attribute
  get_edge_attribute
  source_vertices
  sink_vertices

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

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

report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Nat Goodman

Email: natg at shore.net

Address:

  Institute for Systems Biology
  1441 N 34th St
  Seattle, WA 98103-8904

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Ontology::SimpleGOEngine::GraphAdaptor;

use Graph::Directed;

use strict;

use base qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : $graph = Bio::Ontology::SimpleGOEngine::GraphAdaptor->new()
 Function: Creates a new graph
 Returns : Bio::Ontology::SimpleGOEngine::GraphAdaptor02 or
           Bio::Ontology::SimpleGOEngine::GraphAdaptor05 object,
           depending on which Graph version is available
 Args    : none

=cut

sub new {
  my( $class ) = @_;
  $class = ref $class || $class;

  my $self= bless( {}, $class );
  $self->{_graph}=Graph::Directed->new();
  $self->{_vertex_attributes}={};
  $self->{_edge_attributes}={};
  return $self;
}

# Here are the main methods

sub add_vertex {
  my $self=shift;
  $self->_graph->add_vertex(@_);
}
sub has_vertex {
  my $self=shift;
  $self->_graph->has_vertex(@_);
}
sub add_edge {
  my $self=shift;
  $self->_graph->add_edge(@_);
}
sub has_edge {
  my $self=shift;
  $self->_graph->has_edge(@_);
}
sub vertices {
  my $self=shift;
  $self->_graph->vertices(@_);
}
sub edges {
  my $self=shift;
  $self->_graph->edges(@_);
}
sub edges_at {
  my $self=shift;
  $self->_graph->edges_at(@_);
}
sub predecessors {
  my $self=shift;
  $self->_graph->predecessors(@_);
}
sub successors {
  my $self=shift;
  $self->_graph->successors(@_);
}
sub source_vertices {
  my $self=shift;
  $self->_graph->source_vertices();
}
sub sink_vertices {
  my $self=shift;
  $self->_graph->sink_vertices();
}
# The following methods workaround a performance problem in Graph v0.5x
# when attributes are attached to the graph
sub set_vertex_attribute {
  my($self,$v,$attribute,$value)=@_;
  $self->_vertex2attributes($v)->{$attribute}=$value;
}
sub get_vertex_attribute {
  my($self,$v,$attribute)=@_;
  $self->_vertex2attributes($v)->{$attribute};
}
sub set_edge_attribute {
  my($self,$u,$v,$attribute,$value)=@_;
  $self->_edge2attributes($u,$v)->{$attribute}=$value;
}
sub get_edge_attribute {
  my($self,$u,$v,$attribute)=@_;
  $self->_edge2attributes($u,$v)->{$attribute};
}

=head2 _graph

 Title   : _graph
 Usage   : $self->_graph();
 Function: Internal method to access 'real' graph
 Returns : Graph::Directed object
 Args    : none

=cut

sub _graph {$_[0]->{_graph}; }

=head2 _vertex_attributes

 Title   : _vertex_attributes
 Usage   : $self->vertex_attributes();
 Function: Internal method to access HASH used to store vertex attributes
 Returns : Graph::Directed object
 Args    : none

=cut

sub _vertex_attributes {$_[0]->{_vertex_attributes}; }

=head2 _edge_attributes

 Title   : _edge_attributes
 Usage   : $self->edge_attributes();
 Function: Internal method to access HASH used to store edge attributes
 Returns : Graph::Directed object
 Args    : none

=cut

sub _edge_attributes {$_[0]->{_edge_attributes}; }

=head2 _vertex2attributes

 Title   : _vertex2attributes
 Usage   : $value=$graph->_vertex2attributes($v_->{ATTRIBUTE};
           $graph->_vertex2attributes($v)->{ATTRIBUTE}=$value;
 Function: Internal method to access attributes for a specific vertex
 Returns : HASH
 Args    : none

=cut

sub _vertex2attributes {
  my($self,$vertex)=@_;
  $self->_vertex_attributes->{$vertex} or $self->_vertex_attributes->{$vertex}={};
}

=head2 _edge2attributes

 Title   : _edge2attributes
 Usage   : $value=$graph->_edge2attributes($u,$v)->{ATTRIBUTE};
           $graph->_edge2attributes($u,$v)->{ATTRIBUTE}=$value;
 Function: Internal method to access HASH used to store edge attributes
 Returns : HASH
 Args    : none

=cut

sub _edge2attributes {
  my($self,$u,$v)=@_;
  $self->_edge_attributes->{$u,$v} or $self->_edge_attributes->{$u,$v}={};
}

1;
