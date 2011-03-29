# $Id: GraphAdaptor02.pm 10525 2006-09-26 22:03:22Z sendu $
#
# BioPerl adaptor for old Graph verions (0.2x) for use in 
# Bio::Ontology::SimpleGOEngine
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

Bio::Ontology::SimpleGOEngine::GraphAdaptor02 - Graph adaptor (v02.x) for
Bio::Ontology::SimpleGOEngine

=head1 DESCRIPTION

Internal subclass of Bio::Ontology::SimpleGOEngine::GraphAdaptor for
Graph v0.2x.

Call this via Bio::Ontology::SimpleGOEngine::GraphAdaptor

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

  https://redmine.open-bio.org/projects/bioperl/

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

package Bio::Ontology::SimpleGOEngine::GraphAdaptor02;


use strict;

use base qw(Bio::Ontology::SimpleGOEngine::GraphAdaptor);

# edges
#   v0.2x returns (u0,v0, u1,v1, ...)
#   v0.5x returns ([u0,v0], [u1,v1], ...)
sub edges {
  my $self=shift;
  my @edges02=$self->_graph->edges(@_);
  my @edges;
  while (@edges02) {
    my($u,$v)=(shift @edges02,shift @edges02);
    push(@edges,[$u,$v]);
  }
  @edges;
}

# edges_at
#   v0.2x uses edges() method and returns (u0,v0, u1,v1, ...)
#   v0.5x returns ([u0,v0], [u1,v1], ...)
sub edges_at {
  my $self=shift;
  $self->edges(@_);
}

# set_vertex_attribute
#   v0.2x uses set_attribute($attribute,$v,$value)
sub set_vertex_attribute {
  my($self,$v,$attribute,$value)=@_;
  $self->_graph->set_attribute($attribute,$v,$value);
}

# get_vertex_attribute
#   v0.2x uses get_attribute($attribute,$v)
sub get_vertex_attribute {
  my($self,$v,$attribute)=@_;
  $self->_graph->get_attribute($attribute,$v);
}

# set_edge_attribute
#   v0.2x uses set_attribute($attribute,$u,$v,$value)
sub set_edge_attribute {
  my($self,$u,$v,$attribute,$value)=@_;
  $self->_graph->set_attribute($attribute,$u,$v,$value);
}

# get_edge_attribute
#   v0.2x uses get_attribute($attribute,$u,$v)
sub get_edge_attribute {
  my($self,$u,$v,$attribute)=@_;
  $self->_graph->get_attribute($attribute,$u,$v);
}

1;
