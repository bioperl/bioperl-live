# $Id$
#
# BioPerl module for OntologyEngineI
#
# Cared for by Peter Dimitrov <dimitrov@gnf.org>
#
# (c) Peter Dimitrov
# (c) GNF, Genomics Institute of the Novartis Research Foundation, 2002.
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

OntologyEngineI - DESCRIPTION of Interface

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the interface here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Peter Dimitrov

Email dimitrov@gnf.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Ontology::OntologyEngineI;
use vars qw(@ISA);
use strict;
use Carp;
use Bio::Root::RootI;

@ISA = qw( Bio::Root::Root );

=head2 add_term

 Title   : add_term
 Usage   : add_term(TermI term): TermI
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_term{
  my ($self) = @_;

  $self->throw("Abstract method add_term implementing class did not provide method");
}

=head2 add_relationship

 Title   : add_relationship
 Usage   : add_relationship(RelationshipI relationship)
  add_relatioship(TermI parent, TermI child, TermI relationship_type)
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_relationship{
  my ($self, $relationship) = @_;

  $self->throw("Abstract method add_relationship implementing class did not provide method");

}

=head2 get_relationships

 Title   : get_relationships
 Usage   : get_relationships([TermI term]): RelationshipI[]
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_relationships{
  my ($self) = @_;

  $self->throw("Abstract method get_relationships implementing class did not provide method");
}

=head2 get_relationship_types

 Title   : get_relationship_types
 Usage   : get_relationship_types(): TermI[]
 Function:
 Example :
 Returns :
 Args    :


=cut

sub get_relationship_types{
  my ($self) = @_;

  $self->throw("Abstract method get_relationship_types implementing class did not provide method");
}

=head2 get_child_terms

 Title   : get_child_terms
 Usage   : get_child_terms(TermI term, TermI[] rel_types): TermI[]
 Function:
 Example :

 Returns :
 Args    :


=cut

sub get_child_terms{
  my ($self) = @_;

  $self->throw("Abstract method get_child_terms implementing class did not provide method");
}

=head2 get_descendant_terms

 Title   : get_descendant_terms
 Usage   : get_descendant_terms(TermI term, TermI[] rel_types): TermI[]
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_descendant_terms{
  my ($self) = @_;

  $self->throw("Abstract method get_descendant_terms implementing class did not provide method");
}

=head2 get_parent_terms

 Title   : get_parent_terms
 Usage   : get_parent_terms(TermI term, TermI[] rel_types): TermI[]
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_parent_terms{
  my ($self) = @_;

  $self->throw("Abstract method get_parent_terms implementing class did not provide method");
}

=head2 get_ancestor_terms

 Title   : get_ancestor_terms
 Usage   : get_ancestor_terms(TermI term, TermI[] rel_types): TermI[]
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_ancestor_terms{
  my ($self) = @_;

  $self->throw("Abstract method get_ancestor_terms implementing class did not provide method");
}

=head2 get_leaf_terms

 Title   : get_leaf_terms
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_leaf_terms{
  my ($self) = @_;

  $self->throw("Abstract method get_leaf_terms implementing class did not provide method");
}

=head2 get_root_terms()

 Title   : get_root_terms()
 Usage   : get_root_terms()
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_root_terms(){
  my ($self) = @_;

  $self->throw("Abstract method get_root_terms() implementing class did not provide method");
}

1;
