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

OntologyEngineI - Interface a minimal Ontology implementation should satisfy

=head1 SYNOPSIS

    # see documentation of methods

=head1 DESCRIPTION

This describes the minimal interface an Ontology implementation should
provide.  It intentionally doesn't make explicit references to the
ontology being a DAG, nor does it mandate that the ontology be a
vocabulary. Rather, it tries to generally express what should be
accessible about an ontology.

The idea is to allow for different implementations for different
purposes, which may then differ as to which operations are efficient
and which aren't, and how much richer the functionality is on top of
this minimalistic set of methods. Check modules in the Bio::Ontology
namespace to find out which implementations exist. At the time of
writing, there is a SimpleOntologyEngine (which does not use
Graph.pm), and a Graph.pm-based implementation in SimpleGOEngine.

Ontology parsers in Bio::OntologyIO are required to return an
implementation of this interface.

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
  http://bugzilla.bioperl.org/

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
 Function: Adds TermI object to the ontology engine term store
 Example : $oe->add_term($term)
 Returns : its argument.
 Args    : object of class TermI.


=cut

sub add_term{
    shift->throw_not_implemented();
}

=head2 add_relationship

 Title   : add_relationship
 Usage   : add_relationship(RelationshipI relationship): RelationshipI
  add_relatioship(TermI parent, TermI child, TermI relationship_type)
 Function: Adds a relationship object to the ontology engine.
 Example :
 Returns : Its argument.
 Args    : A RelationshipI object.


=cut

sub add_relationship{
    shift->throw_not_implemented();
}

=head2 get_relationships

 Title   : get_relationships
 Usage   : get_relationships(): RelationshipI[]
 Function: Retrieves all relationship objects.
 Example :
 Returns : Array of RelationshipI objects
 Args    :


=cut

sub get_relationships{
    shift->throw_not_implemented();
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
    shift->throw_not_implemented();
}

=head2 get_child_terms

 Title   : get_child_terms
 Usage   : get_child_terms(TermI term, TermI[] relationship_types): TermI[]
 Function: Retrieves all child terms of a given term, that satisfy a
           relationship among those that are specified in the second
           argument or undef otherwise. get_child_terms is a special
           case of get_descendant_terms, limiting the search to the
           direct descendants.

 Example :
 Returns : Array of TermI objects.
 Args    : First argument is the term of interest, second is the list
           of relationship type terms.


=cut

sub get_child_terms{
    shift->throw_not_implemented();
}

=head2 get_descendant_terms

 Title   : get_descendant_terms
 Usage   : get_descendant_terms(TermI term, TermI[] rel_types): TermI[]
 Function: Retrieves all descendant terms of a given term, that
           satisfy a relationship among those that are specified in
           the second argument or undef otherwise. 
 Example :
 Returns : Array of TermI objects.
 Args    : First argument is the term of interest, second is the list
           of relationship type terms.


=cut

sub get_descendant_terms{
    shift->throw_not_implemented();
}

=head2 get_parent_terms

 Title   : get_parent_terms
 Usage   : get_parent_terms(TermI term, TermI[] relationship_types): TermI[]
 Function: Retrieves all parent terms of a given term, that satisfy a
           relationship among those that are specified in the second
           argument or undef otherwise. get_parent_terms is a special
           case of get_ancestor_terms, limiting the search to the
           direct ancestors.

 Example :
 Returns : Array of TermI objects.
 Args    : First argument is the term of interest, second is the list
           of relationship type terms.


=cut

sub get_parent_terms{
    shift->throw_not_implemented();
}

=head2 get_ancestor_terms

 Title   : get_ancestor_terms
 Usage   : get_ancestor_terms(TermI term, TermI[] relationship_types): TermI[]
 Function: Retrieves all ancestor terms of a given term, that satisfy
           a relationship among those that are specified in the second
           argument or undef otherwise. 

 Example :
 Returns : Array of TermI objects.
 Args    : First argument is the term of interest, second is the list
           of relationship type terms.


=cut

sub get_ancestor_terms{
    shift->throw_not_implemented();
}

=head2 get_leaf_terms

 Title   : get_leaf_terms
 Usage   : get_leaf_terms(): TermI[]
 Function: Retrieves all leaf terms from the ontology. Leaf term is a
           term w/o descendants.

 Example : @leaf_terms = $obj->get_leaf_terms()
 Returns : Array of TermI objects.
 Args    :


=cut

sub get_leaf_terms{
    shift->throw_not_implemented();
}

=head2 get_root_terms()

 Title   : get_root_terms
 Usage   : get_root_terms(): TermI[]
 Function: Retrieves all root terms from the ontology. Root term is a
           term w/o descendants.

 Example : @root_terms = $obj->get_root_terms()
 Returns : Array of TermI objects.
 Args    :


=cut

sub get_root_terms(){
    shift->throw_not_implemented();
}

1;
