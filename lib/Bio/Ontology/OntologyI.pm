#
# BioPerl module for Bio::Ontology::OntologyI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Hilmar Lapp <hlapp at gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

#
# (c) Hilmar Lapp, hlapp at gmx.net, 2003.
# (c) GNF, Genomics Institute of the Novartis Research Foundation, 2003.
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

# POD documentation - main docs before the code

=head1 NAME

Bio::Ontology::OntologyI - Interface for an ontology implementation

=head1 SYNOPSIS

    # see method documentation

=head1 DESCRIPTION

This describes the minimal interface an ontology implementation must
provide. In essence, it represents a namespace with description on top
of the query interface OntologyEngineI.

This interface inherits from L<Bio::Ontology::OntologyEngineI>.

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
email or the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Hilmar Lapp

Email hlapp at gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Ontology::OntologyI;
use strict;


use base qw(Bio::Ontology::OntologyEngineI);

=head1  Methods defined in this interface.

=cut

=head2 name

 Title   : name
 Usage   : $obj->name($newval)
 Function: Get/set the name of this ontology.
 Example : 
 Returns : value of name (a scalar)
 Args    : 

=cut

sub name{
    shift->throw_not_implemented();
}

=head2 authority

 Title   : authority
 Usage   : $auth = $obj->authority()
 Function: Get/set the authority for this ontology, for instance the
           DNS base for the organization granting the name of the
           ontology and identifiers for the terms.

           This attribute is optional and should not generally
           expected by applications to have been set. It is here to
           follow the rules for namespaces, which ontologies serve as
           for terms.

 Example : 
 Returns : value of authority (a scalar)
 Args    : 

=cut

sub authority{
    shift->throw_not_implemented();
}

=head2 identifier

 Title   : identifier
 Usage   : $id = $obj->identifier()
 Function: Get an identifier for this ontology.

           This is primarily intended for look-up purposes. Clients
           should not expect the value to be modifiable, and it may
           not be allowed to set its value from outside. Also, the
           identifier's uniqueness may only hold within the scope of a
           particular application's run time, i.e., it may be a memory
           location.

 Example : 
 Returns : value of identifier (a scalar)
 Args    : 

=cut

sub identifier{
    shift->throw_not_implemented();
}

=head2 definition

 Title   : definition
 Usage   : $def = $obj->definition()
 Function: Get a descriptive definition for this ontology.
 Example : 
 Returns : value of definition (a scalar)
 Args    : 

=cut

sub definition{
    shift->throw_not_implemented();
}

=head2 close

 Title   : close
 Usage   :
 Function: Release any resources this ontology may occupy. In order
           to efficiently release used memory or file handles, you
           should call this method once you are finished with an
           ontology.

 Example :
 Returns : TRUE on success and FALSE otherwise
 Args    : none

=cut

sub close{
    shift->throw_not_implemented();
}

=head1 Methods inherited from L<Bio::Ontology::OntologyEngineI>

Their documentations are copied here for completeness. In most use
cases, you will want to access the query methods of an ontology, not
just the name and description ...

=cut

=head2 add_term

 Title   : add_term
 Usage   : add_term(TermI term): TermI
 Function: Adds TermI object to the ontology engine term store.

           For ease of use, if the ontology property of the term
           object was not set, an implementation is encouraged to set
           it to itself upon adding the term.

 Example : $oe->add_term($term)
 Returns : its argument.
 Args    : object of class TermI.

=cut

=head2 add_relationship

 Title   : add_relationship
 Usage   : add_relationship(RelationshipI relationship): RelationshipI
 Function: Adds a relationship object to the ontology engine.
 Example :
 Returns : Its argument.
 Args    : A RelationshipI object.

=cut

=head2 get_relationships

 Title   : get_relationships
 Usage   : get_relationships(TermI term): RelationshipI
 Function: Retrieves all relationship objects from this ontology engine,
           or all relationships of a term if a term is supplied.
 Example :
 Returns : Array of Bio::Ontology::RelationshipI objects
 Args    : None, or a Bio::Ontology::TermI compliant object for which
           to retrieve the relationships.

=cut

=head2 get_predicate_terms

 Title   : get_predicate_terms
 Usage   : get_predicate_terms(): TermI[]
 Function:
 Example :
 Returns :
 Args    :

=cut

=head2 get_child_terms

 Title   : get_child_terms
 Usage   : get_child_terms(TermI term, TermI predicate_terms): TermI
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

=head2 get_descendant_terms

 Title   : get_descendant_terms
 Usage   : get_descendant_terms(TermI term, TermI rel_types): TermI
 Function: Retrieves all descendant terms of a given term, that
           satisfy a relationship among those that are specified in
           the second argument or undef otherwise.
 Example :
 Returns : Array of TermI objects.
 Args    : First argument is the term of interest, second is the list
           of relationship type terms.

=cut

=head2 get_parent_terms

 Title   : get_parent_terms
 Usage   : get_parent_terms(TermI term, TermI predicate_terms): TermI
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

=head2 get_ancestor_terms

 Title   : get_ancestor_terms
 Usage   : get_ancestor_terms(TermI term, TermI predicate_terms): TermI
 Function: Retrieves all ancestor terms of a given term, that satisfy
           a relationship among those that are specified in the second
           argument or undef otherwise.

 Example :
 Returns : Array of TermI objects.
 Args    : First argument is the term of interest, second is the list
           of relationship type terms.

=cut

=head2 get_leaf_terms

 Title   : get_leaf_terms
 Usage   : get_leaf_terms(): TermI
 Function: Retrieves all leaf terms from the ontology. Leaf term is a
           term w/o descendants.

 Example : @leaf_terms = $obj->get_leaf_terms()
 Returns : Array of TermI objects.
 Args    :

=cut

=head2 get_root_terms()

 Title   : get_root_terms
 Usage   : get_root_terms(): TermI
 Function: Retrieves all root terms from the ontology. Root term is a
           term w/o descendants.

 Example : @root_terms = $obj->get_root_terms()
 Returns : Array of TermI objects.
 Args    :

=cut

=head2 get_all_terms

 Title   : get_all_terms
 Usage   : get_all_terms: TermI
 Function: Retrieves all terms from the ontology.

           We do not mandate an order here in which the terms are
           returned. In fact, the default implementation will return
           them in unpredictable order.

 Example : @terms = $obj->get_all_terms()
 Returns : Array of TermI objects.
 Args    :

=cut


=head2 find_terms

 Title   : find_terms
 Usage   : ($term) = $oe->find_terms(-identifier => "SO:0000263");
 Function: Find term instances matching queries for their attributes.

           An implementation may not support querying for arbitrary
           attributes, but can generally be expected to accept
           -identifier and -name as queries. If both are provided,
           they are implicitly intersected.

 Example :
 Returns : an array of zero or more Bio::Ontology::TermI objects
 Args    : Named parameters. The following parameters should be recognized
           by any implementation:

              -identifier    query by the given identifier
              -name          query by the given name

=cut

=head1 Factory for relationships and terms

=cut

=head2 relationship_factory

 Title   : relationship_factory
 Usage   : $fact = $obj->relationship_factory()
 Function: Get (and set, if the implementation supports it) the object
           factory to be used when relationship objects are created by
           the implementation on-the-fly.

 Example : 
 Returns : value of relationship_factory (a Bio::Factory::ObjectFactoryI
           compliant object)
 Args    : 

=cut

sub relationship_factory{
    return shift->throw_not_implemented();
}

=head2 term_factory

 Title   : term_factory
 Usage   : $fact = $obj->term_factory()
 Function: Get (and set, if the implementation supports it) the object
           factory to be used when term objects are created by
           the implementation on-the-fly.

 Example : 
 Returns : value of term_factory (a Bio::Factory::ObjectFactoryI
           compliant object)
 Args    : 

=cut

sub term_factory{
    return shift->throw_not_implemented();
}

1;
