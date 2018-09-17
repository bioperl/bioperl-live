#
# BioPerl module for OntologyEngineI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

Bio::Ontology::OntologyEngineI - Interface a minimal Ontology implementation should satisfy

=head1 SYNOPSIS

    # see documentation of methods

=head1 DESCRIPTION

This describes the minimal interface an ontology query engine should
provide.  It intentionally does not make explicit references to the
ontology being a DAG, nor does it mandate that the ontology be a
vocabulary. Rather, it tries to generically express what should be
accessible (queriable) about an ontology.

The idea is to allow for different implementations for different
purposes, which may then differ as to which operations are efficient
and which are not, and how much richer the functionality is on top of
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
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Peter Dimitrov

Email dimitrov@gnf.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Ontology::OntologyEngineI;
use strict;
use Carp;

use base qw(Bio::Root::RootI);

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
 Function: Adds a relationship object to the ontology engine.
 Example :
 Returns : Its argument.
 Args    : A RelationshipI object.

=cut

sub add_relationship{
    shift->throw_not_implemented();
}

=head2 add_relationship_type

 Title   : add_relationship_type
 Usage   : add_relationship_type(scalar,OntologyI ontology)
 Function: Adds a relationshiptype object to the ontology engine.
 Example :
 Returns : 1 on success, undef on failure
 Args    : The name(scalar) of the relationshiptype, and the OntologyI 
           it is to be added to.

=cut

sub add_relationship_type{
    shift->throw_not_implemented();
}

=head2 get_relationship_type

 Title   : get_relationship_type
 Usage   : get_relationship_type(scalar): RelationshipTypeI
 Function: Get a relationshiptype object from the ontology engine.
 Example :
 Returns : A RelationshipTypeI object.
 Args    : The name (scalar) of the RelationshipTypeI object desired.

=cut

sub get_relationship_type{
    shift->throw_not_implemented();
}

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

sub get_relationships{
    shift->throw_not_implemented();
}

=head2 get_predicate_terms

 Title   : get_predicate_terms
 Usage   : get_predicate_terms(): TermI
 Function:
 Example :
 Returns :
 Args    :

=cut

sub get_predicate_terms{
    shift->throw_not_implemented();
}

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

sub get_child_terms{
    shift->throw_not_implemented();
}

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

sub get_descendant_terms{
    shift->throw_not_implemented();
}

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

sub get_parent_terms{
    shift->throw_not_implemented();
}

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

sub get_ancestor_terms{
    shift->throw_not_implemented();
}

=head2 get_leaf_terms

 Title   : get_leaf_terms
 Usage   : get_leaf_terms(): TermI
 Function: Retrieves all leaf terms from the ontology. Leaf term is a
           term w/o descendants.

 Example : @leaf_terms = $obj->get_leaf_terms()
 Returns : Array of TermI objects.
 Args    :

=cut

sub get_leaf_terms{
    shift->throw_not_implemented();
}

=head2 get_root_terms

 Title   : get_root_terms
 Usage   : get_root_terms(): TermI
 Function: Retrieves all root terms from the ontology. Root term is a
           term w/o ancestors.

 Example : @root_terms = $obj->get_root_terms()
 Returns : Array of TermI objects.
 Args    :

=cut

sub get_root_terms{
    shift->throw_not_implemented();
}

=head1 Factory for relationships and terms

=cut

=head2 relationship_factory

 Title   : relationship_factory
 Usage   : $fact = $obj->relationship_factory()
 Function: Get (and set, if the implementation supports it) the object
           factory to be used when relationship objects are created by
           the implementation on-the-fly.

 Example : 
 Returns : value of relationship_factory (a Bio::Factory::ObjectFactory
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
 Returns : value of term_factory (a Bio::Factory::ObjectFactory
           compliant object)
 Args    : 

=cut

sub term_factory{
    return shift->throw_not_implemented();
}

=head1 Decorator Methods

 These methods come with a default implementation that uses the
 abstract methods defined for this interface. This may not be very
 efficient, and hence implementors are encouraged to override these
 methods if they can provide more efficient implementations.

=cut

=head2 get_all_terms

 Title   : get_all_terms
 Usage   : get_all_terms: TermI
 Function: Retrieves all terms from the ontology.

           This is more a decorator method. We provide a default
           implementation here that loops over all root terms and gets
           all descendants for each root term. The overall union of
           terms is then made unique by name and ontology.

           We do not mandate an order here in which the terms are
           returned. In fact, the default implementation will return
           them in unpredictable order.

           Engine implementations that can provide a more efficient
           method for obtaining all terms should definitely override
           this.

 Example : @terms = $obj->get_all_terms()
 Returns : Array of TermI objects.
 Args    :

=cut

sub get_all_terms{
    my $self = shift;
    # get all root nodes
    my @roots = $self->get_root_terms();
    # accumulate all descendants for each root term
    my @terms = map { $self->get_descendant_terms($_); } @roots;
    # add on the root terms themselves
    push(@terms, @roots);
    # make unique by name and ontology
    my %name_map = map { ($_->name."@".$_->ontology->name, $_); } @terms;
    # done 
    return values %name_map;
}

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

sub find_terms{
    my $self = shift;
    my %params = @_;
    @params{ map { lc $_; } keys %params } = values %params; # lowercase keys

    my @terms = grep {
	my $ok = exists($params{-identifier}) ?
	    $_->identifier() eq $params{-identifier} : 1;
	$ok && ((! exists($params{-name})) ||
		($_->name() eq $params{-name}));
    } $self->get_all_terms();
    return @terms;
}

=head1 Experimental API method proposals

 Ontologies are a very new domain in bioperl, and we are not sure yet
 what we will want to do on and with ontologies in which
 situation. The methods from here on downwards are solely API
 descriptions to solicit comment and feedback; the chance of any of
 those being actually implemented already is very slim.

 Disclaimer: As long as an API method stays in this section, it is
 subject to change, possibly even radical change or complete
 deletion. If it's not implemented yet (most likely it isn't),
 implement yourself at your own risk.

 So far for the disclaimer. The reason the API description is here,
 however, is to solicit feedback. Please feel encouraged to share your
 opinion, regardless of what it is (a notable difference of this API
 method to others is that there is actually no working code behind it
 - so the defense line is non-existent for practical purposes).

=cut

=head2 common_ancestor_path

 Title   : common_ancestor_path
 Usage   :
 Function: Get the paths from two terms A and B to term C, such that
           there is no other term D to which A and B would have a shorter
           path, provided there is a term C to which both A and B are
           connected by a path.

           Note that the path to the common ancestor between A and A
           exists, has distance zero, and predicate "identity".

           The search for the common ancestor C can be further
           constrained by supplying a predicate term. If supplied, the
           predicates of the two paths (A,C) and (B,C) must have a
           common ancestor identical to the predicate, or that has a
           path to the predicate.

 Example :
 Returns : The path of the first term to the common ancestor in scalar
           context, and both paths in list context. Paths are
           Bio::Ontology::PathI compliant objects.
 Args    : The two terms (Bio::Ontology::TermI objects), and optionally
           a constraining common predicate (Bio::Ontology::TermI object).
           The latter may also be given as a scalar, in which case it
           is treated as a boolean that, if TRUE, means that the two paths
           must have identical predicates in order to be returned.

=cut

sub common_ancestor_path{
    return shift->throw_not_implemented();
}

1;
