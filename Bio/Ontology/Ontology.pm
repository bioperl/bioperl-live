#
# BioPerl module for Bio::Ontology::Ontology
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

Bio::Ontology::Ontology - standard implementation of an Ontology

=head1 SYNOPSIS

    use Bio::Ontology::Ontology;
    use Bio::Ontology::Term;

    # create ontology object
    my $ont = Bio::Ontology::Ontology->new(-name => "OBF");

    # add terms, relationships ...
    my $bp = Bio::Ontology::Term->new(-identifier => '02', -name => "Bioperl");
    my $obf = Bio::Ontology::Term->new(-identifier => '01', -name => "OBF");
    my $partof = Bio::Ontology::RelationshipType->get_instance("PART_OF");
    $ont->add_term($bp);
    $ont->add_term($obf);
    $ont->add_relationship($bp, $obf, $partof);

    # then query
    my @terms = $ont->get_root_terms(); # "OBF"
    my @desc = $ont->get_descendant_terms($terms[0], $partof); # "Bioperl"
    # ... see methods for other ways to query

    # for advanced users, you can re-use the query engine outside of an
    # ontology to let one instance manage multiple ontologies
    my $ont2 = Bio::Ontology::Ontology->new(-name => "Foundations",
                                            -engine => $ont->engine());


=head1 DESCRIPTION

This is a no-frills implementation of L<Bio::Ontology::OntologyI>.

The query functions are implemented by delegation to an
OntologyEngineI implementation.

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

=head1 AUTHOR - Hilmar Lapp

Email hlapp at gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Ontology::Ontology;
use strict;

# Object preamble - inherits from Bio::Root::Root

#use Bio::Ontology::SimpleOntologyEngine; # loaded dynamically now!

use base qw(Bio::Root::Root Bio::Ontology::OntologyI Bio::AnnotatableI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Ontology::Ontology->new();
 Function: Builds a new Bio::Ontology::Ontology object
 Returns : an instance of Bio::Ontology::Ontology
 Args    : any number of named arguments. The following names will be
           recognized by this module:

            -name         the name of the ontology
            -authority    the name of the authority for the ontology
            -identifier   an identifier for the ontology, if any
            -engine       the Bio::Ontology::OntologyEngineI
                          implementation that this instance should use;
                          default is Bio::Ontology::SimpleOntologyEngine

            See the corresponding get/set methods for further documentation
            on individual properties.

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($name,$auth,$def,$id,$engine) =
        $self->_rearrange([qw(NAME
                          AUTHORITY
                          DEFINITION
                          IDENTIFIER
                          ENGINE)
                      ],
                      @args);
  defined($name) && $self->name($name);
  defined($auth) && $self->authority($auth);
  defined($def) && $self->definition($def);
  defined($id) && $self->identifier($id);
  defined($engine) && $self->engine($engine);

  return $self;
}

=head1 Methods from L<Bio::Ontology::OntologyI>

=cut

=head2 name

 Title   : name
 Usage   : $obj->name($newval)
 Function: Get/set the name of the ontology.
 Example :
 Returns : value of name (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub name{
    my $self = shift;

    return $self->{'name'} = shift if @_;
    return $self->{'name'};
}

=head2 authority

 Title   : authority
 Usage   : $obj->authority($newval)
 Function: Get/set the authority for this ontology, for instance the
           DNS base for the organization granting the name of the
           ontology and identifiers for the terms.

           This attribute is optional and should not generally
           expected by applications to have been set. It is here to
           follow the rules for namespaces, which ontologies serve as
           for terms.

 Example :
 Returns : value of authority (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub authority{
    my $self = shift;

    return $self->{'authority'} = shift if @_;
    return $self->{'authority'};
}

=head2 definition

 Title   : definition
 Usage   : $obj->definition($newval)
 Function: Get/set a descriptive definition of the ontology.
 Example :
 Returns : value of definition (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub definition{
    my $self = shift;

    return $self->{'definition'} = shift if @_;
    return $self->{'definition'};
}

=head2 identifier

 Title   : identifier
 Usage   : $id = $obj->identifier()
 Function: Get an identifier for this ontology.

           This is primarily intended for look-up purposes. The value
           is not modifiable and is determined automatically by the
           implementation.  Also, the identifier's uniqueness will only
           hold within the scope of a particular application's run
           time since it is derived from a memory location.

 Example :
 Returns : value of identifier (a scalar)
 Args    :


=cut

sub identifier{
    my $self = shift;

    if(@_) {
        $self->throw("cannot modify identifier for ".ref($self))
            if exists($self->{'identifier'});
        my $id = shift;
        $self->{'identifier'} = $id if $id;
    }
    if(! exists($self->{'identifier'})) {
        ($self->{'identifier'}) = "$self" =~ /(0x[0-9a-fA-F]+)/;
    }
    return $self->{'identifier'};
}

=head2 close

 Title   : close
 Usage   :
 Function: Release any resources this ontology may occupy. In order
           to efficiently release unused memory or file handles, you
           should call this method once you are finished with an
           ontology.

 Example :
 Returns : TRUE on success and FALSE otherwise
 Args    : none


=cut

sub close{
    my $self = shift;

    # if it is in the ontology store, remove it from there
    my $store = Bio::Ontology::OntologyStore->get_instance();
    $store->remove_ontology($self);
    # essentially we need to dis-associate from the engine here
    $self->engine(undef);
    return 1;
}

=head1 Implementation-specific public methods

=cut

=head2 engine

 Title   : engine
 Usage   : $engine = $obj->engine()
 Function: Get/set the ontology engine to which all the query methods
           delegate.
 Example :
 Returns : an object implementing Bio::Ontology::OntologyEngineI
 Args    : on set, new value (an object implementing
           Bio::Ontology::OntologyEngineI, or  undef)

See L<Bio::Ontology::OntologyEngineI>.

=cut

sub engine{
    my $self = shift;

    if (@_) {
        my $engine = shift;
        if($engine &&
           (! (ref($engine) &&
               $engine->isa("Bio::Ontology::OntologyEngineI")))) {
            $self->throw("object of class ".ref($engine)." does not implement".
                         " Bio::Ontology::OntologyEngineI. Bummer!");
        }
        $self->{'engine'} = $engine;
    } elsif (! exists($self->{'engine'})) {
        # instantiate on demand
        eval {
            # this introduces a dependency on Graph.pm, so load dynamically
            require Bio::Ontology::SimpleOntologyEngine;
        };
        if ($@) {
            $self->throw("failed to load SimpleOntologyEngine, possibly "
                         ."Graph.pm is not installed; either install or supply "
                         ."another OntologyEngineI implementation:\n"
                         .$@);
        }
        $self->{'engine'} = Bio::Ontology::SimpleOntologyEngine->new();
    }
    return $self->{'engine'};
}

=head1 Methods defined in L<Bio::Ontology::OntologyEngineI>

=cut

=head2 add_term

 Title   : add_term
 Usage   : add_term(TermI term): TermI
 Function: Adds TermI object to the ontology engine term store

           If the ontology property of the term object was not set,
           this implementation will set it to itself upon adding the
           term.

 Example : $oe->add_term($term)
 Returns : its argument.
 Args    : object of class TermI.


=cut

sub add_term{
    my $self = shift;
    my $term = shift;

    # set ontology if not set already
    $term->ontology($self) if $term && (! $term->ontology());
    return $self->engine->add_term($term,@_);
}

=head2 add_relationship

 Title   : add_relationship
 Usage   : add_relationship(RelationshipI relationship): RelationshipI
           add_relatioship(TermI subject, TermI predicate, TermI object)
 Function: Adds a relationship object to the ontology engine.
 Example :
 Returns : Its argument.
 Args    : A RelationshipI object.


=cut

sub add_relationship {
  my $self = shift;
  my $rel = shift;

  if($rel && $rel->isa("Bio::Ontology::TermI")) {
    # we need to construct the relationship object on the fly
    my ($predicate,$object) = @_;
    $rel = Bio::Ontology::Relationship->new(
                                            -subject_term   => $rel,
                                            -object_term    => $object,
                                            -predicate_term => $predicate,
                                            -ontology       => $self,
                                           );
  }
  # set ontology if not set already
  $rel->ontology($self) unless $rel->ontology();
  return $self->engine->add_relationship($rel);
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
    my $self = shift;
    return $self->engine->get_relationship_type(@_);
}

=head2 get_relationships

 Title   : get_relationships
 Usage   : get_relationships(TermI term): RelationshipI[]
 Function: Retrieves all relationship objects in the ontology, or all
           relationships of a given term.
 Example :
 Returns : Array of Bio::Ontology::RelationshipI objects
 Args    : Optionally, a Bio::Ontology::TermI compliant object


=cut

sub get_relationships {
  my $self = shift;
  my $term = shift;
  if($term) {
        # we don't need to filter in this case
        return $self->engine->get_relationships($term);
  }
  # else we need to filter by ontology
  return grep { my $ont = $_->ontology;
                # the first condition is a superset of the second, but
                # we add it here for efficiency reasons, as many times
                # it will short-cut to true and is supposedly faster than
                # string comparison
                ($ont == $self) || ($ont->name eq $self->name);
              } $self->engine->get_relationships(@_);
}

=head2 get_predicate_terms

 Title   : get_predicate_terms
 Usage   : get_predicate_terms(): TermI
 Function: Retrieves all relationship types.
 Example :
 Returns : Array of TermI objects
 Args    :


=cut

sub get_predicate_terms{
    my $self = shift;
    
    # skipped Bio::Ontology::Relationship w/o defined Ontology (bug 2573)
    return grep { $_->ontology && ($_->ontology->name eq $self->name)
              } $self->engine->get_predicate_terms(@_);
}

=head2 get_child_terms

 Title   : get_child_terms
 Usage   : get_child_terms(TermI term, TermI predicate_terms): TermI
 Function: Retrieves all child terms of a given term, that satisfy a
           relationship among those that are specified in the second
           argument or undef otherwise. get_child_terms is a special
           case of get_descendant_terms, limiting the search to the
           direct descendants.

           Note that a returned term may possibly be in another
           ontology than this one, because the underlying engine may
           manage multiple ontologies and the relationships of terms
           between them. If you only want descendants within this
           ontology, you need to filter the returned array.

 Example :
 Returns : Array of TermI objects.
 Args    : First argument is the term of interest, second is the list
           of relationship type terms.


=cut

sub get_child_terms{
    return shift->engine->get_child_terms(@_);
}

=head2 get_descendant_terms

 Title   : get_descendant_terms
 Usage   : get_descendant_terms(TermI term, TermI rel_types): TermI
 Function: Retrieves all descendant terms of a given term, that
           satisfy a relationship among those that are specified in
           the second argument or undef otherwise.

           Note that a returned term may possibly be in another
           ontology than this one, because the underlying engine may
           manage multiple ontologies and the relationships of terms
           between them. If you only want descendants within this
           ontology, you need to filter the returned array.

 Example :
 Returns : Array of TermI objects.
 Args    : First argument is the term of interest, second is the list
           of relationship type terms.


=cut

sub get_descendant_terms{
    return shift->engine->get_descendant_terms(@_);
}

=head2 get_parent_terms

 Title   : get_parent_terms
 Usage   : get_parent_terms(TermI term, TermI predicate_terms): TermI
 Function: Retrieves all parent terms of a given term, that satisfy a
           relationship among those that are specified in the second
           argument or undef otherwise. get_parent_terms is a special
           case of get_ancestor_terms, limiting the search to the
           direct ancestors.

           Note that a returned term may possibly be in another
           ontology than this one, because the underlying engine may
           manage multiple ontologies and the relationships of terms
           between them. If you only want descendants within this
           ontology, you need to filter the returned array.

 Example :
 Returns : Array of TermI objects.
 Args    : First argument is the term of interest, second is the list
           of relationship type terms.


=cut

sub get_parent_terms{
    return shift->engine->get_parent_terms(@_);
}

=head2 get_ancestor_terms

 Title   : get_ancestor_terms
 Usage   : get_ancestor_terms(TermI term, TermI predicate_terms): TermI
 Function: Retrieves all ancestor terms of a given term, that satisfy
           a relationship among those that are specified in the second
           argument or undef otherwise.

           Note that a returned term may possibly be in another
           ontology than this one, because the underlying engine may
           manage multiple ontologies and the relationships of terms
           between them. If you only want descendants within this
           ontology, you need to filter the returned array.

 Example :
 Returns : Array of TermI objects.
 Args    : First argument is the term of interest, second is the list
           of relationship type terms.


=cut

sub get_ancestor_terms{
    return shift->engine->get_ancestor_terms(@_);
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
    my $self = shift;
    return grep { my $ont = $_->ontology;
                  # the first condition is a superset of the second, but
                  # we add it here for efficiency reasons, as many times
                  # it will short-cut to true and is supposedly faster than
                  # string comparison
                  ($ont == $self) || ($ont->name eq $self->name);
              } $self->engine->get_leaf_terms(@_);
}

=head2 get_root_terms()

 Title   : get_root_terms
 Usage   : get_root_terms(): TermI
 Function: Retrieves all root terms from the ontology. Root term is a
           term w/o parents.

 Example : @root_terms = $obj->get_root_terms()
 Returns : Array of TermI objects.
 Args    :

=cut

sub get_root_terms{
    my $self = shift;
    return grep { my $ont = $_->ontology;
                  # the first condition is a superset of the second, but
                  # we add it here for efficiency reasons, as many times
                  # it will short-cut to true and is supposedly faster than
                  # string comparison
                  ($ont == $self) || ($ont->name eq $self->name);
              } $self->engine->get_root_terms(@_);
}

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

sub get_all_terms{
    my $self = shift;
    return grep { my $ont = $_->ontology;
                  # the first condition is a superset of the second, but
                  # we add it here for efficiency reasons, as many times
                  # it will short-cut to true and is supposedly faster than
                  # string comparison
                  ($ont == $self) || ($ont->name eq $self->name);
              } $self->engine->get_all_terms(@_);
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
           by any implementations:

              -identifier    query by the given identifier
              -name          query by the given name

=cut

sub find_terms{
    my $self = shift;
    return grep { $_->ontology->name eq $self->name;
              } $self->engine->find_terms(@_);
}

=head2 find_identical_terms

 Title   : find_identical_terms
 Usage   : ($term) = $oe->find_identical_terms($term0);
 Function: Find term instances where name or synonym
           matches the query exactly
 Example :
 Returns : an array of zero or more Bio::Ontology::TermI objects
 Args    : a Bio::Ontology::TermI object

=cut

sub find_identical_terms{
    my $self = shift;
    return grep { $_->ontology->name eq $self->name;
              } $self->engine->find_identical_terms(@_);
}


=head2 find_similar_terms

 Title   : find_similar_terms
 Usage   : ($term) = $oe->find_similar_terms($term0);
 Function: Find term instances where name or synonym, or part of one,
           matches the query.
 Example :
 Returns : an array of zero or more Bio::Ontology::TermI objects
 Args    : a Bio::Ontology::TermI object

=cut

sub find_similar_terms{
    my $self = shift;
    return grep { $_->ontology->name eq $self->name;
              } $self->engine->find_similar_terms(@_);
}

=head2 find_identically_named_terms

 Title   : find_identically_named_terms
 Usage   : ($term) = $oe->find_identically_named_terms($term0);
 Function: Find term instances where names match the query term
           name exactly
 Example :
 Returns : an array of zero or more Bio::Ontology::TermI objects
 Args    : a Bio::Ontology::TermI object

=cut

sub find_identically_named_terms{
    my $self = shift;
    return grep { $_->ontology->name eq $self->name
              } $self->engine->find_identically_named_terms(@_);
}

=head1 Factory for relationships and terms

=cut

=head2 relationship_factory

 Title   : relationship_factory
 Usage   : $fact = $obj->relationship_factory()
 Function: Get (and set, if the engine supports it) the object
           factory to be used when relationship objects are created by
           the implementation on-the-fly.

 Example :
 Returns : value of relationship_factory (a Bio::Factory::ObjectFactoryI
           compliant object)
 Args    :

=cut

sub relationship_factory{
    return shift->engine->relationship_factory(@_);
}

=head2 term_factory

 Title   : term_factory
 Usage   : $fact = $obj->term_factory()
 Function: Get (and set, if the engine supports it) the object
           factory to be used when term objects are created by
           the implementation on-the-fly.

 Example :
 Returns : value of term_factory (a Bio::Factory::ObjectFactoryI
           compliant object)
 Args    :

=cut

sub term_factory{
    return shift->engine->term_factory(@_);
}


=head2 annotation

 Title   : annotation
 Usage   : $annos = $obj->annotation()
 Function: Get/Set the Bio::Annotation::Collection object
           The collection contains Bio::Annotation::SimpleValue
           objects to store header information like the version
           and date present in the header section of an Ontology
           file.

 Example :
 Returns : value of annotation (a Bio::Annotation::Collection
           compliant object)
 Args    : A Bio::Annotation::Collection object (Optional)

=cut

sub annotation{
    my $self = shift;
    $self->{'annotation'} = shift if @_;
    return $self->{'annotation'};
}


#################################################################
# aliases
#################################################################

*get_relationship_types = \&get_predicate_terms;

1;
