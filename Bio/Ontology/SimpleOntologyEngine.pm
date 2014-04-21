#
# BioPerl module for Bio::Ontology::SimpleOntologyEngine
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Peter Dimitrov <dimitrov@gnf.org>
#
# Copyright Peter Dimitrov
# (c) Peter Dimitrov, dimitrov@gnf.org, 2002.
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
# POD documentation - main docs before the code

=head1 NAME

Bio::Ontology::SimpleOntologyEngine - Implementation of OntologyEngineI interface

=head1 SYNOPSIS

  my $soe = Bio::Ontology::SimpleOntologyEngine->new;

=head1 DESCRIPTION

This is a "simple" implementation of Bio::Ontology::OntologyEngineI.

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

=head1 CONTRIBUTORS

Hilmar Lapp, hlapp at gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Ontology::SimpleOntologyEngine;
use strict;
use Carp;
use Bio::Ontology::RelationshipFactory;
use Data::Dumper;

use base qw(Bio::Root::Root Bio::Ontology::OntologyEngineI);

=head2 new

 Title   : new
 Usage   : $soe = Bio::Ontology::SimpleOntologyEngine->new;
 Function: Initializes the ontology engine.
 Example : $soe = Bio::Ontology::SimpleOntologyEngine->new;
 Returns : Object of class SimpleOntologyEngine.
 Args    :

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    #   my %param = @args;

    $self->_term_store(                  {} );
    $self->_relationship_store(          {} );
    $self->_inverted_relationship_store( {} );
    $self->_relationship_type_store(     {} );
    $self->_instantiated_terms_store(    {} );

    # set defaults for the factories
    $self->relationship_factory(
        Bio::Ontology::RelationshipFactory->new( -type => "Bio::Ontology::Relationship" ) );
    return $self;
}

=head2 _instantiated_terms_store

 Title   : _instantiated_terms_store
 Usage   : $obj->_instantiated_terms_store($newval)
 Function:
 Example :
 Returns : hash
 Args    : empty hash

=cut

sub _instantiated_terms_store {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{'_instantiated_terms_store'} = $value;
    }
    return $self->{'_instantiated_terms_store'};
}

=head2 mark_instantiated

 Title   : mark_instantiated
 Usage   : $self->mark_instantiated(TermI terms): TermI
 Function: Marks TermI objects as fully instantiated,
           allowing for proper counting of the number of terms in the term store.
           The TermI objects has to be already stored in the term store in order
           to be marked.
 Example : $self->mark_instantiated($term);
 Returns : its argument or throws an exception if a term is not
           in the term store.
 Args    : array of objects of class TermI.

=cut

sub mark_instantiated {
    my ( $self, @terms ) = @_;

    foreach my $term (@terms) {
        $self->throw( "term " . $term->identifier . " not in the term store\n" )
            if !defined $self->_term_store->{ $term->identifier };
        $self->_instantiated_terms_store->{ $term->identifier } = 1;
    }

    return @terms;
}

=head2 mark_uninstantiated

 Title   : mark_uninstantiated
 Usage   : $self->mark_uninstantiated(TermI terms): TermI
 Function: Marks TermI objects as not fully instantiated,
 Example : $self->mark_uninstantiated($term);
 Returns : its argument or throws an exception if a term is not
           in the term store(if the term is not marked it does nothing).
 Args    : array of objects of class TermI.

=cut

sub mark_uninstantiated {
    my ( $self, @terms ) = @_;

    foreach my $term (@terms) {
        $self->throw( "term " . $term->identifier . " not in the term store\n" )
            if !defined $self->_term_store->{ $term->identifier };
        delete $self->_instantiated_terms_store->{ $term->identifier }
            if defined $self->_instantiated_terms_store->{ $term->identifier };
    }

    return @terms;
}

=head2 _term_store

 Title   : term_store
 Usage   : $obj->_term_store($newval)
 Function:
 Example :
 Returns : reference to an array of Bio::Ontology::TermI objects
 Args    : reference to an array of Bio::Ontology::TermI objects

=cut

sub _term_store {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        if ( defined $self->{'_term_store'} ) {
            $self->throw("_term_store already defined\n");
        } else {
            $self->{'_term_store'} = $value;
        }
    }

    return $self->{'_term_store'};
}

=head2 add_term

 Title   : add_term
 Usage   : add_term(TermI term): TermI
 Function: Adds TermI object to the ontology engine term store.
 Marks the term fully instantiated by default.
 Example : $soe->add_term($term)
 Returns : its argument.
 Args    : object of class TermI.

=cut

sub add_term {
    my ( $self, $term ) = @_;
    my $term_store = $self->_term_store;

    if ( defined $term_store->{ $term->identifier } && $self->_instantiated_terms_store->{ $term->identifier }) {
        $self->throw( "term " . $term->identifier . " already defined\n" );
    } else {
        $term_store->{ $term->identifier } = $term;
        $self->_instantiated_terms_store->{ $term->identifier } = 1;
    }

    return $term;
}

=head2 get_term_by_identifier

 Title   : get_term_by_identifier
 Usage   : get_term_by_identifier(String id): TermI
 Function: Retrieves terms from the term store by their identifier
           field, or an empty list if not there.
 Example : $term = $soe->get_term_by_identifier("IPR000001");
 Returns : An array of zero or more Bio::Ontology::TermI objects.
 Args    : An array of identifier strings

=cut

sub get_term_by_identifier {
    my ( $self, @ids ) = @_;
    my @ans = ();

    foreach my $id (@ids) {
        my $term = $self->_term_store->{$id};
        push @ans, $term if defined $term;
    }

    return @ans;
}

=head2 _get_number_rels

 Title   : get_number_rels
 Usage   :
 Function:
 Example :
 Returns :
 Args    :

=cut

sub _get_number_rels {
    my ($self) = @_;
    my $num_rels = 0;

    foreach my $entry ( $self->_relationship_store ) {
        $num_rels += scalar keys %$entry;
    }
    return $num_rels;
}

=head2 _get_number_terms

 Title   : _get_number_terms
 Usage   :
 Function:
 Example :
 Returns :
 Args    :

=cut

sub _get_number_terms {
    my ($self) = @_;

    return scalar $self->_filter_unmarked( values %{ $self->_term_store } );

}

=head2 _relationship_store

 Title   : _storerelationship_store
 Usage   : $obj->relationship_store($newval)
 Function:
 Example :
 Returns : reference to an array of Bio::Ontology::TermI objects
 Args    : reference to an array of Bio::Ontology::TermI objects

=cut

sub _relationship_store {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        if ( defined $self->{'_relationship_store'} ) {
            $self->throw("_relationship_store already defined\n");
        } else {
            $self->{'_relationship_store'} = $value;
        }
    }

    return $self->{'_relationship_store'};
}

=head2 _inverted_relationship_store

 Title   : _inverted_relationship_store
 Usage   :
 Function:
 Example :
 Returns : reference to an array of Bio::Ontology::TermI objects
 Args    : reference to an array of Bio::Ontology::TermI objects

=cut

sub _inverted_relationship_store {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        if ( defined $self->{'_inverted_relationship_store'} ) {
            $self->throw("_inverted_relationship_store already defined\n");
        } else {
            $self->{'_inverted_relationship_store'} = $value;
        }
    }

    return $self->{'_inverted_relationship_store'};
}

=head2 _relationship_type_store

 Title   : _relationship_type_store
 Usage   : $obj->_relationship_type_store($newval)
 Function:
 Example :
 Returns : reference to an array of Bio::Ontology::RelationshipType objects
 Args    : reference to an array of Bio::Ontology::RelationshipType objects

=cut

sub _relationship_type_store {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        if ( defined $self->{'_relationship_type_store'} ) {
            $self->throw("_relationship_type_store already defined\n");
        } else {
            $self->{'_relationship_type_store'} = $value;
        }
    }

    return $self->{'_relationship_type_store'};
}

=head2 _add_relationship_simple

 Title   : _add_relationship_simple
 Usage   :
 Function:
 Example :
 Returns :
 Args    :

=cut

sub _add_relationship_simple {
    my ( $self, $store, $rel, $inverted ) = @_;

    my $subject = $rel->subject_term
        or $self->throw('cannot add relationship, relationship has no subject_term');
    my $object  = $rel->object_term
        or $self->throw('cannot add relationship, relationship has no object_term');

    my ( $parent_id, $child_id ) = ( $object->identifier, $subject->identifier );
    ( $parent_id, $child_id ) = ( $child_id, $parent_id ) if $inverted;

    if (     defined $store->{$parent_id}
         &&  defined $store->{$parent_id}->{$child_id}
         && $store->{$parent_id}->{$child_id}->name ne $rel->predicate_term->name
     ) {
        $self->throw( "relationship "
                . $rel->predicate_term->name
                . " between "
                . $parent_id . " and "
                . $child_id
                . " already defined as "
                . $store->{$parent_id}->{$child_id}->name
                . "\n" );
    }

    # all is well if we get here
    $store->{$parent_id}->{$child_id} = $rel->predicate_term;
}

=head2 add_relationship

 Title   : add_relationship
 Usage   : add_relationship(RelationshipI relationship): RelationshipI
 Function: Adds a relationship object to the ontology engine.
 Example :
 Returns : Its argument.
 Args    : A RelationshipI object.

=cut

sub add_relationship {
    my ( $self, $rel ) = @_;

    $self->_add_relationship_simple( $self->_relationship_store,          $rel, 0 );
    $self->_add_relationship_simple( $self->_inverted_relationship_store, $rel, 1 );
    $self->_relationship_type_store->{ $self->_unique_termid( $rel->predicate_term ) } =
        $rel->predicate_term;

    return $rel;
}

=head2 get_relationships

 Title   : get_relationships
 Usage   : get_relationships(): RelationshipI
 Function: Retrieves all relationship objects.
 Example :
 Returns : Array of RelationshipI objects
 Args    :

=cut

sub get_relationships {
    my $self = shift;
    my $term = shift;
    my @rels;
    my $store   = $self->_relationship_store;
    my $relfact = $self->relationship_factory();

    my @parent_ids = $term
        ?

        # if a term is supplied then only get the term's parents
        ( map { $_->identifier(); } $self->get_parent_terms($term) )
        :

        # otherwise use all parent ids
        ( keys %{$store} );

    # add the term as a parent too if one is supplied
    push( @parent_ids, $term->identifier ) if $term;

    foreach my $parent_id (@parent_ids) {
        my $parent_entry = $store->{$parent_id};

        # if a term is supplied, add a relationship for the parent to the term
        # except if the parent is the term itself (we added that one before)
        if ( $term && ( $parent_id ne $term->identifier() ) ) {
            my @parent_terms = $self->get_term_by_identifier($parent_id);
            foreach my $parent_term (@parent_terms) {
                push(
                    @rels,
                    $relfact->create_object(
                        -object_term    => $parent_term,
                        -subject_term   => $term,
                        -predicate_term => $parent_entry->{ $term->identifier },
                        -ontology       => $term->ontology()
                    )
                );
            }

        } else {

            # otherwise, i.e., no term supplied, or the parent equals the
            # supplied term
            my @parent_terms = $term ? ($term) : $self->get_term_by_identifier($parent_id);
            foreach my $child_id ( keys %$parent_entry ) {
                my $rel_info = $parent_entry->{$child_id};
                my ($subj_term) = $self->get_term_by_identifier($child_id);

                foreach my $parent_term (@parent_terms) {
                    push(
                        @rels,
                        $relfact->create_object(
                            -object_term    => $parent_term,
                            -subject_term   => $subj_term,
                            -predicate_term => $rel_info,
                            -ontology       => $parent_term->ontology
                        )
                    );
                }
            }
        }
    }

    return @rels;
}

=head2 get_all_relationships

 Title   : get_all_relationships
 Usage   : get_all_relationships(): RelationshipI
 Function: Retrieves all relationship objects.
 Example :
 Returns : Array of RelationshipI objects
 Args    :

=cut

sub get_all_relationships {
    return shift->get_relationships();
}

=head2 get_predicate_terms

 Title   : get_predicate_terms
 Usage   : get_predicate_terms(): TermI
 Function: Retrives all relationship types stored in the engine
 Example :
 Returns : reference to an array of Bio::Ontology::RelationshipType objects
 Args    :

=cut

sub get_predicate_terms {
    my ($self) = @_;

    return values %{ $self->_relationship_type_store };
}

=head2 _is_rel_type

 Title   : _is_rel_type
 Usage   :
 Function:
 Example :
 Returns :
 Args    :

=cut

sub _is_rel_type {
    my ( $self, $term, @rel_types ) = @_;

    foreach my $rel_type (@rel_types) {
        if ( $rel_type->identifier || $term->identifier ) {
            return 1 if $rel_type->identifier eq $term->identifier;
        } else {
            return 1 if $rel_type->name eq $term->name;
        }
    }

    return 0;
}

=head2 _typed_traversal

 Title   : _typed_traversal
 Usage   :
 Function:
 Example :
 Returns :
 Args    :

=cut

sub _typed_traversal {
    my ( $self, $rel_store, $level, $term_id, @rel_types ) = @_;
    return if !defined( $rel_store->{$term_id} );
    my %parent_entry = %{ $rel_store->{$term_id} };
    my @children     = keys %parent_entry;

    my @ans;

    if ( @rel_types > 0 ) {
        @ans = ();

        foreach my $child_id (@children) {
            push @ans, $child_id
                if $self->_is_rel_type( $rel_store->{$term_id}->{$child_id}, @rel_types );
        }
    } else {
        @ans = @children;
    }
    if ( $level < 1 ) {
        my @ans1 = ();

        foreach my $child_id (@ans) {
            push @ans1, $self->_typed_traversal( $rel_store, $level - 1, $child_id, @rel_types )
                if defined $rel_store->{$child_id};
        }
        push @ans, @ans1;
    }

    return @ans;
}

=head2 get_child_terms

 Title   : get_child_terms
 Usage   : get_child_terms(TermI term, TermI predicate_terms): TermI
           get_child_terms(TermI term, RelationshipType predicate_terms): TermI
 Function: Retrieves all child terms of a given term, that satisfy a
           relationship among those that are specified in the second
           argument or undef otherwise. get_child_terms is a special
           case of get_descendant_terms, limiting the search to the
           direct descendants.
 Example :
 Returns : Array of TermI objects.
 Args    : First argument is the term of interest, second is the list of
           relationship type terms.

=cut

sub get_child_terms {
    my ( $self, $term, @relationship_types ) = @_;

    $self->throw("must provide TermI compliant object")
        unless defined($term) && $term->isa("Bio::Ontology::TermI");

    return $self->_filter_unmarked(
        $self->get_term_by_identifier(
            $self->_typed_traversal(
                $self->_relationship_store, 1, $term->identifier, @relationship_types
            )
        )
    );
}

=head2 get_descendant_terms

 Title   : get_descendant_terms
 Usage   : get_descendant_terms(TermI term, TermI rel_types): TermI
           get_child_terms(TermI term, RelationshipType predicate_terms): TermI
 Function: Retrieves all descendant terms of a given term, that
           satisfy a relationship among those that are specified in
           the second argument or undef otherwise. Uses
           _typed_traversal to find all descendants.

 Example :
 Returns : Array of TermI objects.
 Args    : First argument is the term of interest, second is the list of
           relationship type terms.

=cut

sub get_descendant_terms {
    my ( $self, $term, @relationship_types ) = @_;

    $self->throw("must provide TermI compliant object")
        unless defined($term) && $term->isa("Bio::Ontology::TermI");

    return $self->_filter_unmarked(
        $self->_filter_repeated(
            $self->get_term_by_identifier(
                $self->_typed_traversal(
                    $self->_relationship_store, 0, $term->identifier, @relationship_types
                )
            )
        )
    );
}

=head2 get_parent_terms

 Title   : get_parent_terms
 Usage   : get_parent_terms(TermI term, TermI predicate_terms): TermI
           get_child_terms(TermI term, RelationshipType predicate_terms): TermI
 Function: Retrieves all parent terms of a given term, that satisfy a
           relationship among those that are specified in the second
           argument or undef otherwise. get_parent_terms is a special
           case of get_ancestor_terms, limiting the search to the
           direct ancestors.

 Example :
 Returns : Array of TermI objects.
 Args    : First argument is the term of interest, second is the list of relationship type terms.

=cut

sub get_parent_terms {
    my ( $self, $term, @relationship_types ) = @_;
    $self->throw("term must be a valid object, not undef") unless defined $term;

    return $self->_filter_unmarked(
        $self->get_term_by_identifier(
            $self->_typed_traversal(
                $self->_inverted_relationship_store,
                1, $term->identifier, @relationship_types
            )
        )
    );
}

=head2 get_ancestor_terms

 Title   : get_ancestor_terms
 Usage   : get_ancestor_terms(TermI term, TermI predicate_terms): TermI
           get_child_terms(TermI term, RelationshipType predicate_terms): TermI
 Function: Retrieves all ancestor terms of a given term, that satisfy
           a relationship among those that are specified in the second
           argument or undef otherwise. Uses _typed_traversal to find
           all ancestors.

 Example :
 Returns : Array of TermI objects.
 Args    : First argument is the term of interest, second is the list
           of relationship type terms.

=cut

sub get_ancestor_terms {
    my ( $self, $term, @relationship_types ) = @_;
    $self->throw("term must be a valid object, not undef") unless defined $term;

    return $self->_filter_unmarked(
        $self->_filter_repeated(
            $self->get_term_by_identifier(
                $self->_typed_traversal(
                    $self->_inverted_relationship_store, 0,
                    $term->identifier,                   @relationship_types
                )
            )
        )
    );
}

=head2 get_leaf_terms

 Title   : get_leaf_terms
 Usage   : get_leaf_terms(): TermI
 Function: Retrieves all leaf terms from the ontology. Leaf term is a term w/o descendants.
 Example : @leaf_terms = $obj->get_leaf_terms()
 Returns : Array of TermI objects.
 Args    :

=cut

sub get_leaf_terms {
    my ($self) = @_;
    my @leaf_terms;

    foreach my $term ( values %{ $self->_term_store } ) {
        push @leaf_terms, $term
            if !defined $self->_relationship_store->{ $term->identifier }
                && defined $self->_instantiated_terms_store->{ $term->identifier };
    }

    return @leaf_terms;
}

=head2 get_root_terms

 Title   : get_root_terms
 Usage   : get_root_terms(): TermI
 Function: Retrieves all root terms from the ontology. Root term is a term w/o descendants.
 Example : @root_terms = $obj->get_root_terms()
 Returns : Array of TermI objects.
 Args    :

=cut

sub get_root_terms {
    my ($self) = @_;
    my @root_terms;

    foreach my $term ( values %{ $self->_term_store } ) {
        push @root_terms, $term
            if !defined $self->_inverted_relationship_store->{ $term->identifier }
                && defined $self->_instantiated_terms_store->{ $term->identifier };
    }

    return @root_terms;
}

=head2 _filter_repeated

 Title   : _filter_repeated
 Usage   : @lst = $self->_filter_repeated(@old_lst);
 Function: Removes repeated terms
 Example :
 Returns : List of unique TermI objects
 Args    : List of TermI objects

=cut

sub _filter_repeated {
    my ( $self, @args ) = @_;
    my %h;

    foreach my $element (@args) {
        $h{ $element->identifier } = $element if !defined $h{ $element->identifier };
    }

    return values %h;
}

=head2 get_all_terms

 Title   : get_all_terms
 Usage   : get_all_terms(): TermI
 Function: Retrieves all terms currently stored in the ontology.
 Example : @all_terms = $obj->get_all_terms()
 Returns : Array of TermI objects.
 Args    :

=cut

sub get_all_terms {
    my ($self) = @_;

    return $self->_filter_unmarked( values %{ $self->_term_store } );
}

=head2 find_terms

 Title   : find_terms
 Usage   : ($term) = $oe->find_terms(-identifier => "SO:0000263");
 Function: Find term instances matching queries for their attributes.

           This implementation can efficiently resolve queries by
           identifier.

 Example :
 Returns : an array of zero or more Bio::Ontology::TermI objects
 Args    : Named parameters. The following parameters should be recognized
           by any implementations:

              -identifier    query by the given identifier
              -name          query by the given name

=cut

sub find_terms {
    my ( $self, @args ) = @_;
    my @terms;

    my ( $id, $name ) = $self->_rearrange( [qw(IDENTIFIER NAME)], @args );

    if ( defined($id) ) {
        @terms = $self->get_term_by_identifier($id);
    } else {
        @terms = $self->get_all_terms();
    }
    if ( defined($name) ) {
        @terms = grep { $_->name() eq $name; } @terms;
    }
    return @terms;
}

=head2 relationship_factory

 Title   : relationship_factory
 Usage   : $fact = $obj->relationship_factory()
 Function: Get/set the object factory to be used when relationship
           objects are created by the implementation on-the-fly.

 Example :
 Returns : value of relationship_factory (a Bio::Factory::ObjectFactoryI
           compliant object)
 Args    : on set, a Bio::Factory::ObjectFactoryI compliant object

=cut

sub relationship_factory {
    my $self = shift;

    return $self->{'relationship_factory'} = shift if @_;
    return $self->{'relationship_factory'};
}

=head2 term_factory

 Title   : term_factory
 Usage   : $fact = $obj->term_factory()
 Function: Get/set the object factory to be used when term objects are
           created by the implementation on-the-fly.

           Note that this ontology engine implementation does not
           create term objects on the fly, and therefore setting this
           attribute is meaningless.

 Example :
 Returns : value of term_factory (a Bio::Factory::ObjectFactoryI
           compliant object)
 Args    : on set, a Bio::Factory::ObjectFactoryI compliant object

=cut

sub term_factory {
    my $self = shift;

    if (@_) {
        $self->warn(
            "setting term factory, but " . ref($self) . " does not create terms on-the-fly" );
        return $self->{'term_factory'} = shift;
    }
    return $self->{'term_factory'};
}

=head2 _filter_unmarked

 Title   : _filter_unmarked
 Usage   : _filter_unmarked(TermI terms): TermI
 Function: Removes the uninstantiated terms from the list of terms
 Example :
 Returns : array of fully instantiated TermI objects
 Args    : array of TermI objects

=cut

sub _filter_unmarked {
    my ( $self, @terms ) = @_;
    my @filtered_terms = ();

    if ( scalar(@terms) >= 1 ) {
        foreach my $term (@terms) {
            push @filtered_terms, $term
                if defined $self->_instantiated_terms_store->{ $term->identifier };
        }
    }

    return @filtered_terms;
}

=head2 remove_term_by_id

 Title   : remove_term_by_id
 Usage   : remove_term_by_id(String id): TermI
 Function: Removes TermI object from the ontology engine using the
           string id as an identifier. Current implementation does not
           enforce consistency of the relationships using that term.
 Example : $term = $soe->remove_term_by_id($id);
 Returns : Object of class TermI or undef if not found.
 Args    : The string identifier of a term.

=cut

sub remove_term_by_id {
    my ( $self, $id ) = @_;

    if ( $self->get_term_by_identifier($id) ) {
        my $term = $self->{_term_store}->{$id};
        delete $self->{_term_store}->{$id};
        return $term;
    } else {
        $self->warn("Term with id '$id' is not in the term store");
        return;
    }
}

=head2 to_string

 Title   : to_string
 Usage   : print $sv->to_string();
 Function: Currently returns formatted string containing the number of
           terms and number of relationships from the ontology engine.
 Example : print $sv->to_string();
 Returns :
 Args    :

=cut

sub to_string {
    my ($self) = @_;
    my $s = "";

    $s .= "-- # Terms:\n";
    $s .= scalar( $self->get_all_terms ) . "\n";
    $s .= "-- # Relationships:\n";
    $s .= $self->_get_number_rels . "\n";

    return $s;
}

=head2 _unique_termid

 Title   : _unique_termid
 Usage   :
 Function: Returns a string that can be used as ID using fail-over
           approaches.

           If the identifier attribute is not set, it uses the
           combination of name and ontology name, provided both are
           set. If they are not, it returns the name alone.

           Note that this is a private method. Call from inheriting
           classes but not from outside.

 Example :
 Returns : a string
 Args    : a Bio::Ontology::TermI compliant object

=cut

sub _unique_termid {
    my $self = shift;
    my $term = shift;

    return $term->identifier() if $term->identifier();
    my $id = $term->ontology->name() if $term->ontology();
    if ($id) {
        $id .= '|';
    } else {
        $id = '';
    }
    $id .= $term->name();
}

#################################################################
# aliases
#################################################################

*get_relationship_types = \&get_predicate_terms;

1;
