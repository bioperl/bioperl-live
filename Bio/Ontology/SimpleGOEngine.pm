# $Id$
#
# BioPerl module for Bio::Ontology::SimpleGOEngine
#
# Cared for by Christian M. Zmasek <czmasek@gnf.org> or <cmzmasek@hotmail.com>
#
# (c) Christian M. Zmasek, czmasek@gnf.org, 2002.
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

SimpleGOEngine - a Ontology Engine for GO implementing OntologyEngineI

=head1 SYNOPSIS

  use Bio::Ontology::simpleGOparser;

  my $parser = Bio::Ontology::simpleGOparser->new
	( -go_defs_file_name    => "/home/czmasek/GO/GO.defs",
	  -components_file_name => "/home/czmasek/GO/component.ontology",
	  -functions_file_name  => "/home/czmasek/GO/function.ontology",
	  -processes_file_name  => "/home/czmasek/GO/process.ontology" );

  my $engine = $parser->parse();

  my $IS_A    = Bio::Ontology::RelationshipType->get_instance( "IS_A" );
  my $PART_OF = Bio::Ontology::RelationshipType->get_instance( "PART_OF" );


=head1 DESCRIPTION

Needs Graph.pm from CPAN.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                         - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.  Bug reports can be submitted via
 email or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR

Christian M. Zmasek

Email: czmasek@gnf.org  or  cmzmasek@hotmail.com

WWW:   http://www.genetics.wustl.edu/eddy/people/zmasek/

Address:

  Genomics Institute of the Novartis Research Foundation
  10675 John Jay Hopkins Drive
  San Diego, CA 92121

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...



package Bio::Ontology::SimpleGOEngine;

use Graph::Directed;

use vars qw( @ISA );
use strict;
use Bio::Root::Root;
use Bio::Ontology::RelationshipType;
use Bio::Ontology::Relationship;
use Bio::Ontology::OntologyEngineI;

use constant TRUE    => 1;
use constant FALSE   => 0;
use constant IS_A    => "IS_A";
use constant PART_OF => "PART_OF";
use constant TERM    => "TERM";
use constant TYPE    => "TYPE";

@ISA = qw( Bio::Ontology::OntologyEngineI );



=head2 new

 Title   : new
 Usage   : $engine = Bio::Ontology::SimpleGOEngine->new()
 Function: Creates a new SimpleGOEngine
 Returns : A new SimpleGOEngine object
 Args    :

=cut

sub new {
    my( $class, @args ) = @_;

    my $self = $class->SUPER::new( @args );

    $self->init();

    return $self;
} # new



=head2 init

 Title   : init()
 Usage   : $engine->init();
 Function: Initializes this Engine.
 Returns :
 Args    :

=cut

sub init {
    my ( $self ) = @_;

    $self->{ "_is_a_relationship" }    = Bio::Ontology::RelationshipType->get_instance( IS_A );
    $self->{ "_part_of_relationship" } = Bio::Ontology::RelationshipType->get_instance( PART_OF );

    $self->graph( Graph::Directed->new() );


} # init



=head2 is_a_relationship

 Title   : is_a_relationship()
 Usage   : $IS_A = $engine->is_a_relationship();
 Function: Returns a Bio::Ontology::RelationshipType object for "is-a"
           relationships
 Returns : Bio::Ontology::RelationshipType set to "IS_A"
 Args    :

=cut

sub is_a_relationship {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->throw( "Attempted to change immutable field" );
    }

    return $self->{ "_is_a_relationship" };
} # is_a_relationship



=head2 part_of_relationship

 Title   : part_of_relationship()
 Usage   : $PART_OF = $engine->part_of_relationship();
 Function: Returns a Bio::Ontology::RelationshipType object for "part-of"
           relationships
 Returns : Bio::Ontology::RelationshipType set to "PART_OF"
 Args    :

=cut

sub part_of_relationship {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->throw( "Attempted to change immutable field" );
    }

    return $self->{ "_part_of_relationship" };
} # part_of_relationship




=head2 add_term

 Title   : add_term
 Usage   : $engine->add_term( $GOterm_obj );
 Function: Adds a Bio::Ontology::GOterm to this engine
 Returns : true
 Args    : Bio::Ontology::GOterm


=cut

sub add_term {
    my ( $self, $term ) = @_;

    $self->_check_class( $term, "Bio::Ontology::GOterm" );

    my $goid = $term->GO_id();

    if ( $self->graph()->has_vertex( $term ) ) {
        $self->throw( "Ontology already contains a GO term with an identifier of \"$goid\"" );
    }

    $self->graph()->add_vertex( $goid );
    $self->graph()->set_attribute( TERM, $goid, $term );

    return TRUE;

} # add_term



=head2 has_term

 Title   : has_term
 Usage   : $engine->has_term( $term );
 Function: Checks whether this engine contains a particular GO term
 Returns : true or false
 Args    : Bio::Ontology::GOterm
           or
           GO term identifier (e.g. "GO:0012345")


=cut

sub has_term {
    my ( $self, $term ) = @_;

    $term = $self->_get_id( $term );

    if ( $self->graph()->has_vertex( $term ) ) {
        return TRUE;
    }
    else {
        return FALSE;
    }

} # has_term



=head2 add_relationship

 Title   : add_relationship
 Usage   : $engine->add_relationship( $relationship );
           $engine->add_relatioship( $parent_obj, $child_obj, $relationship_type );
           $engine->add_relatioship( $parent_id, $child_id, $relationship_type);
 Function: Adds a relationship to this engine
 Returns : true if successfully added, false otherwise
 Args    : GO id, GO id, Bio::Ontology::RelationshipType
           or
           Bio::Ontology::GOterm, Bio::Ontology::GOterm, Bio::Ontology::RelationshipType
           or
           Bio::Ontology::RelationshipI

=cut

# term objs or term ids
sub add_relationship {
    my ( $self, $parent, $child, $type ) = @_;

    if ( scalar( @_ ) == 2 ) {
        $self->_check_class( $parent, "Bio::Ontology::RelationshipI" );
        $child = $parent->child_term();
        $type = $parent->relationship_type();
        $parent = $parent->parent_term();
    }


    $self->_check_class( $type, "Bio::Ontology::RelationshipType" );

    $parent = $self->_get_id( $parent );
    $child = $self->_get_id( $child );

    if ( ! $self->graph()->has_vertex( $child ) ) {
        $self->throw( "Ontology does not contain a GO term with an identifier of \"$child\"" );
    }
    if ( ! $self->graph()->has_vertex( $parent ) ) {
        $self->throw( "Ontology does not contain a GO term with an identifier of \"$parent\"" );
    }

    # This prevents multi graphs.
    if ( $self->graph()->has_edge( $parent, $child ) ) {
        return FALSE;
    }

    $self->graph()->add_edge( $parent, $child );
    $self->graph()->set_attribute( TYPE, $parent, $child, $type );

    return TRUE;

} # add_relationship




=head2 get_relationships


 Title   : get_relationships
 Usage   : $engine->get_relationships( $term );
 Function: Returns all relationships of a term
 Returns : Relationship[]
 Args    : GO id
           or
           Bio::Ontology::GOterm

=cut

sub get_relationships {
    my ( $self, $term ) = @_;

    $term = $self->_get_id( $term );

    if ( ! $self->graph()->has_vertex( $term ) ) {
        $self->throw( "Ontology does not contain a GO term with an identifier of \"$term\"" );
    }

    my @childs  = $self->get_child_terms( $term );
    my @parents = $self->get_parent_terms( $term );

    my @rels = ();

    foreach my $child ( @childs ) {
        my $rel = Bio::Ontology::Relationship->new();
        $rel->parent_term( $self->get_terms( $term ) );
        $rel->child_term( $child );
        $rel->relationship_type( $self->graph()->get_attribute( TYPE, $term, $child->GO_id() ) );
        push( @rels, $rel );
    }
    foreach my $parent ( @parents ) {
        my $rel = Bio::Ontology::Relationship->new();
        $rel->parent_term( $parent );
        $rel->child_term( $self->get_terms( $term ) );
        $rel->relationship_type( $self->graph()->get_attribute( TYPE, $parent->GO_id(), $term ) );
        push( @rels, $rel );
    }

    return @rels;

} # get_relationships



=head2 get_relationship_types

 Title   : get_relationship_types
 Usage   : $engine->get_relationship_types();
 Function: Returns the types of relationships this engine contains
 Returns : Bio::Ontology::RelationshipType[]
 Args    :


=cut

sub get_relationship_types {
    my ( $self ) = @_;

    my @a = ( $self->is_a_relationship(),
              $self->part_of_relationship() );

    return @a;
} # get_relationship_types




=head2 get_child_terms

 Title   : get_child_terms
 Usage   : $engine->get_child_terms( $term_obj, @rel_types );
           $engine->get_child_terms( $term_id, @rel_types );
 Function: Returns the children of this term
 Returns : Bio::Ontology::GOterm[]
 Args    : Bio::Ontology::GOterm, Bio::Ontology::RelationshipType[]
           or
           GO id, Bio::Ontology::RelationshipType[]

           if NO Bio::Ontology::RelationshipType[] is indicated: children
           of ALL types are returned

=cut

sub get_child_terms {
    my ( $self, $term, @types ) = @_;

    return $self->_get_child_parent_terms_helper( $term, TRUE, @types );

} # get_child_terms





=head2 get_descendant_terms

 Title   : get_descendant_terms
 Usage   : $engine->get_descendant_terms( $term_obj, @rel_types );
           $engine->get_descendant_terms( $term_id, @rel_types );
 Function: Returns the descendants of this term
 Returns : Bio::Ontology::GOterm[]
 Args    : Bio::Ontology::GOterm, Bio::Ontology::RelationshipType[]
           or
           GO id, Bio::Ontology::RelationshipType[]

           if NO Bio::Ontology::RelationshipType[] is indicated: descendants
           of ALL types are returned

=cut

sub get_descendant_terms {
    my ( $self, $term, @types ) = @_;

    my %ids = ();
    my @ids = ();

    $term = $self->_get_id( $term );

    if ( ! $self->graph()->has_vertex( $term ) ) {
        $self->throw( "Ontology does not contain a GO term with an identifier of \"$term\"" );
    }

    $self->_get_descendant_terms_helper( $term, \%ids, \@types );

    while( ( my $id ) = each ( %ids ) ) {
        push( @ids, $id );
    }

    return $self->get_terms( @ids );

} # get_descendant_terms




=head2 get_parent_terms

 Title   : get_parent_terms
 Usage   : $engine->get_parent_terms( $term_obj, @rel_types );
           $engine->get_parent_terms( $term_id, @rel_types );
 Function: Returns the parents of this term
 Returns : Bio::Ontology::GOterm[]
 Args    : Bio::Ontology::GOterm, Bio::Ontology::RelationshipType[]
           or
           GO id, Bio::Ontology::RelationshipType[]

           if NO Bio::Ontology::RelationshipType[] is indicated: parents
           of ALL types are returned

=cut

sub get_parent_terms {
    my ( $self, $term, @types ) = @_;

    return $self->_get_child_parent_terms_helper( $term, FALSE, @types );

} # get_parent_terms



=head2 get_ancestor_terms

 Title   : get_ancestor_terms
 Usage   : $engine->get_ancestor_terms( $term_obj, @rel_types );
           $engine->get_ancestor_terms( $term_id, @rel_types );
 Function: Returns the ancestors of this term
 Returns : Bio::Ontology::GOterm[]
 Args    : Bio::Ontology::GOterm, Bio::Ontology::RelationshipType[]
           or
           GO id, Bio::Ontology::RelationshipType[]

           if NO Bio::Ontology::RelationshipType[] is indicated: ancestors
           of ALL types are returned

=cut

sub get_ancestor_terms {
    my ( $self, $term, @types ) = @_;

    my %ids = ();
    my @ids = ();

    $term = $self->_get_id( $term );

    if ( ! $self->graph()->has_vertex( $term ) ) {
        $self->throw( "Ontology does not contain a GO term with an identifier of \"$term\"" );
    }

    $self->_get_ancestor_terms_helper( $term, \%ids, \@types );

    while( ( my $id ) = each ( %ids ) ) {
        push( @ids, $id );
    }

    return $self->get_terms( @ids );

} # get_ancestor_terms





=head2 get_leaf_terms

 Title   : get_leaf_terms
 Usage   : $engine->get_leaf_terms();
 Function: Returns the leaf terms
 Returns : Bio::Ontology::GOterm[]
 Args    :

=cut

sub get_leaf_terms {
    my ( $self ) = @_;

    my @a = $self->graph()->sink_vertices();

    return $self->get_terms( @a );

}



=head2 get_root_terms()

 Title   : get_root_terms
 Usage   : $engine->get_root_terms();
 Function: Returns the root terms
 Returns : Bio::Ontology::GOterm[]
 Args    :

=cut

sub get_root_terms {
    my ( $self ) = @_;


    my @a = $self->graph()->source_vertices();

    return $self->get_terms( @a );

}



=head2 get_term

 Title   : get_term
 Usage   : $engine->get_term( "GO:1234567" );
 Function: Returns a GO term with a given identifier
 Returns : Bio::Ontology::GOterm is present, false otherwise
 Args    : GO id


=cut

sub get_term {
    my ( $self, $id ) = @_;
    if ( $self->graph()->has_vertex( $id ) ) {
        return( $self->graph()->get_attribute( TERM, $id ) );
    }

} # get_term




=head2 get_terms

 Title   : get_terms
 Usage   : $engine->get_term( "GO:1234567", "GO:2234567" );
 Function: Returns a GO terms with given identifiers
 Returns : Bio::Ontology::GOterm[]
 Args    : GO id[]


=cut

sub get_terms {
    my ( $self, @ids ) = @_;

    my @terms = ();

    foreach my $id ( @ids ) {
        if ( $self->graph()->has_vertex( $id ) ) {
            push( @terms, $self->graph()->get_attribute( TERM, $id ) );
        }
    }

    return @terms;

} # get_terms




=head2 each_term

 Title   : each_term
 Usage   : $engine->each_term();
 Function: Returns all terms in this engine
 Returns : Bio::Ontology::GOterm[]
 Args    :

=cut

sub each_term {
    my ( $self ) = @_;

    return( $self->get_terms( $self->graph()->vertices() ) );

} # each_term




=head2 graph

 Title   : graph()
 Usage   : $engine->graph();
 Function: Returns the Graph this engine is based on
 Returns : Graph
 Args    :

=cut

sub graph {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->_check_class( $value, "Graph::Directed" );
        $self->{ "_graph" } = $value;
    }

    return $self->{ "_graph" };
} # graph



# Internal methods
# ----------------


# Checks the correct format of a GO id
# Gets the GO out of a GOterm
sub _get_id {
    my ( $self, $term ) = @_;
    if ( ref( $term ) ) {
        if ( $term->isa( "Bio::Ontology::GOterm" ) ) {
            return $term->GO_id();
        }
        else {
            $self->throw( "Found [" . ref( $term ) . "] where [Bio::Ontology::GOterm] expected" );
        }
    }
    else {
        if ( $term =~ /^GO:\d{7}$/ ) {
            return $term;
        }
        elsif ( $term =~ /^\d{7}$/ ) {
            return "GO:" . $term;
        }
        else {
            $self->throw( "Illegal format for GO identifier [$term]" );
        }
    }
} # _get_id


# Helper for getting children and parent terms
sub _get_child_parent_terms_helper {
    my ( $self, $term, $do_get_child_terms, @types ) = @_;

    foreach my $type ( @types ) {
        $self->_check_class( $type, "Bio::Ontology::RelationshipType" );
    }

    my @relative_terms = ();

    $term = $self->_get_id( $term );
    if ( ! $self->graph()->has_vertex( $term ) ) {
        $self->throw( "Ontology does not contain a GO term with an identifier of \"$term\"" );
    }

    my @all_relative_terms = ();
    if ( $do_get_child_terms ) {
        @all_relative_terms = $self->graph()->successors( $term );
    }
    else {
        @all_relative_terms = $self->graph()->predecessors( $term );
    }

    foreach my $relative ( @all_relative_terms ) {
        if ( scalar( @types ) > 0 ) {
            foreach my $type ( @types ) {
                my $relative_type;
                if ( $do_get_child_terms ) {
                    $relative_type = $self->graph()->get_attribute( TYPE, $term, $relative );
                }
                else {
                    $relative_type = $self->graph()->get_attribute( TYPE, $relative, $term );
                }
                if ( $relative_type->equals( $type ) ) {
                    push( @relative_terms, $relative );
                }
            }
        }
        else {
            push( @relative_terms, $relative );
        }
    }

    return $self->get_terms( @relative_terms );

} # get_child_terms


# Recursive helper
sub _get_descendant_terms_helper {
    my ( $self, $term, $ids_ref, $types_ref ) = @_;

    my @child_terms = $self->get_child_terms( $term, @$types_ref );

    if ( scalar( @child_terms ) < 1 ) {
        return;
    }

    foreach my $child_term ( @child_terms ) {
        my $child_term_id = $child_term->GO_id();
        $ids_ref->{ $child_term_id } = 0;
        $self->_get_descendant_terms_helper( $child_term_id, $ids_ref, $types_ref );
    }

} # _get_descendant_terms_helper


# Recursive helper
sub _get_ancestor_terms_helper {
    my ( $self, $term, $ids_ref, $types_ref ) = @_;

    my @parent_terms = $self->get_parent_terms( $term, @$types_ref );

    if ( scalar( @parent_terms ) < 1 ) {
        return;
    }

    foreach my $parent_term ( @parent_terms ) {
        my $parent_term_id = $parent_term->GO_id();
        $ids_ref->{ $parent_term_id } = 0;
        $self->_get_ancestor_terms_helper( $parent_term_id, $ids_ref, $types_ref );
    }

} # get_ancestor_terms_helper



sub _check_class {
    my ( $self, $value, $expected_class ) = @_;

    if ( ! defined( $value ) ) {
        $self->throw( "Found [undef] where [$expected_class] expected" );
    }
    elsif ( ! ref( $value ) ) {
        $self->throw( "Found [scalar] where [$expected_class] expected" );
    }
    elsif ( ! $value->isa( $expected_class ) ) {
        $self->throw( "Found [" . ref( $value ) . "] where [$expected_class] expected" );
    }

} # _check_class


1;
