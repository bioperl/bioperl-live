# $Id$
#
# BioPerl module for Bio::Ontology::RelationshipType  
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

RelationshipType - a relationship type for an ontology


=head1 SYNOPSIS


=head1 DESCRIPTION

This class can be used to model various types of relationships
(such as "IS_A", "PART_OF", "CONTAINS", "FOUND_IN").
 
This class extends Bio::Ontology::Term. 

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
  http://bio.perl.org/bioperl-bugs/

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

package Bio::Ontology::RelationshipType;
use vars qw( @ISA $IS_A $PART_OF $CONTAINS $FOUND_IN );
use strict;
use Bio::Ontology::Term;


use constant PART_OF  => "PART_OF";
use constant IS_A     => "IS_A";
use constant CONTAINS => "CONTAINS";
use constant FOUND_IN => "FOUND_IN";


@ISA = qw( Bio::Ontology::Term );


$IS_A     = Bio::Ontology::RelationshipType->get_instance( IS_A );
$PART_OF  = Bio::Ontology::RelationshipType->get_instance( PART_OF );
$CONTAINS = Bio::Ontology::RelationshipType->get_instance( CONTAINS );
$FOUND_IN = Bio::Ontology::RelationshipType->get_instance( FOUND_IN );



=head2 get_instance

 Title   : get_instance
 Usage   : $IS_A     = Bio::Ontology::RelationshipType->get_instance( "IS_A" );
           $PART_OF  = Bio::Ontology::RelationshipType->get_instance( "PART_OF" );
           $CONTAINS = Bio::Ontology::RelationshipType->get_instance( "CONTAINS" );
           $FOUND_IN = Bio::Ontology::RelationshipType->get_instance( "FOUND_IN" );   
 Function: Factory method to create instances of RelationshipType
 Returns : [Bio::Ontology::RelationshipType]
 Args    : "IS_A" or "PART_OF" or "CONTAINS" or "FOUND_IN" [scalar]

=cut

sub get_instance {
    my ( $class, $value ) = @_;

    my $instance = $class->new();
    $instance->category( "relationship type" );

    if ( $value eq IS_A ) {
        $instance->identifier( IS_A );
        $instance->name( IS_A );
        $instance->definition( IS_A . " relationship type" );
    }
    elsif ( $value eq PART_OF ) {
        $instance->identifier( PART_OF );
        $instance->name( PART_OF );
        $instance->definition( PART_OF . " relationship type" );
    }
    elsif ( $value eq CONTAINS ) {
        $instance->identifier( CONTAINS );
        $instance->name( CONTAINS );
        $instance->definition( CONTAINS . " relationship type" );
    }
    elsif ( $value eq FOUND_IN ) {
        $instance->identifier( FOUND_IN );
        $instance->name( FOUND_IN );
        $instance->definition( FOUND_IN . " relationship type" );
    }
    else {
        my $msg = "Found unknown type of relationship: [" . $value . "]\n";
        $msg .= "Known types are: [" . IS_A . "], [" . PART_OF . "], [" . CONTAINS . "], [" . FOUND_IN . "]";
        $class->throw( $msg );
    }
    
    return $instance;
    
} # get_instance



=head2 init

 Title   : init()
 Usage   : $type->init();   
 Function: Initializes this to all undef and empty lists.
 Returns : 
 Args    :

=cut

sub init {
   my( $self ) = @_;

    $self->{ "_identifier" }  = undef;
    $self->{ "_name" }        = undef;
    $self->{ "_definition" }  = undef;
    $self->{ "_version" }     = undef;
    $self->{ "_is_obsolete" } = undef;
    $self->{ "_comment" }     = undef;
    $self->remove_synonyms();
  
} # init



=head2 equals

 Title   : equals
 Usage   : if ( $type->equals( $other_type ) ) { ...   
 Function: Compares this type to another one, based on string "eq" of
           the "identifier" field
 Returns : true or false
 Args    : [Bio::Ontology::RelationshipType]

=cut

sub equals {
    my( $self, $type ) = @_;
    
    $self->_check_class( $type, "Bio::Ontology::RelationshipType" );
    
    unless ( $self->identifier() && $type->identifier() ) {
        $self->throw( "Cannot compare RelationshipType with a undef identifier" );
    } 
    
    return( $self->identifier() eq $type->identifier() ); 
  
} # equals



=head2 identifier

 Title   : identifier
 Usage   : $term->identifier( "IS_A" );
           or
           print $term->identifier();
 Function: Set/get for the immutable identifier of this Type.
 Returns : The identifier [scalar].
 Args    : The identifier [scalar] (optional).

=cut

sub identifier {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        if ( defined $self->{ "_identifier" } ) {
            $self->throw( "Attempted to change field in immutable object" );
        }
        $self->{ "_identifier" } = $value;
    }

    return $self->{ "_identifier" };

} # identifier




=head2 name

 Title   : name
 Usage   : $term->name( "is a type" );
           or
           print $term->name();
 Function: Set/get for the immutable name of this Type.
 Returns : The name [scalar].
 Args    : The name [scalar] (optional).

=cut

sub name {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        if ( defined $self->{ "_name" } ) {
            $self->throw( "Attempted to change field in immutable object" );
        }
        $self->{ "_name" } = $value;
    }

    return $self->{ "_name" };

} # name





=head2 definition

 Title   : definition
 Usage   : $term->definition( "" );
           or
           print $term->definition();
 Function: Set/get for the immutable definition of this Type.
 Returns : The definition [scalar].
 Args    : The definition [scalar] (optional).

=cut

sub definition {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        if ( defined $self->{ "_definition" } ) {
            $self->throw( "Attempted to change field in immutable object" );
        }
        $self->{ "_definition" } = $value;
    }

    return $self->{ "_definition" };

} # definition



=head2 category

 Title   : category
 Usage   : $term->category( $top );
           or 
           $top = $term->category();
 Function: Set/get for a immutable  relationship between this Type and
           another Term (e.g. the top level of the ontology).
 Returns : The category [TermI].
 Args    : The category [TermI or scalar -- which
           becomes the name of the catagory term] (optional).

=cut

sub category {
    my ( $self, $value ) = @_;
    
    if ( defined $value ) {
        if ( defined $self->{ "_category" } ) {
            $self->throw( "Attempted to change field in immutable object" );
        }
        if ( ! ref( $value ) ) {
            my $term = $self->new();
            $term->name( $value );
            $self->{ "_category" } = $term; 
        }
        elsif ( $value->isa( "Bio::Ontology::TermI" ) ) {
            $self->{ "_category" } = $value; 
        } 
        else {
            $self->throw( "Found [". ref( $value ) 
            . "] where [Bio::Ontology::TermI] or [scalar] expected" );
        }
    }
    
    return $self->{ "_category" };
    
} # category



=head2 version

 Title   : version
 Usage   : $term->version( "1.00" );
           or 
           print $term->version();
 Function: Set/get for immutable version information.
 Returns : The version [scalar].
 Args    : The version [scalar] (optional).

=cut

sub version {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        if ( defined $self->{ "_version" } ) {
            $self->throw( "Attempted to change field in immutable object" );
        }
        $self->{ "_version" } = $value;
    }

    return $self->{ "_version" };
    
} # version



=head2 is_obsolete

 Title   : is_obsolete
 Usage   : $term->is_obsolete( 1 );
           or
           if ( $term->is_obsolete() )
 Function: Set/get for the immutable obsoleteness of this Type.
 Returns : the obsoleteness [0 or 1].
 Args    : the obsoleteness [0 or 1] (optional).

=cut

sub is_obsolete {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
            if ( defined $self->{ "_is_obsolete" } ) {
            $self->throw( "Attempted to change field in immutable object" );
        }
        $self->_is_true_or_false( $value );
        $self->{ "_is_obsolete" } = $value;
    }

    return $self->{ "_is_obsolete" };

} # is_obsolete





=head2 comment

 Title   : comment
 Usage   : $term->comment( "..." );
           or 
           print $term->comment();
 Function: Set/get for an arbitrary immutable comment about this Type.
 Returns : A comment.
 Args    : A comment (optional).

=cut

sub comment {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        if ( defined $self->{ "_comment" } ) {
            $self->throw( "Attempted to change field in immutable object" );
        }
        $self->{ "_comment" } = $value;
    }
   
    return $self->{ "_comment" };
    
} # comment



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

} # _check_type


1;
