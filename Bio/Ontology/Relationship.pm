# $Id$
#
# BioPerl module for Relationship
#
# Cared for by Christian M. Zmasek <czmasek@gnf.org> or <cmzmasek@yahoo.com>
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

Relationship - a relationship for an ontology

=head1 SYNOPSIS

  $rel = Bio::Ontology::Relationship->new( -identifier        => "16847",
                                           -parent_term       => $parent,
                                           -child_term        => $child,
                                           -relationship_type => $type );

=head1 DESCRIPTION

This is a basic implementation of Bio::Ontology::RelationshipI. 

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

Email: czmasek@gnf.org  or  cmzmasek@yahoo.com

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


package Bio::Ontology::Relationship;
use vars qw( @ISA );
use strict;
use Bio::Root::Root;
use Bio::Ontology::TermI;
use Bio::Ontology::RelationshipI;

@ISA = qw( Bio::Ontology::RelationshipI );





=head2 new

 Title   : new
 Usage   : $rel = Bio::Ontology::Relationship->new( -identifier        => "16847",
                                                    -parent_term       => $parent,
                                                    -child_term        => $child,
                                                    -relationship_type => $type );                   
 Function: Creates a new Bio::Ontology::Relationship.
 Returns : A new Bio::Ontology::Relationship object.
 Args    : -identifier            => the identifier of this relationship [scalar]
           -parent_term           => the parent term [Bio::Ontology::TermI]
           -child_term            => the child term [Bio::Ontology::TermI]  
           -relationship_type     => the relationship type [Bio::Ontology::TermI]  

=cut

sub new {

    my( $class, @args ) = @_;
    
    my $self = $class->SUPER::new( @args );
   
    my ( $identifier,
         $parent_term,
         $child_term,
         $relationship_type )
    = $self->_rearrange( [ qw( IDENTIFIER
                               PARENT_TERM
                               CHILD_TERM
                               RELATIONSHIP_TYPE ) ], @args );
   
    $self->init(); 
    
    $identifier        && $self->identifier( $identifier );
    $parent_term       && $self->parent_term( $parent_term );
    $child_term        && $self->child_term( $child_term );
    $relationship_type && $self->relationship_type( $relationship_type );   
                                                    
    return $self;
    
} # new



=head2 init

 Title   : init()
 Usage   : $rel->init();   
 Function: Initializes this Relationship to all undef.
 Returns : 
 Args    :

=cut

sub init {
   my( $self ) = @_;

    $self->{ "_identifier" }        = undef;
    $self->{ "_parent_term" }       = undef;
    $self->{ "_child_term" }        = undef;
    $self->{ "_relationship_type" } = undef;
  
} # init



=head2 identifier

 Title   : identifier
 Usage   : $rel->identifier( "100050" );
           or
           print $rel->identifier();
 Function: Set/get for the identifier of this Relationship.
 Returns : The identifier [scalar].
 Args    : The identifier [scalar] (optional).

=cut

sub identifier {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_identifier" } = $value;
    }

    return $self->{ "_identifier" };
} # identifier




=head2 parent_term

 Title   : parent_term
 Usage   : $rel->parent_term( $parent );
           or
           $parent = $rel->parent_term();
 Function: Set/get for the parent term of this Relationship.
 Returns : The parent term [Bio::Ontology::TermI].
 Args    : The parent term [Bio::Ontology::TermI] (optional).

=cut

sub parent_term {
    my ( $self, $term ) = @_;
  
    if ( defined $term ) {
        $self->_check_class( $term, "Bio::Ontology::TermI" );
        $self->{ "_parent_term" } = $term;
    }

    return $self->{ "_parent_term" };
    
} # parent_term



=head2 child_term

 Title   : child_term
 Usage   : $rel->child_term( $child );
           or
           $child = $rel->child_term();
 Function: Set/get for the child term of this Relationship.
 Returns : The child term [Bio::Ontology::TermI].
 Args    : The child term [Bio::Ontology::TermI] (optional).

=cut

sub child_term {
    my ( $self, $term ) = @_;
  
    if ( defined $term ) {
        $self->_check_class( $term, "Bio::Ontology::TermI" );
        $self->{ "_child_term" } = $term;
    }

    return $self->{ "_child_term" };
}



=head2 relationship_type

 Title   : relationship_type
 Usage   : $rel->relationship_type( $type );
           or
           $type = $rel->relationship_type();
 Function: Set/get for the relationship type of this relationship.
 Returns : The relationship type [Bio::Ontology::TermI].
 Args    : The relationship type [Bio::Ontology::TermI] (optional).

=cut

sub relationship_type {
    my ( $self, $term ) = @_;
  
    if ( defined $term ) {
        $self->_check_class( $term, "Bio::Ontology::TermI" );
        $self->{ "_relationship_type" } = $term;
    }

    return $self->{ "_relationship_type" };
}



=head2 to_string

 Title   : to_string()
 Usage   : print $rel->to_string();
 Function: to_string method for Relationship.
 Returns : A string representation of this Relationship.
 Args    :

=cut

sub to_string {
    my( $self ) = @_;
    
    local $^W = 0;

    my $s = "";

    $s .= "-- Identifier:\n";
    $s .= $self->identifier()."\n";
    $s .= "-- Parent Term Identifier:\n";
    $s .= $self->parent_term()->identifier()."\n";
    $s .= "-- Child Term Identifier:\n";
    $s .= $self->child_term()->identifier()."\n";
    $s .= "-- Relationship Type Identifier:\n";
    $s .= $self->relationship_type()->identifier();
    
    return $s;
    
} # to_string



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
