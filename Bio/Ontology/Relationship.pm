#
# BioPerl module for Relationship
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Christian M. Zmasek <czmasek-at-burnham.org> or <cmzmasek@yahoo.com>
#
# (c) Christian M. Zmasek, czmasek-at-burnham.org, 2002.
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

Bio::Ontology::Relationship - a relationship for an ontology

=head1 SYNOPSIS

  $rel = Bio::Ontology::Relationship->new( -identifier     => "16847",
                                           -subject_term   => $subj,
                                           -object_term    => $obj,
                                           -predicate_term => $pred );

=head1 DESCRIPTION

This is a basic implementation of Bio::Ontology::RelationshipI. 

The terminology we use here is the one commonly used for ontologies,
namely the triple of (subject, predicate, object), which in addition
is scoped in a namespace (ontology). It is called triple because it is
a tuple of three ontology terms.

There are other terminologies in use for expressing relationships. For
those who it helps to better understand the concept, the triple of
(child, relationship type, parent) would be equivalent to the
terminology chosen here, disregarding the question whether the notion
of parent and child is sensible in the context of the relationship
type or not. Especially in the case of ontologies with a wide variety
of predicates the parent/child terminology and similar ones can
quickly become ambiguous (e.g., A synthesises B), meaningless (e.g., A
binds B), or even conflicting (e.g., A is-parent-of B), and are
therefore strongly discouraged.

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

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Christian M. Zmasek

Email: czmasek-at-burnham.org  or  cmzmasek@yahoo.com

WWW:   http://monochrome-effect.net/

Address: 

  Genomics Institute of the Novartis Research Foundation
  10675 John Jay Hopkins Drive
  San Diego, CA 92121

=head1 CONTRIBUTORS

 Hilmar Lapp, email: hlapp at gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Ontology::Relationship;
use strict;
use Bio::Ontology::TermI;

use base qw(Bio::Root::Root Bio::Ontology::RelationshipI);




=head2 new

 Title   : new
 Usage   : $rel = Bio::Ontology::Relationship->new(-identifier   => "16847",
                                                   -subject_term => $subject,
                                                   -object_term  => $object,
                                                   -predicate_term => $type );
 Function: Creates a new Bio::Ontology::Relationship.
 Returns : A new Bio::Ontology::Relationship object.
 Args    : -identifier     => the identifier of this relationship [scalar]
           -subject_term   => the subject term [Bio::Ontology::TermI]
           -object_term    => the object term [Bio::Ontology::TermI]  
           -predicate_term => the predicate term [Bio::Ontology::TermI]

=cut

sub new {

    my( $class, @args ) = @_;
    
    my $self = $class->SUPER::new( @args );
   
    my ( $identifier,
         $subject_term,
			$child,        # for backwards compatibility
         $object_term,
			$parent,       # for backwards compatibility
         $predicate_term,
			$reltype,      # for backwards compatibility
			$ont)
	= $self->_rearrange( [qw( IDENTIFIER
				  SUBJECT_TERM
				  CHILD_TERM
				  OBJECT_TERM
				  PARENT_TERM
				  PREDICATE_TERM
				  RELATIONSHIP_TYPE
				  ONTOLOGY)
			      ], @args );
   
    $self->init(); 
    
    $self->identifier( $identifier );
    $subject_term = $child unless $subject_term;
    $object_term = $parent unless $object_term;
    $predicate_term = $reltype unless $predicate_term;
    $self->subject_term( $subject_term) if $subject_term;
    $self->object_term( $object_term) if $object_term;
    $self->predicate_term( $predicate_term ) if $predicate_term;
    $self->ontology($ont) if $ont;
                                                    
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
    
    $self->{ "_identifier" }     = undef;
    $self->{ "_subject_term" }   = undef;
    $self->{ "_object_term" }    = undef;
    $self->{ "_predicate_term" } = undef;
    $self->ontology(undef);
   
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




=head2 subject_term

 Title   : subject_term
 Usage   : $rel->subject_term( $subject );
           or
           $subject = $rel->subject_term();
 Function: Set/get for the subject term of this Relationship.

           The common convention for ontologies is to express
           relationships between terms as triples (subject, predicate,
           object).

 Returns : The subject term [Bio::Ontology::TermI].
 Args    : The subject term [Bio::Ontology::TermI] (optional).

=cut

sub subject_term {
    my ( $self, $term ) = @_;
  
    if ( defined $term ) {
        $self->_check_class( $term, "Bio::Ontology::TermI" );
        $self->{ "_subject_term" } = $term;
    }

    return $self->{ "_subject_term" };
    
} # subject_term



=head2 object_term

 Title   : object_term
 Usage   : $rel->object_term( $object );
           or
           $object = $rel->object_term();
 Function: Set/get for the object term of this Relationship.

           The common convention for ontologies is to express
           relationships between terms as triples (subject, predicate,
           object).

 Returns : The object term [Bio::Ontology::TermI].
 Args    : The object term [Bio::Ontology::TermI] (optional).

=cut

sub object_term {
    my ( $self, $term ) = @_;
  
    if ( defined $term ) {
        $self->_check_class( $term, "Bio::Ontology::TermI" );
        $self->{ "_object_term" } = $term;
    }

    return $self->{ "_object_term" };
}



=head2 predicate_term

 Title   : predicate_term
 Usage   : $rel->predicate_term( $type );
           or
           $type = $rel->predicate_term();
 Function: Set/get for the predicate (relationship type) of this
           relationship.

           The common convention for ontologies is to express
           relationships between terms as triples (subject, predicate,
           object).

 Returns : The predicate term [Bio::Ontology::TermI].
 Args    : The predicate term [Bio::Ontology::TermI] (optional).

=cut

sub predicate_term {
    my ( $self, $term ) = @_;
  
    if ( defined $term ) {
        $self->_check_class( $term, "Bio::Ontology::TermI" );
        $self->{ "_predicate_term" } = $term;
    }

    return $self->{ "_predicate_term" };
}


=head2 ontology

 Title   : ontology
 Usage   : $ont = $obj->ontology()
 Function: Get/set the ontology that defined this relationship.
 Example : 
 Returns : an object implementing L<Bio::Ontology::OntologyI>
 Args    : on set, undef or an object implementing 
           Bio::Ontology::OntologyI (optional)

See L<Bio::Ontology::OntologyI>.

=cut

sub ontology{
    my $self = shift;
    my $ont;

    if(@_) {
	$ont = shift;
	if($ont) {
	    $ont = Bio::Ontology::Ontology->new(-name => $ont) if ! ref($ont);
	    if(! $ont->isa("Bio::Ontology::OntologyI")) {
		$self->throw(ref($ont)." does not implement ".
			     "Bio::Ontology::OntologyI. Bummer.");
	    }
	} 
	return $self->{"_ontology"} = $ont;
    } 
    return $self->{"_ontology"};
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
    $s .= "-- Subject Term Identifier:\n";
    $s .= $self->subject_term()->identifier()."\n";
    $s .= "-- Object Term Identifier:\n";
    $s .= $self->object_term()->identifier()."\n";
    $s .= "-- Relationship Type Identifier:\n";
    $s .= $self->predicate_term()->identifier();
    
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

#################################################################
# aliases for backwards compatibility
#################################################################

=head1 Deprecated Methods

  These methods are deprecated and defined here solely to preserve
  backwards compatibility.

=cut

*child_term        = \&subject_term;
*parent_term       = \&object_term;
*relationship_type = \&predicate_term;

1;
