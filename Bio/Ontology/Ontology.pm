# $Id$
#
# BioPerl module for Bio::Ontology::Ontology
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

    # create ontology object
    my $ont = Bio::Ontology::Ontology->new(-name => "OBF");

    # add terms, relationships ...
    my $bp = Bio::Ontology::Term->new(-name => "Bioperl");
    my $obf = Bio::Ontology::Term->new(-name => "OBF");
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

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Hilmar Lapp

Email hlapp at gmx.net

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Ontology::Ontology;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::Ontology::OntologyI;
use Bio::Ontology::SimpleOntologyEngine;

@ISA = qw(Bio::Root::Root Bio::Ontology::OntologyI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Ontology::Ontology();
 Function: Builds a new Bio::Ontology::Ontology object 
 Returns : an instance of Bio::Ontology::Ontology
 Args    :


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
    $engine = Bio::Ontology::SimpleOntologyEngine->new() unless $engine;
    $self->engine($engine);

    return $self;
}

=head1

  Methods from L<Bio::Ontology::OntologyI>

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

=head1

  Implementation-specific public methods

=cut

=head2 engine

 Title   : engine
 Usage   : $engine = $obj->engine()
 Function: Get/set the ontology engine to which all the query methods
           delegate.
 Example : 
 Returns : an object implementing L<Bio::Ontology::OntologyEngineI>
 Args    : on set, new value (an object implementing 
           L<Bio::Ontology::OntologyEngineI>, or  undef)


=cut

sub engine{
    my $self = shift;
    
    if(@_) {
	my $engine = shift;
	if($engine && (! (ref($engine) &&
			  $engine->isa("Bio::Ontology::OntologyEngineI")))) {
	    $self->throw("object of class ".ref($engine)." does not implement".
			 " Bio::Ontology::OntologyEngineI. Bummer!");
	}
	$self->{'engine'} = $engine;
    }
    return $self->{'engine'};
}

=head1

  Methods defined in L<Bio::Ontology::OntologyEngineI>

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
  add_relatioship(TermI parent, TermI child, TermI relationship_type)
 Function: Adds a relationship object to the ontology engine.
 Example :
 Returns : Its argument.
 Args    : A RelationshipI object.


=cut

sub add_relationship{
    my $self = shift;
    my $rel = shift;

    # set ontology if not set already, and if it's a RelationshipI object
    $rel->ontology($self)
	if( $rel && ref($rel) &&
	    $rel->isa("Bio::Ontology::RelationshipI") && (! $rel->ontology()));
    return $self->engine->add_relationship($rel,@_);
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
    return shift->engine->get_relationships(@_);
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
    return shift->engine->get_relationship_types(@_);
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
    return shift->engine->get_child_terms(@_);
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
    return shift->engine->get_descendant_terms(@_);
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
    return shift->engine->get_parent_terms(@_);
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
    return shift->engine->get_ancestor_terms(@_);
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
    return shift->engine->get_leaf_terms(@_);
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

sub get_root_terms{
    return shift->engine->get_root_terms(@_);
}

1;
