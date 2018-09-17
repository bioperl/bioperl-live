#
# BioPerl module for Bio::Das::FeatureTypeI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Lincoln Stein <lstein@cshl.org>
#
# Copyright Lincoln Stein
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Das::FeatureTypeI - Simple interface to Sequence Ontology feature types

=head1 SYNOPSIS

  # Get a Bio::Das::FeatureTypeI object from somewhere
  $term = $db->fetch....

  # Get the name of the term
  $definition = $term->name;

  # Get the accession of the term
  $accession = $term->accession;

  # Get the definition of the term
  $definition = $term->definition;

  # Get the parents of the term, optionally filtered by relationship
  @parents = $term->parents($relationship);

  # Get the children of the term, optionally filtered by relationship
  @children = $term->children($relationship);

  # Given a parent and child, returns their relationship, or undef if
  # not directly related
  $relationship = $parent->relationship($child);

  # Return true if two terms are identical
  $match = $term1->equals($term2);

  # Return true if $term2 is a descendent of $term1, optionally
  # filtering by relationship ("isa" assumed)
  $match = $term1->is_descendent($term2,$relationship);

  # Return true if $term2 is a parent of $term1, optionally
  # filtering by relationship ("isa" assumed)
  $match = $term1->is_parent($term2,$relationship);

  # Return true if $term2 is equal to $term1 or if $term2 descends
  # from term 1 via the "isa" relationship
  $match = $term1->match($term2);

  # Create a new term de novo
  $term = Bio::Das::FeatureTypeI->new(-name       => $name,
                                      -accession  => $accession,
                                      -definition => $definition);

  # Add a child to a term
  $term1->add_child($term2,$relationship);

  # Delete a child from a term
  $term1->delete_child($term2);

=head1 DESCRIPTION

Bio::Das::FeatureTypeI is an interface to the Gene Ontology
Consortium's Sequence Ontology (SO).  The SO, like other ontologies,
is a directed acyclic graph in which a child node may have multiple
parents.  The relationship between parent and child is one of a list
of relationships.  The SO currently recognizes two relationships "isa"
and "partof".

The intent of this interface is to interoperate with older software
that uses bare strings to represent feature types.  For this reason,
the interface overloads the stringify ("") and string equals (eq)
operations.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bio.perl.org

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Lincoln Stein

Email lstein@cshl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...

package Bio::Das::FeatureTypeI;
use strict;

use overload '""'     => 'name',
             eq       => 'match',
             fallback => 1;

# Object preamble - inherits from Bio::Root::RootI;

=pod

this is somehow FUBAR, implementation classes cannot successfully inherit from Bio::Das::FeatureTypeI

=cut

use base qw(Bio::Root::RootI);

=head2 name

 Title   : name
 Usage   : $string = $term->name
 Function: return the term for the type
 Returns : a string
 Args    : none
 Status  : Public

=cut

sub name { shift->throw_not_implemented }

=head2 accession

 Title   : accession
 Usage   : $string = $term->accession
 Function: return the accession number for the term
 Returns : a string
 Args    : none
 Status  : Public

=cut

sub accession  { shift->throw_not_implemented }

=head2 definition

 Title   : definition
 Usage   : $string = $term->definition
 Function: return the human-readable definition for the term
 Returns : a string
 Args    : none
 Status  : Public

=cut

sub definition  { shift->throw_not_implemented  }

=head2 parents

 Title   : parents
 Usage   : @terms = $term->parents($relationship)
 Function: return parent terms
 Returns : list of Bio::Das::FeatureTypeI
 Args    : none
 Status  : Public

Returns the parents for the current term, empty if there are none.  An
optional relationship argument will return those parents
that are related via the specified relationship type.

The relationship is one of "isa" or "partof".

=cut

sub parents { shift->throw_not_implemented; }

=head2 children

 Title   : children
 Usage   : @terms = $term->children($relationship)
 Function: return children terms
 Returns : list of Bio::Das::FeatureTypeI
 Args    : none
 Status  : Public

Returns the children for the current term, empty if there are none.  An
optional relationship argument will return those children
that are related via the specified relationship type.

The relationship is one of "isa" or "partof".

=cut

sub children { shift->throw_not_implemented; }

=head2 relationship

 Title   : relationship
 Usage   : $relationship = $parent->relationship($child)
 Function: return the relationship between a parent and a child
 Returns : one of "isa" or "partof"
 Args    : none
 Status  : Public

This method returns the relationship between a parent and one of its
immediate descendents.  It can return "isa", "partof", or undef if
there is not a direct parent/child relationship (kissing cousins are
*not* recognized).

=cut

sub relationship { shift->throw_not_implemented }

=head2 equals

 Title   : equals
 Usage   : $boolean = $term1->equals($term2)
 Function: return true if $term1 and $term2 are the same
 Returns : boolean
 Args    : second term
 Status  : Public

The two terms must be identical.  In practice, this means that if
term2 is a Bio::Das::FeatureI object, then its accession number must
match the first term's accession number.  Otherwise, if term2 is a
bare string, then it must equal (in a case insensitive manner)
the name of term1.

NOTE TO IMPLEMENTORS: This method is defined in terms of other
methods, so does not need to be implemented.

=cut

#'
sub equals {
  my $self = shift;
  my $term2 = shift;
  if ($term2->isa('Bio::Das::FeatureTypeI')) {
    return $self->accession eq $term2->accession;
  } else {
    return lc $self->name eq lc $term2;
  }
}

=head2 is_descendent

 Title   : is_descendent
 Usage   : $boolean = $term1->is_descendent($term2 [,$relationship])
 Function: return true of $term2 is a descendent of $term1
 Returns : boolean
 Args    : second term
 Status  : Public

This method returns true if $term2 descends from $term1.  The
operation traverses the tree.  The traversal can be limited to the
relationship type ("isa" or "partof") if desired.  $term2 can be a
bare string, in which case the term names will be used as the basis
for term matching (see equals()).

NOTE TO IMPLEMENTORS: this method is defined as the inverse of
is_parent().  Do not implement it directly, but do implement
is_parent().

=cut

sub is_descendent {
  my $self = shift;
  my ($term,$relationship) = @_;
  $self->throw("$term is not a Bio::Das::FeatureTypeI")
    unless $term->isa('Bio::Das::FeatureTypeI');
  $term->is_parent($self,$relationship);
}

=head2 is_parent

 Title   : is_parent
 Usage   : $boolean = $term1->is_parent($term2 [,$relationship])
 Function: return true of $term2 is a parent of $term1
 Returns : boolean
 Args    : second term
 Status  : Public

This method returns true if $term2 is a parent of $term1.  The
operation traverses the tree.  The traversal can be limited to the
relationship type ("isa" or "partof") if desired.  $term2 can be a
bare string, in which case the term names will be used as the basis
for term matching (see equals()).

NOTE TO IMPLEMENTORS: Implementing this method will also implement
is_descendent().

=cut

sub is_parent { shift->throw_not_implemented }

=head2 match

 Title   : match
 Usage   : $boolean = $term1->match($term2)
 Function: return true if $term1 equals $term2 or if $term2 is an "isa" descendent
 Returns : boolean
 Args    : second term
 Status  : Public

This method combines equals() and is_descendent() in such a way that
the two terms will match if they are the same or if the second term is
an instance of the first one.  This is also the basis of the operator
overloading of eq.

NOTE TO IMPLEMENTORS: This method is defined in terms of other methods
and does not need to be implemented.

=cut

sub match {
  my $self  = shift;
  my $term2 = shift;
  return 1 if $self->equals($term2);
  return $self->is_descendent($term2,'isa');
}

=head2 new

 Title   : new
 Usage   : $term = Bio::Das::FeatureTypeI->new(@args)
 Function: create a new term
 Returns : new term
 Args    : see below
 Status  : Public

This method creates a new Bio::Das::FeatureTypeI.  Arguments:

  Argument    Description
  --------   ------------

  -name       Name of this term

  -accession  Accession number for the term

  -definition Definition of the term

=cut

sub new { shift->throw_not_implemented }

=head2 add_child

 Title   : add_child
 Usage   : $boolean = $term->add_child($term2,$relationship)
 Function: add a child to a term
 Returns : a boolean indicating success
 Args    : new child
 Throws  : a "cycle detected" exception
 Status  : Public

This method adds a new child to the indicated node.  It may detect a
cycle in the DAG and throw a "cycle detected" exception.

=cut

sub add_child { shift->throw_not_implemented }


=head2 delete_child

 Title   : delete_child
 Usage   : $boolean = $term->delete_child($term2);
 Function: delete a child of the term
 Returns : a boolean indicating success
 Args    : child to be deleted
 Throws  : a "not a child" exception
 Status  : Public

This method deletes a new child from the indicated node.  It will
throw an exception if the indicated child is not a direct descendent.

=cut

sub delete_child { shift->throw_not_implemented }

1;
