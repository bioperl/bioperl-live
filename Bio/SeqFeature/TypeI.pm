# $Id$

=head1 NAME

Bio::SeqFeature::TypeI - A feature type.

=head1 SYNOPSIS

=head1 DESCRIPTION

  A TypeI is the type of a SeqFeature object.  It overloads the
  stringification operator to behave as a string for backward
  compatibility.  It is a Bio::Overload::TermI and participates in a
  directed acyclic graph of ISA relationships (eg. the type 'SINE' isa
  'transposable element').  Note that these relationships are not
  (necessarily) reflected in a Perl hierarchy, and that the ISA
  relationship is not the only possible relationship between types
  (the HASA relationship may also be represented, for example), but is
  the relationship implied in the methods referring to parents,
  children, ancestors, and descendents.  Again: the ISA relating TypeI
  objects has B<nothing> to do with the @ISA relating Perl classes!

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
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR

Paul Edlefsen E<lt>paul@systemsbiology.orgE<gt>.

Copyright (c) 2003 Institute for Systems Biology

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=head1 CONTRIBUTORS

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqFeature::TypeI;

use strict;
use vars qw( $VERSION @ISA );
use overload 
  '""' => 'toString',
  cmp  => '_cmp';

use Bio::Root::RootI;
use Bio::Ontology::TermI;
use Bio::LocallyIdentifiableI;

$VERSION = '0.01';
@ISA = qw( Bio::Root::RootI Bio::Ontology::TermI Bio::LocallyIdentifiableI );

=head2 parents

 Title   : parents
 Usage   : @parents = $feature_type->parents()
 Function: retrieves the parents of this feature type (ISA relationship parents)
 Returns : a list of Bio::SeqFeature::TypeI objects, or undef if this is a root
 Args    : none
 Status  : Public

=cut

sub parents {
  shift->throw_not_implemented( @_ );
}

=head2 children

 Title   : children
 Usage   : @children = $feature_type->children()
 Function: retrieves the children of this feature type
           (ISA relationship children)
 Returns : a list of Bio::SeqFeature::TypeI objects, or undef if this is a leaf
 Args    : none
 Status  : Public

=cut

sub children {
  shift->throw_not_implemented( @_ );
}

=head2 ancestors

 Title   : ancestors
 Usage   : @ancestors = $feature_type->ancestors()
 Function: retrieves the ancestors of this feature type
           (ISA relationship ancestors)
 Returns : a list of Bio::SeqFeature::TypeI objects, or undef if this is a root;
           the list will contain the ancestors in a breadth first traversal
           ordering (no element will appear after its ancestor)
 Args    : none
 Status  : Public

=cut

sub ancestors {
  shift->throw_not_implemented( @_ );
}

=head2 descendents

 Title   : descendents
 Usage   : @descendents = $feature_type->descendents()
 Function: retrieves the descendents of this feature type
           (ISA relationship descendents)
 Returns : a list of Bio::SeqFeature::TypeI objects, or undef if this is a leaf;
           the list will contain the descendents in a breadth first traversal
           ordering (no element will appear after its descendent)
 Args    : none
 Status  : Public

=cut

sub descendents {
  shift->throw_not_implemented( @_ );
}

=head2 is_parent

 Title   : is_parent
 Usage   : if( $feature_type->is_parent( $potential_parent_type ) ) { .. }
 Function: Accesses the relationship between this type and another
 Returns : true iff the given TypeI is a parent of this one
 Args    : a Bio::SeqFeature::TypeI object, or a (string) unique_id value
 Status  : Public

=cut

sub is_parent {
  shift->throw_not_implemented( @_ );
}

=head2 is_ancestor

 Title   : is_ancestor
 Usage   : if( $feature_type->is_ancestor( $potential_ancestor_type ) ) { .. }
 Function: Accesses the relationship between this type and another
 Returns : true iff the given TypeI is an ancestor of this one
 Args    : a Bio::SeqFeature::TypeI object, or a (string) unique_id value
 Status  : Public

=cut

sub is_ancestor {
  shift->throw_not_implemented( @_ );
}

=head2 is_child

 Title   : is_child
 Usage   : if( $feature_type->is_child( $potential_child_type ) ) { .. }
 Function: Accesses the relationship between this type and another
 Returns : true iff the given TypeI is a child of this one
 Args    : a Bio::SeqFeature::TypeI object, or a (string) unique_id value
 Status  : Public

=cut

sub is_child {
  shift->throw_not_implemented( @_ );
}

=head2 is_descendent

 Title   : is_descendent
 Usage   : if( $feature_type->is_descendent( $potential_descendent_type ) ) { .. }
 Function: Accesses the relationship between this type and another
 Returns : true iff the given TypeI is a descendent of this one
 Args    : a Bio::SeqFeature::TypeI object, or a (string) unique_id value
 Status  : Public

=cut

sub is_descendent {
  shift->throw_not_implemented( @_ );
}

=head2 id

 Title   : id
 Usage   : $id = $feature_type->id( [new_id] )
 Function: gets/sets the unique id of this feature type
 Returns : a String
 Args    : a new value (optional)
 Status  : Public

  This is an alias for identifier().

=cut

sub id {
  shift->identifier( @_ );
}

=head2 accession

 Title   : accession
 Usage   : $id = $feature_type->accession( [new_id] )
 Function: gets/sets the unique id of this feature type
 Returns : a String
 Args    : a new value (optional)
 Status  : Public

  This is an alias for identifier().

=cut

sub accession {
  shift->identifier( @_ );
}

=head2 label

 Title   : label
 Usage   : $name = $feature_type->label( [new_name] )
 Function: gets/sets the name of this feature type
 Returns : a String
 Args    : a new value (optional)
 Status  : Public

  This is an alias for name().

=cut

sub label {
  shift->name( @_ );
}

=head2 toString

 Title   : toString
 Usage   : $str_val = $feature_type->toString()
 Function: returns the string value of this feature type
           (name || identifier || StrVal)
 Returns : a String
 Args    : None
 Status  : Public

  This method is defined in the interface and need not be implemented elsewhere.

=cut

sub toString {
  my $self = shift;
  $self->name() || $self->identifier() || overload::StrVal( $self );
} # toString()

## method for overload for comparing two TypeI objects.
## Uses unique_id() || toString().
sub _cmp {
  my $self = shift;
  my ( $b, $reversed ) = @_;
  my ( $a_string, $b_string );
  if( ref( $b ) && $b->isa( 'Bio::SeqFeature::TypeI' ) ) {
    $a_string = $self->unique_id();
    $b_string;
    if( $a_string ) {
      $b_string = $b->unique_id();
    } else {
      $a_string = $self->toString();
      $b_string = $b->toString();
    }
  } else { # $b is not a TypeI.  We'll treat it like a string, then.
    $a_string = $self->unique_id() || $self->toString();
    $b_string = $b;
  }
  ( $a_string, $b_string ) = ( $b_string, $a_string ) if $reversed;
  return ( $a_string cmp $b_string );
} # _cmp(..)

1;
