# $Id $
# A simple implementation of the Bio::SeqFeature::TypeI interface.

=head1 NAME

Bio::SeqFeature::SimpleType - A feature type.

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

  SimpleType is a trivial implementation of TypeI.  Parents and
  children are stored in lists, so is_parent(..) and is_child(..) take
  some time.  Ancestors and descendents are determined by a real-time
  breadth-first traversal of the ISA graph starting at this node.  The
  graph is represented only by the '_parents' and '_children' lists in
  the SimpleType nodes participating, although mixing with other
  TypeIs is allowed (subject to some constraints; see the add_child()
  and _add_parents() methods).  The graph is alterable only via the
  public add_child() and remove_child() methods of participants.

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

Robert Hubley E<lt>rhubley@systemsbiology.orgE<gt>.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqFeature::SimpleType;
use vars qw( $VERSION @ISA );
use strict;
$VERSION = '0.01';

use Bio::DB::GFF::Util::Rearrange; # for 'rearrange'

use Bio::Root::Root;
use Bio::SeqFeature::TypeI;

@ISA = qw( Bio::Root::Root Bio::SeqFeature::TypeI );

=head2 new

 Title   : new
 Usage   : my $simple_type =
              new Bio::SeqFeature:SimpleType( @args );
 Function: Builds a new Bio::SeqFeature::SimpleType object·
 Returns : The new L<Bio::SeqFeature::SimpleType> object
 Args    : See below
 Status  : Public

  Argument       Description
  --------       ------------

  -unique_id     The unique identifier of this Type

  -name          The (not necessarily unique) name of this Type

  -description   The description of this Type, aka the definition of this Term

  -category      The category of this Term (see L<Bio::Ontology::TermI>)

  -version       The version of this Term (see L<Bio::Ontology::TermI>)

  -comment       Some additional commentary on this Type/Term

  -synonyms      An array ref of synonyms of this term

=cut

sub new {
  my $caller = shift;
  my $self = $caller->SUPER::new( @_ );

  my ( $unique_id, $name, $description, $category, $version, $comment, $synonyms );
  if( scalar( @_ ) && $_[ 0 ] =~ /^-/ ) {
    ( $unique_id, $name, $description, $category, $version, $comment, $synonyms ) =
      rearrange(
        [ [ qw( IDENTIFIER ID UNIQUE_ID UNIQUEID ) ],
          qw( NAME ),
          [ qw( DESCRIPTION DEFINITION ) ],
          qw( CATEGORY ),
          [ qw( VERSION V ) ],
          qw( COMMENT ),
          [ qw( SYNONYM SYNONYMS AKA ) ]
        ],
        @_
      );
  } else {
    ## TODO: Do we accept this form of argument list?
    if( $_[ 0 ] || ( scalar( @_ ) > 1 ) ) {
      $self->warn( "arguments to the new() method of Bio::SeqFeature::SimpleType have been ignored: Bio::SeqFeature::SimpleType->new( " . join( ", ", @_ ) . " )" );
    }
  }

  $unique_id && $self->unique_id( $unique_id );
  $name && $self->name( $name );
  $description && $self->description( $description );
  $category && $self->category( $category );
  $version && $self->version( $version );
  $comment && $self->comment( $comment );

  if( defined $synonyms ) {
    if( ref $synonyms eq 'ARRAY' ) {
      $self->add_synonyms( @$synonyms );
    } else {
      $self->add_synonyms( $synonyms );
    }
  }

  return $self;
} # new(..)

=head2 unique_id

 Title   : unique_id
 Usage   : my $unique_id = $self->unique_id( [$new_id] );
 Function: Get/Set the B<unique> identifier of this Type
 Returns : The current (or former, if used as a set method) guaranteed-to-be
           unique identifier or undef if there is none.
 Args    : [optional] A new identifier
 Status  : Public

=cut

sub unique_id {
  my ( $self, $identifier ) = @_;
  my $current_val = $self->{ '_unique_id' };
  if( defined $identifier ) {
    $self->{ '_unique_id' } = $identifier;
  }
  return $current_val;
} # unique_id(..)

=head2 identifier

 Title   : identifier
 Usage   : my $unique_id = $self->idenifier( [$new_id] );
 Function: Get/Set the B<unique> identifier of this Type
 Returns : A guaranteed-to-be unique identifier or undef if there is none.
 Args    : [optional] A new identifier
 Status  : Public

   This is implemented as a (glob ref) alias for unique_id(..)

=cut

  *identifier = \&unique_id;

=head2 name

 Title   : name
 Usage   : my $name= $self->name( [$new_name] );
 Function: Get/Set the name of this Type
 Returns : The current (or former, if used as a set method) value of the name.
 Args    : [optional] A new name
 Status  : Public

=cut

sub name {
  my ( $self, $name ) = @_;
  my $current_val = $self->{'_name'};
  if( defined $name ) { 
      $self->{'_name'} = $name;
  }
  return $current_val;
} # name(..)

=head2 description

 Title   : description
 Usage   : my $description = $self->description( [$new_description] );
 Function: Get/Set the description of this Type
 Returns : The current (or former, if used as a set method) value of the
           description.
 Args    : [optional] A new description
 Status  : Public

=cut

sub description {
  my ( $self, $description ) = @_;
  my $current_val = $self->{ '_description' };
  if( defined $description ) {
    $self->{ '_description' } = $description;
  }
  return $current_val;
} # description(..)

=head2 definition

 Title   : definition
 Usage   : my $description = $self->definition( [$new_description] );
 Function: Get/Set the description of this Type
 Returns : The current (or former, if used as a set method) value of the
           description.
 Args    : [optional] A new description
 Status  : Public

  This is implemented as a (glob ref) alias for description(..).

=cut

  *definition = \&description;

=head2 category

 Title   : category
 Usage   : my $category = $self->category( [$new_category] );
 Function: Get/Set the category of this Type
 Returns : The current (or former, if used as a set method) value of the
           category.
 Args    : [optional] A new category
 Status  : Public

=cut

sub category {
  my ( $self, $category ) = @_;
  my $current_val = $self->{ '_category' };
  if( defined $category ) {
    $self->{ '_category' } = $category;
  }
  return $current_val;
} # category(..)

=head2 version

 Title   : version
 Usage   : my $version = $self->version( [$new_version] );
 Function: Get/Set the version of this Type
 Returns : The current (or former, if used as a set method) value of the
           version.
 Args    : [optional] A new version
 Status  : Public

=cut

sub version {
  my ( $self, $version ) = @_;
  my $current_val = $self->{ '_version' };
  if( defined $version ) {
    $self->{ '_version' } = $version;
  }
  return $current_val;
} # version(..)

=head2 comment

 Title   : comment
 Usage   : my $comment = $self->comment( [$new_comment] );
 Function: Get/Set the comment of this Type
 Returns : The current (or former, if used as a set method) value of the
           comment.
 Args    : [optional] A new comment
 Status  : Public

=cut

sub comment {
  my ( $self, $comment ) = @_;
  my $current_val = $self->{ '_comment' };
  if( defined $comment ) {
    $self->{ '_comment' } = $comment;
  }
  return $current_val;
} # comment(..)

=head2 each_synonym

 Title   : each_synonym
 Usage   : my @aliases = $simple_type->each_synonym();
 Function: Returns a list of aliases of this SimpleType.
 Returns : A list of aliases, which may be any type but are usually strings,
           or undef if there are none.
 Args    : none
 Status  : Public

=cut

sub each_synonym {
  my $self = shift;

  if( $self->{ '_synonyms' } ) {
    return @{ $self->{ '_synonyms' } };
  } else {
    return undef;
  }
} # each_synonym()

=head2 add_synonyms

 Title   : add_synonyms
 Usage   : $simple_type->add_synonyms( @synonyms );
 Function: Adds (pushes) one or more synonyms into the list of synonyms.
 Returns : Nothing
 Args    : A list of new synonyms (synonyms may be any scalar but are generally
           strings)
 Status  : Public

  Note that there is no attempt to ensure that a synonym appears only once.

=cut

sub add_synonyms {
  my $self = shift;
  push( @{ $self->{ '_synonyms' } }, @_ );
} # add_synonyms

=head2 remove_synonyms

 Title   : remove_synonyms
 Usage   : $simple_type->remove_synonyms();
 Function: Deletes (and returns) the synonyms of this SimpleType.
 Returns : The former list of synonyms, as with each_synonym().
 Args    : none
 Status  : Public

=cut

sub remove_synonyms {
    my ( $self ) = @_;

    my @a = $self->each_synonym();
    $self->{ "_synonyms" } = [];
    return @a;

} # remove_synonyms

=head2 add_child

 Title   : add_child
 Usage   : $simple_type->add_child( $another_simple_type );
 Function: Add a child type to the children list
 Returns : Nothing
 Args    : A L<Bio::SeqFeature::SimpleType> object (see below)
 Status  : Public

  Note that the new child must be a SimpleType object, or at least
  implement the _add_parent(..) method, to maintain consistency.
  Otherwise an exception will be thrown.

=cut

sub add_child {
  my $self = shift;
  my $new_child = shift;

  unless( ref( $new_child ) && $new_child->can( '_add_parent' ) ) {
    $self->throw( "Unable to add a child type that is not a SimpleType object" );
  }
  $new_child->_add_parent( $self );
  push( @{ $self->{ '_children' } }, $new_child );
} # add_child(..)

=head2 _add_parent

 Title   : _add_parent
 Usage   : $simple_type->_add_parent( $another_simple_type );
 Function: Add a parent type to the parent list
 Returns : Nothing
 Args    : A L<Bio::SeqFeature::SimpleType> object
 Status  : Protected (to be called by other SimpleType objects only)

  Note that the new parent must be a SimpleType object, to maintain
  consistency.  Otherwise an exception will be thrown.

=cut

sub _add_parent {
  my $self = shift;
  my $new_parent = shift;

  unless( ref( $new_parent ) &&
          $new_parent->isa( 'Bio::SeqFeature::SimpleType' ) ) {
    $self->throw( "Unable to add a parent type that is not a SimpleType object" );
  }
  push( @{ $self->{ '_parents' } }, $new_parent );
} # _add_parent(..)

=head2 remove_child

 Title   : remove_child
 Usage   : $simple_type->remove_child( $another_simple_type );
 Function: Remove a child type from the children list
 Returns : The removed child, iff the child was previously present and
           is now removed.
 Args    : A L<Bio::SeqFeature::SimpleType> object
 Status  : Public

=cut

sub remove_child {
  my $self = shift;
  my $child_to_remove = shift;

  my $child;
  for( my $i = 0; $i < scalar( @{ $self->{ '_children' } } ); $i++ ) {
    $child = $self->{ '_children' }->[ $i ];
    if( $child eq $child_to_remove ) {
      # Note that we silently ignore a false result.  This is laziness.
      $child->_remove_parent( $self );
      splice( @{ $self->{ '_children' } }, $i, 1 );
      return $child;
    }
  }
  return undef;
} # remove_child(..)

=head2 _remove_parent

 Title   : _remove_parent
 Usage   : $simple_type->_remove_parent( $another_simple_type );
 Function: Remove a parent type from the parent list
 Returns : The removed parent, iff the parent was previously present
           and is now removed.
 Args    : A L<Bio::SeqFeature::SimpleType> object
 Status  : Protected (to be called by other SimpleType objects only)

=cut

sub _remove_parent {
  my $self = shift;
  my $parent_to_remove = shift;

  my $parent;
  for( my $i = 0; $i < scalar( @{ $self->{ '_parents' } } ); $i++ ) {
    $parent = $self->{ '_parents' }->[ $i ];
    if( $parent eq $parent_to_remove ) {
      splice( @{ $self->{ '_parents' } }, $i, 1 );
      return $parent;
    }
  }
  return undef;
} # _remove_parent(..)

=head2 parents

 Title   : parents
 Usage   : @parents = $feature_type->parents()
 Function: retrieves the parents of this feature type (ISA relationship parents)
 Returns : a list of Bio::SeqFeature::SimpleType objects, or undef if this is a
           root
 Args    : none
 Status  : Public

=cut

sub parents {
  my $self = shift;
  unless( $self->{ '_parents' } && @{ $self->{ '_parents' } } ) {
    return undef;
  }
  return @{$self->{ '_parents' }};
} # parents()

=head2 children

 Title   : children
 Usage   : @children = $feature_type->children()
 Function: retrieves the children of this feature type
           (ISA relationship children)
 Returns : a list of Bio::SeqFeature::SimpleType objects, or undef if
           this is a leaf
 Args    : none
 Status  : Public

=cut

sub children {
  my $self = shift;
  unless( $self->{ '_children' } && @{ $self->{ '_children' } } ) {
    return undef;
  }
  return @{$self->{ '_children' }};
} # children()

=head2 ancestors

 Title   : ancestors
 Usage   : @ancestors = $feature_type->ancestors()
 Function: retrieves the ancestors of this feature type
           (ISA relationship ancestors)
 Returns : a list of Bio::SeqFeature::SimpleType objects, or undef if
           this is a root; the list will contain the ancestors in a
           breadth first traversal ordering (no element will appear
           after its ancestor)
 Args    : none
 Status  : Public

=cut

sub ancestors {
  my $self = shift;
  my ( @ancestors, @queue, %seen );

  push( @queue, $self );
  $seen{ $self } = 1;

  my $type_object;
  while ( @queue ) {
    $type_object = shift @queue; # dequeue
    foreach my $type_objects_parent ( $type_object->parents() ) {
      next unless defined( $type_objects_parent );
      next if $seen{ $type_objects_parent };
      $seen{ $type_objects_parent } = 1;
      push( @ancestors, $type_objects_parent ); # save for returning
      push( @queue, $type_objects_parent ); # enqueue
    }
  }
  if( @ancestors ) {
    return @ancestors;
  } else {
    return undef;
  }
} # ancestors()

=head2 descendents

 Title   : descendents
 Usage   : @descendents = $feature_type->descendents()
 Function: retrieves the descendents of this feature type
           (ISA relationship descendents)
 Returns : a list of Bio::SeqFeature::SimpleType objects, or undef if this is a
           leaf; the list will contain the descendents in a breadth first
           traversal ordering (no element will appear after its descendent)
 Args    : none
 Status  : Public

=cut

sub descendents {
  my $self = shift;
  my ( @descendents, @queue, %seen );

  push( @queue, $self );
  $seen{ $self } = 1;

  my $type_object;
  while ( @queue ) {
    $type_object = shift @queue; # dequeue
    foreach my $type_objects_child ( $type_object->children() ) {
      next unless defined( $type_objects_child );
      next if $seen{ $type_objects_child };
      $seen{ $type_objects_child } = 1;
      push( @descendents, $type_objects_child ); # save for returning
      push( @queue, $type_objects_child ); # enqueue
    }
  }
  if( @descendents ) {
    return @descendents;
  } else {
    return undef;
  }
} # descendents()

=head2 is_parent

 Title   : is_parent
 Usage   : if( $feature_type->is_parent( $potential_parent_type ) ) { .. }
 Function: Accesses the relationship between this type and another
 Returns : true iff the given type is a parent of this one
 Args    : a Bio::SeqFeature::TypeI object, or a (string) unique_id value
 Status  : Public

=cut

sub is_parent {
  my $self = shift;
  my $query_type_object = shift;
  foreach my $parent ( @{ $self->{ '_parents' } } ) {
    if( $parent eq $query_type_object ) {
      return 1;
    }
  }
  return 0;
} # is_parent(..)

=head2 is_ancestor

 Title   : is_ancestor
 Usage   : if( $feature_type->is_ancestor( $potential_ancestor_type ) ) { .. }
 Function: Accesses the relationship between this type and another
 Returns : true iff the given type is an ancestor of this one
 Args    : a Bio::SeqFeature::TypeI object, or a (string) unique_id value
 Status  : Public

=cut

sub is_ancestor {
  my $self = shift;
  my $query_type_object = shift;

  my ( @queue, %seen );

  push( @queue, $self );
  $seen{ $self } = 1;

  my $type_object;
  while ( @queue ) {
    $type_object = shift @queue; # dequeue
    foreach my $type_objects_parent ( $type_object->parents() ) {
      next unless defined( $type_objects_parent );
      next if $seen{ $type_objects_parent };
      if( $type_objects_parent eq $query_type_object ) {
        return 1;
      }
      $seen{ $type_objects_parent } = 1;
      push( @queue, $type_objects_parent ); # enqueue
    }
  }
  return 0;
} # is_ancestor(..)

=head2 is_child

 Title   : is_child
 Usage   : if( $feature_type->is_child( $potential_child_type ) ) { .. }
 Function: Accesses the relationship between this type and another
 Returns : true iff the given type is a child of this one
 Args    : a Bio::SeqFeature::TypeI object, or a (string) unique_id value
 Status  : Public

=cut

sub is_child {
  my $self = shift;
  my $query_type_object = shift;

  foreach my $child ( @{$self->{'_children'}} ) {
    if( $child eq $query_type_object ) {
      return 1;
    }
  }
  return 0;
} # is_child(..)

=head2 is_descendent

 Title   : is_descendent
 Usage   : if( $feature_type->is_descendent( $potential_descendent_type ) ) { .. }
 Function: Accesses the relationship between this type and another
 Returns : true iff the given type is a descendent of this one
 Args    : a Bio::SeqFeature::TypeI object, or a (string) unique_id value
 Status  : Public

=cut

sub is_descendent {
  my $self = shift;
  my $query_type_object = shift;

  my ( @queue, %seen );

  push( @queue, $self );
  $seen{ $self } = 1;

  my $type_object;
  while( @queue ) {
    $type_object = shift @queue; # dequeue
    foreach my $type_objects_child ( $type_object->children() ) {
      next unless defined( $type_objects_child );
      next if $seen{ $type_objects_child };
      $seen{ $type_objects_child } = 1;
      if( $type_objects_child eq $query_type_object ) {
        return 1;
      }
      push( @queue, $type_objects_child ); # enqueue
    }
  }
  return 0;
} # is_descendent(..)

1;

__END__
