#
# BioPerl module for Bio::DB::Tree::AnnotatedImpl
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Tree::AnnotatedImpl - Methods to implement get/set methods for tag/value pairs for Database backed Tree and Node objects

=head1 SYNOPSIS

DO NOT USE THIS MODULE DIRECTLY

=head1 DESCRIPTION

This module is just a collection of implementations for
Bio::DB::Tree::Tree and Bio::DB::Tree::Node objects for storing
annotations in the form of key/value pairs.


=cut

package Bio::DB::Tree::AnnotatedImpl;
use strict;
use base qw(Bio::Root::RootI);



=head2 set_flatAnnotations

 Title   : set_flatAnnotations
 Usage   : $tree->_setAnnotations($flatannotations)
 Function: Directly set the tag/values hash for the object based on a FlatAnnotation format (see L<Bio::DB::Tree::Store::DBI::SQLite> for example)
 Returns : none
 Args    : flatannotation in the form of "key1=value1,value2;key2=value3"

=cut


sub set_flatAnnotations {
  my $self = shift;
  return unless @_;
  my $flatannot = shift;
  return unless defined $flatannot;
  # should this be implemented in the SQLite module?
  for my $pair ( split(';',$flatannot) ) {
    if ( $pair !~ /=/ ) {
      $self->warn("improperly formed flat annotation field $pair\n");
      next;
    }
    my ($key,$values) = split('=',$pair,2);
    $self->{'_tags'}->{$key} = [split(',',$values)];
    # perf test here to see if there is advantage in direct set vs calling the method
    # $self->set_tag_values($key,split(',',$values));
  }
  return;
}

=head2 Get/Set Tag/Value pairs

These methods associate tag/value pairs with a Tree

=head2 set_tag_value

 Title   : set_tag_value
 Usage   : $tree->set_tag_value($tag,$value)
           $tree->set_tag_value($tag,@values)
 Function: Sets a tag value(s) to a tree. Replaces old values.
 Returns : number of values stored for this tag
 Args    : $tag   - tag name
           $value - value to store for the tag

=cut

sub set_tag_value{
    my ($self,$tag,@values) = @_;
    if( ! defined $tag || ! scalar @values ) {
      $self->warn("cannot call set_tag_value with an undefined value");
    }
    $self->remove_tag ($tag);
    map { push @{$self->{'_tags'}->{$tag}}, $_ } @values;
    $self->_dirty(1);
    return scalar @{$self->{'_tags'}->{$tag}};
}

=head2 add_tag_value

 Title   : add_tag_value
 Usage   : $tree->add_tag_value($tag,$value)
 Function: Adds a tag value to a tree 
 Returns : number of values stored for this tag
 Args    : $tag   - tag name
           $value - value to store for the tag

=cut

sub add_tag_value{
    my ($self,$tag,$value) = @_;
    if( ! defined $tag || ! defined $value ) {
      $self->warn("cannot call add_tag_value with an undefined value");
    }
    push @{$self->{'_tags'}->{$tag}}, $value;
    $self->_dirty(1);
    return scalar @{$self->{'_tags'}->{$tag}};
}

=head2 remove_tag

 Title   : remove_tag
 Usage   : $tree->remove_tag($tag)
 Function: Remove the tag and all values for this tag
 Returns : boolean representing success (0 if tag does not exist)
 Args    : $tag - tagname to remove


=cut

sub remove_tag {
   my ($self,$tag) = @_;
   if( exists $self->{'_tags'}->{$tag} ) {
     $self->{'_tags'}->{$tag} = undef;
     delete $self->{'_tags'}->{$tag};
     $self->_dirty(1);
     return 1;
   }
   return 0;
}

=head2 remove_all_tags

 Title   : remove_all_tags
 Usage   : $tree->remove_all_tags()
 Function: Removes all tags 
 Returns : None
 Args    : None

=cut

sub remove_all_tags{
   my ($self) = @_;
   $self->{'_tags'} = {};
   $self->_dirty(1);
   return;
}

=head2 get_all_tags

 Title   : get_all_tags
 Usage   : my @tags = $tree->get_all_tags()
 Function: Gets all the tag names for this Tree
 Returns : Array of tagnames
 Args    : None

=cut

sub get_all_tags{
   my ($self) = @_;
   if ( ! exists $self->{'_tags'} ) {
     $self->_load_from_db;
   }
   my @tags = sort keys %{$self->{'_tags'} || {}};
   return @tags;
}

=head2 get_tag_values

 Title   : get_tag_values
 Usage   : my @values = $tree->get_tag_values($tag)
 Function: Gets the values for given tag ($tag)
 Returns : Array of values or empty list if tag does not exist
 Args    : $tag - tag name

=cut

sub get_tag_values{
   my ($self,$tag) = @_;
   if ( ! exists $self->{'_tags'} ) {
     $self->_load_from_db;
   }
   return wantarray ? @{$self->{'_tags'}->{$tag} || []} :
                     (@{$self->{'_tags'}->{$tag} || []})[0];
}

=head2 has_tag

 Title   : has_tag
 Usage   : $tree->has_tag($tag)
 Function: Boolean test if tag exists in the Tree
 Returns : Boolean
 Args    : $tag - tagname

=cut

sub has_tag {
   my ($self,$tag) = @_;
   if ( ! exists $self->{'_tags'} ) {
     $self->_load_from_db;
   }
   return exists $self->{'_tags'}->{$tag};
}

1;
