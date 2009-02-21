# BioPerl module for Bio::Annotation::Tree
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Weigang Qiu <weigang at genectr.hunter.cuny.edu>
#
# Based on the Bio::Annotation::DBLink by Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Annotation::Tree - Provide a tree as an annotation to a Bio::AnnotatableI
object

=head1 SYNOPSIS

   # Read a tree and an alignment

   $treeio=Bio::TreeIO->new(-file=>'foo.dnd', -format=>'newic');
   $tree=$treeio->next_tree;
   $alnio=Bio::AlignIO->new(-file=>'foo.aln', -format=>'clustalw');
   $aln=$alnio->next_aln;

   # Construct a tree annotation
   $ann_tree = Bio::Annotation::Tree->new (-tree_id  => 'mytree',
                                           -tree_obj     => $tree,
                                            );

   # Add the tree annotation to AlignI
   $ac = Bio::Annotation::Collection->new();
   $ac->add_Annotation('tree', $ann_tree);
   $aln->annotation($ac);

   # NOTE & TODO: 
   # The above procedures are sensible only if 
   # the tree is generated from the alignment.  However, 
   # currently no effort has been made to check the consistency
   # between the tree OTU names and the sequence names

=head1 DESCRIPTION

Provides a Bio::AnnotationI object which contains a Bio::Tree::TreeI, which can
be added to a Bio::AnnotationCollectionI, which in turn be attached to a
Bio::AnnotatableI (typically a Bio::AlignI object)

=head1 AUTHOR

Weigang Qiu - weigang at genectr.hunter.cuny.edu

=head1 CONTRIBUTORS

Aaron Mackey
Jason Stajich

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a '_'

=cut

# Let the code begin...

package Bio::Annotation::Tree;
use strict;

use base qw(Bio::Root::Root Bio::AnnotationI Bio::TreeIO);


sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  my ($tree_id, $tree_obj, $tag) =
      $self->_rearrange([ qw(
                              TREE_ID
                              TREE_OBJ
                              TAGNAME
                         ) ], @args);

  defined $tag                        && $self->tagname($tag);
  defined $tree_id                    && $self->tree_id($tree_id);
  defined $tree_obj                   && $self->tree($tree_obj);
  return $self;

# other possible variables to store
#                              TREE_PROGRAM
#                              TREE_METHOD
#                              TREE_FREQUENCY
#  defined $program                    && $self->program($program);
#  defined $method                     && $self->method($method);
#  defined $freq                       && $self->freq($tree_freq);

}

=head1 AnnotationI implementing functions

=cut

=head2 as_text

 Title   : as_text
 Usage   : $ann_tree->as_text();
 Function: output tree as a string
 Returns : a newic tree file
 Args    : None

=cut

sub as_text{
  my ($self) = @_;

  my $tree = $self->tree || $self->throw("Tree object absent");
  my $treeio = Bio::TreeIO->new();
  $treeio->write_tree($tree);
}

=head2 display_text

 Title   : display_text
 Usage   : my $str = $ann->display_text();
 Function: returns a string. Unlike as_text(), this method returns a string
           formatted as would be expected for te specific implementation.

           One can pass a callback as an argument which allows custom text
           generation; the callback is passed the current instance and any text
           returned
 Example :
 Returns : a string
 Args    : [optional] callback

=cut

{
  my $DEFAULT_CB = sub { $_[0]->as_text || ''};

  sub display_text {
    my ($self, $cb) = @_;
    $cb ||= $DEFAULT_CB;
    $self->throw("Callback must be a code reference") if ref $cb ne 'CODE';
    return $cb->($self);
  }

}

=head2 hash_tree

 Title   : hash_tree
 Usage   : my $hashtree = $value->hash_tree
 Function: For supporting the AnnotationI interface just returns the value
           as a hashref with the key 'value' pointing to the value
 Returns : hashrf to tree
 Args    : none

=cut

sub hash_tree{
    my $self = shift;
    my $h = {};
    $h->{'value'} = $self->tree();
    return $h;
}

=head2 tagname

 Title   : tagname
 Usage   : $obj->tagname($newval)
 Function: Get/set the tagname for this annotation value.
           Setting this is optional. If set, it obviates the need to
           provide a tag to Bio::AnnotationCollectionI when adding
           this object. When obtaining an AnnotationI object from the
           collection, the collection will set the value to the tag
           under which it was stored unless the object has a tag
           stored already.
 Returns : value of tagname (a scalar)
 Args    : new value (a scalar, optional)


=cut

sub tagname{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'tagname'} = $value;
    }
    return $self->{'tagname'};
}

=head1 Specific accessors for Tree

=head2 tree_id

 Title   : tree_id
 Usage   : $obj->tree_id($newval)
 Function: Get/set a name for the tree
 Returns : value of tagname (a scalar)
 Args    : new value (a scalar, optional)


=cut

sub tree_id {
    my $self = shift;
    return $self->{'tree_id'} = shift if defined($_[0]);
    return $self->{'tree_id'};
}

=head2 tree

 Title   : tree
 Usage   : $obj->tree($newval)
 Function: Get/set tree
 Returns : tree ref
 Args    : new value (a tree ref, optional)


=cut

sub tree {
    my $self = shift;
    return $self->{'tree'} = shift if defined($_[0]);
    return $self->{'tree'};
}

1;

