# BioPerl module for Bio::Tree::AnnotatableNode
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Mira Han <mirhan@indiana.edu>
#
# Copyright Mira Han
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::AnnotatableNode - A Tree Node with support for annotation

=head1 SYNOPSIS

    use Bio::Tree::AnnotatableNode;
    my $nodeA = Bio::Tree::AnnotatableNode->new();
    my $nodeL = Bio::Tree::AnnotatableNode->new();
    my $nodeR = Bio::Tree::AnnotatableNode->new();

    my $node = Bio::Tree::AnnotatableNode->new();
    $node->add_Descendents($nodeL);
    $node->add_Descendents($nodeR);

    print "node is not a leaf \n" if( $node->is_leaf);

    # $node is-a Bio::AnnotatableI, hence:
    my $ann_coll = $node->annotation();
    # $ann_coll is-a Bio::AnnotationCollectionI, hence:
    my @all_anns = $ann_coll->get_Annotations();
    # do something with the annotation objects

=head1 DESCRIPTION

Makes a Tree Node with Annotations, suitable for building a Tree.  See
L<Bio::Tree::Node> for a full list of functionality.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Mira Han

Email mirhan@indiana.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Tree::AnnotatableNode;
use strict;

use Bio::Annotation::Collection;
use Bio::Seq;
use base qw(Bio::Tree::Node Bio::AnnotatableI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tree::AnnotatableNode->new();
 Function: Builds a new Bio::Tree::AnnotatableNode object
 Returns : Bio::Tree::AnnotatableNode
 Args    : -tostring => code reference to the tostring callback function (optional)

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my $to_string_cb = $self->_rearrange([qw(TOSTRING)], @args);
  if ($to_string_cb) {
    $self->to_string_callback($to_string_cb);
  }
  return $self;
}

sub DESTROY {
    my ($self) = @_;
    # try to insure that everything is cleaned up
    $self->SUPER::DESTROY();
}

=head1 Methods for implementing Bio::AnnotatableI

=cut

=head2 annotation

 Title   : annotation
 Usage   : $ann = $node->annotation or 
           $node->annotation($ann)
 Function: Gets or sets the annotation
 Returns : Bio::AnnotationCollectionI object
 Args    : None or Bio::AnnotationCollectionI object
See L<Bio::AnnotationCollectionI> and L<Bio::Annotation::Collection>
for more information

=cut

sub annotation 
{
  my ($self,$value) = @_;
  if( defined $value ) {
    $self->throw("object of class ".ref($value)." does not implement ".
        "Bio::AnnotationCollectionI. Too bad.")      unless $value->isa("Bio::AnnotationCollectionI");
    $self->{'_annotation'} = $value;
  } 
  elsif( ! defined $self->{'_annotation'}) 
  {
    $self->{'_annotation'} = Bio::Annotation::Collection->new();
  }
  return $self->{'_annotation'};
}


=head1 Methods for implementing tag access through Annotation::SimpleValue

=cut

=head2 add_tag_value

 Title   : add_tag_value
 Usage   : $node->add_tag_value($tag,$value)
 Function: Adds a tag value to a node 
 Returns : number of values stored for this tag
 Args    : $tag   - tag name
           $value - value to store for the tag

=cut

sub add_tag_value 
{
  my ($self,$tag,$value) = @_;
  if( ! defined $tag || ! defined $value ) {
    $self->warn("cannot call add_tag_value with an undefined value");
  }
  my $ac = $self->annotation();
  my $sv = Bio::Annotation::SimpleValue->new(-value => $value);
  $ac->add_Annotation($tag, $sv); 
  return scalar $ac->get_Annotations($tag);
}

=head2 remove_tag

 Title   : remove_tag
 Usage   : $node->remove_tag($tag)
 Function: Remove the tag and all values for this tag
 Returns : boolean representing success (0 if tag does not exist)
 Args    : $tag - tagname to remove


=cut

sub remove_tag 
{
  my ($self,$tag) = @_;
  my $ac = $self->annotation();
  if( @{$ac->get_Annotations($tag)} ) {
    $ac->remove_Annotations($tag);
    return 1;
  }
  return 0;
}

=head2 remove_all_tags

 Title   : remove_all_tags
 Usage   : $node->remove_all_tags()
 Function: Removes all tags 
 Returns : None
 Args    : None

=cut

sub remove_all_tags
{
  my ($self) = @_;
  my $ac = $self->annotation();
  $ac->remove_Annotations();
  return;
}

=head2 get_all_tags

 Title   : get_all_tags
 Usage   : my @tags = $node->get_all_tags()
 Function: Gets all the tag names for this Node
 Returns : Array of tagnames
 Args    : None

=cut

sub get_all_tags{
  my ($self) = @_;
  my $ac = $self->annotation();
  my @tags = sort $ac->get_all_annotation_keys(); 
  # how to restrict it to SimpleValues?
  return @tags;
}

=head2 get_tag_values

 Title   : get_tag_values
 Usage   : my @values = $node->get_tag_value($tag)
 Function: Gets the values for given tag ($tag)
 Returns : Array of values or empty list if tag does not exist
 Args    : $tag - tag name

=cut

sub get_tag_values{
  my ($self,$tag) = @_;
  my $ac = $self->annotation();
  my @values = map {$_->value()} $ac->get_Annotations($tag);
  return @values;
}

=head2 has_tag

 Title   : has_tag
 Usage   : $node->has_tag($tag)
 Function: Boolean test if tag exists in the Node
 Returns : Boolean
 Args    : $tag - tagname

=cut

sub has_tag {
  my ($self,$tag) = @_;
  my $ac = $self->annotation();
  return ( scalar $ac->get_Annotations($tag) > 0);
}


=head1 Methods for implementing to_string

=cut

=head2 to_string_callback

 Title   : to_string_callback
 Usage   : $node->to_string_callback(\&func)
 Function: get/set callback for to_string
 Returns : code reference for the to_string callback function
 Args    : \&func - code reference to be set as the callback function

=cut

sub to_string_callback {
    # get/set callback, using $DEFAULT_CB if nothing is set
    my ($self, $foo) = @_;
    if ($foo) {
      # $foo is callback code ref, self as first arg (so you have access to object data)
      $self->{'_to_string_cb'} = $foo;
    }
    else {
      if (! defined $self->{'_to_string_cb'}) {
        $self->{'_to_string_cb'} = \&Bio::Tree::NodeI::to_string;
      }
    }
    return $self->{'_to_string_cb'};
}

sub to_string {
  my ($self) = @_;
  my $cb = $self->to_string_callback();
  return $cb->($self);
}

=head1 Methods for accessing Bio::Seq

=cut

=head2 sequence

 Title   : sequence
 Usage   : $ann = $node->sequence or 
           $node->sequence($seq)
 Function: Gets or sets the sequence
 Returns : array reference of Bio::SeqI objects
 Args    : None or Bio::SeqI object
See L<Bio::SeqI> and L<Bio::Seq>
for more information

=cut

sub sequence
{
  my ($self,$value) = @_;
  if( defined $value ) {
    $self->throw("object of class ".ref($value)." does not implement ".
        "Bio::SeqI. Too bad.")      unless $value->isa("Bio::SeqI");
    push (@{$self->{'_sequence'}}, $value);
  } 
  #elsif( ! defined $self->{'_sequence'}) 
  #{
  #  $self->{'_sequence'} = Bio::Seq->new();
  #}
  return $self->{'_sequence'};
}

=head2 has_sequence

 Title   : has_sequence
 Usage   : if( $node->has_sequence) { # do something } 
 Function: tells if node has sequence attached
 Returns : Boolean for whether or not node has Bio::SeqI attached.
 Args    : None 

=cut

sub has_sequence
{
  my ($self) = @_;
  return $self->{'_sequence'} && @{$self->{'_sequence'}};
}


1;
