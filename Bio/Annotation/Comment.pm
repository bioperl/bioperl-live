#
# BioPerl module for Bio::Annotation::Comment
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Annotation::Comment - A comment object, holding text

=head1 SYNOPSIS


    $comment = Bio::Annotation::Comment->new();
    $comment->text("This is the text of this comment");
    $annotation->add_Annotation('comment', $comment);


=head1 DESCRIPTION

A holder for comments in annotations, just plain text. This is a very simple
object, and justifiably so.

=head1 AUTHOR - Ewan Birney 

Email birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Annotation::Comment;
use strict;

use base qw(Bio::Root::Root Bio::AnnotationI);

=head2 new

 Title   : new
 Usage   : $comment = Bio::Annotation::Comment->new( '-text' => 'some text for this comment');
 Function: This returns a new comment object, optionally with
           text filed
 Example :
 Returns : a Bio::Annotation::Comment object
 Args    : a hash with -text optionally set


=cut


sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($text,$tag, $type) = $self->_rearrange([qw(TEXT TAGNAME TYPE)], @args);

  defined $text && $self->text($text);
  defined $tag && $self->tagname($tag);
  defined $type && $self->type($type);
  return $self;
}

=head1 AnnotationI implementing functions

=cut

=head2 as_text

 Title   : as_text
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub as_text{
   my ($self) = @_;

   return "Comment: ".$self->text;
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
  my $DEFAULT_CB = sub {$_[0]->text || ''};

  sub display_text {
    my ($self, $cb) = @_;
    $cb ||= $DEFAULT_CB;
    $self->throw("Callback must be a code reference") if ref $cb ne 'CODE';
    return $cb->($self);
  }

}

=head2 hash_tree

 Title   : hash_tree
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub hash_tree{
    my $self = shift;
   
    my $h = {};
    $h->{'text'} = $self->text;
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

 Example : 
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

=head1 Specific accessors for Comments

=cut


=head2 text

 Title   : text
 Usage   : $value = $self->text($newval)
 Function: get/set for the text field. A comment object
           just holds a single string which is accessible through
           this method
 Example : 
 Returns : value of text
 Args    : newvalue (optional)


=cut

sub text{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'text'} = $value;
    }
    return $self->{'text'};

}

=head2 value

 Title   : value
 Usage   : $value = $self->value($newval)
 Function: Alias of the 'text' method
 Example :
 Returns : value of text
 Args    : newvalue (optional)


=cut


*value = \&text;

=head2 type

 Title   : type
 Usage   : $value = $self->type($newval)
 Function: get/set for the comment type field.  The comment type
           is normally found as a subfield within comment sections
           in some files, such as SwissProt
 Example : 
 Returns : value of text
 Args    : newvalue (optional)


=cut

sub type {
   my ($self,$type) = @_;
   if( defined $type) {
      $self->{'type'} = $type;
    }
    return $self->{'type'};

}

1;
