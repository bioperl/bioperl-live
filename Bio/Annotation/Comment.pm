
#
# BioPerl module for Bio::Pfam::Annotation::Comment
#
# Cared for by Ewan Birney <pfam@sanger.ac.uk>
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
    $annotation->add_Comment($comment);


=head1 DESCRIPTION

A holder for comments in annotations, just plain text. This is a very simple
object, and justifably so.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Annotation::Comment;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::AnnotationI;

@ISA = qw(Bio::AnnotationI Bio::Root::Root);

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
  my ($text) = $self->_rearrange([qw( TEXT )], @args);

  defined $text && $self->text($text);

  return $self;
}

=head2 AnnotationI implementing functions

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

=head2 hash_tree

 Title   : hash_tree
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub hash_tree{
   my ($self) = @_;
   
   my $h = {};
   $h->{'text'} = $self->text;
}

=head2 Specific accessors for Comments

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



1;
