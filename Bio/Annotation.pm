
#
# BioPerl module for Bio::Annotation
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Kevin Howe, Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Annotation - A generic object for annotations

=head1 SYNOPSIS

    # get an annotation somehow in $ann

    # description is a simple, one line description 
    print "Description is ",$ann->description "\n";


    foreach $comment ( $ann->each_Comment ) {
       # $comment is a Bio::Annotation::Comment object
          print "Comment: ", $comment->text(), "\n"
       }
    }

    foreach $link ( $ann->each_DBLink ) {
       # link is a Bio::Annotation::DBLink object
       print "Link to ",$link->primary_id, " in database", $link->db "\n";
    }

    foreach $ref ( $ann->each_Reference ) {
       # link is a Bio::Annotation::Reference object
       print "Reference title ", $ref->title , "\n";
    }

    #
    # Making an annotation object from scratch
    #

    $ann = Bio::Pfam::Annotation->new();

    $ann->description("Description text");
    print "Annotation description is ", $ann->description, "\n";
   

=head1 DESCRIPTION

The object represents generic biological annotation of an object. It
has the ability to provide 

    a brief, one line description
    free text comments
    links to other biological objects
    references to literature

It does not have the following abilities

    The basis (experimental/non experimental/homology) 
       of the annotation. This is considered to be part of
       the object which owns the annotation. This is 
       because the type of relevant basis is usually 
       dependent on the object

    The previous revisions of the object
       This should be a property of whatever database this
       object comes from


=head1 CONTACT

Mail birney@sanger.ac.uk with any queries

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Annotation;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;

# we don't really need these object but we should do.
# declare them here to prevent tears later.

use Bio::Annotation::Reference;
use Bio::Annotation::DBLink;
use Bio::Annotation::Comment;

@ISA = qw(Bio::Root::Object);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self, %params) = @_;
  my( $text ) = ( $params{'-DESCRIPTION'}||$params{'-description'} );

  my $make = $self->SUPER::_initialize;

  $self->{ 'description' } = $text;
  $self->{ 'refs' } = [];
  $self->{ 'comment' } = [];
  $self->{ 'link' } = [];
  return $make; # success - we hope!
}


=head2 description

 Title   : description
 Usage   : $obj->description($newval)
 Function: 
 Example : 
 Returns : value of description
 Args    : newvalue (optional)


=cut

sub description{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'description'} = $value;
    }
    return $self->{'description'};

}

=head2 gene_name

 Title   : gene_name
 Usage   : $obj->gene_name($newval)
 Function: 
 Example : 
 Returns : value of gene name
 Args    : newvalue (optional)


=cut

sub gene_name{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'gene_name'} = $value;
    }
    return $self->{'gene_name'};

}


=head2 add_Reference

 Title   : add_Reference
 Usage   : $self->add_Reference($ref)
 Function: adds a reference object
 Example :
 Returns : 
 Args    :


=cut

sub add_Reference{
   my ($self) = shift;
   foreach my $ref ( @_ ) {
       push(@{$self->{'refs'}},$ref);
   }
}

=head2 each_Reference

 Title   : each_Reference
 Usage   : foreach $ref ( $self->each_Reference() )
 Function: gets an array of reference
 Example :
 Returns : 
 Args    :


=cut

sub each_Reference{
   my ($self,@args) = @_;
   
   return @{$self->{'refs'}}; 
}



=head2 add_Comment

 Title   : add_Comment
 Usage   : $self->add_Comment($ref)
 Function: adds a Comment object
 Example :
 Returns : 
 Args    :


=cut

sub add_Comment{
   my ($self) = shift;
   foreach my $com ( @_ ) {
       if( ! $com->isa('Bio::Annotation::Comment') ) {
	   $self->throw("Is not a comment object but a  [$com]");
       }

       push(@{$self->{'comment'}},$com);
   }
}

=head2 each_Comment

 Title   : each_Comment
 Usage   : foreach $ref ( $self->each_Comment() )
 Function: gets an array of Comment of objects
 Example :
 Returns : 
 Args    :


=cut

sub each_Comment{
   my ($self) = @_;
   
   return @{$self->{'comment'}}; 
}


=head2 add_DBLink

 Title   : add_DBLink
 Usage   : $self->add_DBLink($ref)
 Function: adds a link object
 Example :
 Returns : 
 Args    :


=cut

sub add_DBLink{
   my ($self,$com) = @_;
   push(@{$self->{'link'}},$com);
}

=head2 each_DBLink

 Title   : each_DBLink
 Usage   : foreach $ref ( $self->each_DBlink() )
 Function: gets an array of DBlink of objects
 Example :
 Returns : 
 Args    :


=cut

sub each_DBLink{
   my ($self) = @_;
   
   return @{$self->{'link'}}; 
}


1;





