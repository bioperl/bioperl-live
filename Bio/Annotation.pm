#
# $Id$
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
    print "Description is ",$ann->description, "\n";

    foreach $genename ( $self->each_gene_name() ) {
	print "gene name: $genename\n";
    }

    foreach $comment ( $ann->each_Comment ) {
       # $comment is a Bio::Annotation::Comment object
       print "Comment: ", $comment->text(), "\n";
    }

    foreach $link ( $ann->each_DBLink ) {
       # link is a Bio::Annotation::DBLink object
       print "Link to ",$link->primary_id, " in database", $link->database, "\n";
    }

    foreach $ref ( $ann->each_Reference ) {
       # link is a Bio::Annotation::Reference object
       print "Reference title ", $ref->title , "\n";
    }

    #
    # Making an annotation object from scratch
    #

    $ann = Bio::Annotation->new();

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

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/


=head1 AUTHOR - Ewan Birney 

Email birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Annotation;
use vars qw(@ISA);
use strict;

use Bio::Root::RootI;

# we don't really need these object but we should 
# declare them here to prevent tears later.

use Bio::Annotation::Reference;
use Bio::Annotation::DBLink;
use Bio::Annotation::Comment;

@ISA = qw(Bio::Root::RootI);

=head2 new

 Title   : new
 Usage   : $annotation = Bio::Annotation->new( '-description' => 'a description line');
 Function: Makes a new Annotation object. The main thing 
           you will want to do with this is add comment objects and
           dblink objects, with calls like

            $annotation->add_Comment($comment);
            $annotation->add_DBLink($dblink);

 Example :
 Returns : a new Bio::Annotation Object
 Args    : hash, potentially with one field, -description


=cut


sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($text ) = $self->_rearrange([qw(DESCRIPTION)], @args);
  
  defined $text && $self->description($text);
  $self->{ 'refs' } = [];
  $self->{ 'comment' } = [];
  $self->{ 'link' } = [];
  $self->{ '_names' } = [];

  return $self; 
}


=head2 description

 Title   : description
 Usage   : $obj->description($newval)
 Function: 
 Example : 
 Returns : value of description
 Args    : newvalue (optional)


=cut

sub description {
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'description'} = $value;
    }
    return $self->{'description'};
}

#  =head2 gene_name

#   Title   : gene_name
#   Usage   : $obj->gene_name($newval)
#   Function: Get/set the primary gene name.

#             Use of this method is deprecated. Use add_gene_name()/
#             each_gene_name() instead.
#   Example : 
#   Returns : value of gene name
#   Args    : newvalue (optional)


#  =cut

sub gene_name{
    my ($self,$value) = @_;

    $self->warn("gene_name() is deprecated. Use add_gene_name/each_gene_name instead.");
    if( defined $value) {
	$self->add_gene_name($value);
    }
    ($value) = $self->each_gene_name();
    return $value;
}


=head2 add_gene_name

 Title   : add_gene_name
 Usage   : $self->add_gene_name($name1[,$name2,...])
 Function: adds a reference object
 Example :
 Returns : 
 Args    : a string, or a list of strings


=cut

sub add_gene_name{
    my ($self) = shift;
    foreach my $name ( @_ ) {
	push(@{$self->{'_names'}},$name);
    }
}

=head2 remove_gene_name

 Title   : remove_gene_name
 Usage   : $self->remove_gene_name($index)
 Function: removes a particular gene name
 Example :
 Returns : 
 Args    : index of the name to remove


=cut

sub remove_gene_name{
   my ($self,$idx) = @_;
   splice @{$self->{'_names'}}, $idx, 1;
}


=head2 each_gene_name

 Title   : each_gene_name
 Usage   : foreach $genename ( $self->each_gene_name() ) {
               print "seq has gene name $genename\n";
           }
 Function: gets the array of gene names
 Example :
 Returns : an array of strings
 Args    :


=cut

sub each_gene_name{
   my ($self,@args) = @_;
   
   return @{$self->{'_names'}}; 
}

=head2 add_Reference

 Title   : add_Reference
 Usage   : $self->add_Reference($ref1[,$ref2,...])
 Function: adds a reference object
 Example :
 Returns : 
 Args    : a Bio::Annotation::Reference or derived object


=cut

sub add_Reference{
   my ($self) = shift;
   foreach my $ref ( @_ ) {
       if( ! $ref->isa('Bio::Annotation::Reference') ) {
	   $self->throw("Is not a Bio::Annotation::Reference object but a [$ref]");
       }
       push(@{$self->{'refs'}},$ref);
   }
}

=head2 remove_Reference

 Title   : remove_Reference
 Usage   : $self->remove_Reference($index)
 Function: removes a reference object
 Example :
 Returns : 
 Args    : index number from references array


=cut

sub remove_Reference{
   my ($self,$idx) = @_;
   splice @{$self->{'refs'}}, $idx, 1;
}


=head2 each_Reference

 Title   : each_Reference
 Usage   : foreach $ref ( $self->each_Reference() )
 Function: gets an array of reference
 Example :
 Returns : an array of Bio::Annotation::Reference or derived objects
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
 Args    : a Bio::Annotation::Comment or derived object


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

=head2 remove_Comment

 Title   : remove_Comment
 Usage   : $self->remove_Comment($index)
 Function: removes a comment object
 Example :
 Returns : 
 Args    : index number from comments array


=cut

sub remove_Comment{
   my ($self,$idx) = @_;
   splice @{$self->{'comment'}}, $idx, 1;
}

=head2 each_Comment

 Title   : each_Comment
 Usage   : foreach $ref ( $self->each_Comment() )
 Function: gets an array of Comment of objects
 Example :
 Returns : an array of Bio::Annotation::Comment or derived objects
 Args    : none


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
 Args    : a Bio::Annotation::DBLink or derived object


=cut

sub add_DBLink{
   my ($self,$com) = @_;
   if( ! $com->isa('Bio::Annotation::DBLink') ) {
       $self->throw("Is not a link object but a  [$com]");
   }
   push(@{$self->{'link'}},$com);
}

=head2 remove_DBLink

 Title   : remove_DBLink
 Usage   : $self->remove_DBLink($index)
 Function: removes a DBLink object
 Example :
 Returns : 
 Args    : index number from links array


=cut

sub remove_DBLink{
   my ($self,$idx) = @_;
   splice @{$self->{'link'}}, $idx, 1;
}


=head2 each_DBLink

 Title   : each_DBLink
 Usage   : foreach $ref ( $self->each_DBlink() )
 Function: gets an array of DBlink of objects
 Example :
 Returns : an array of Bio::Annotation::DBlink or derived objects
 Args    :


=cut

sub each_DBLink{
   my ($self) = @_;
   
   return @{$self->{'link'}}; 
}


1;







