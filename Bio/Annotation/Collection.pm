# $Id$

#
# BioPerl module for Bio::Annotation::Collection.pm
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Annotation::Collection - Default Perl implementation of AnnotationCollectionI

=head1 SYNOPSIS

   # get an AnnotationCollectionI somehow, eg

   $ac = $seq->annotation();

   foreach $key ( $ac->get_all_annotation_keys() ) {
       @values = $ac->get_Annotations($key);
       foreach $value ( @values ) {
          # value is an Bio::AnnotationI, and defines a "as_text" method
          print "Annotation ",$key," stringified value ",$value->as_text,"\n";

          # also defined hash_tree method, which allows data orientated
          # access into this object
          $hash = $value->hash_tree();
       }
   }

=head1 DESCRIPTION

Bioperl implementation for Bio::AnnotationCollecitonI 

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bio.perl.org/MailList.html  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bioperl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Annotation::Collection;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::AnnotationCollectionI;
use Bio::Root::Root;
use Bio::Annotation::TypeManager;
use Bio::Annotation::SimpleValue;


@ISA = qw(Bio::Root::Root Bio::AnnotationCollectionI);


=head2 new

 Title   : new
 Usage   : $coll = Bio::Annotation::Collection->new()
 Function: Makes a new Annotation::Collection object. 
 Returns : Bio::Annotation::Collection
 Args    : none

=cut

sub new{
   my ($class,@args) = @_;

   my $self = $class->SUPER::new(@args);

   $self->{'_annotation'} = {};
   $self->_typemap(Bio::Annotation::TypeManager->new());

   return $self;
}


=head2 Bio::Annotation::CollectionI implementing methods

=cut

=head2 get_all_annotation_keys

 Title   : get_all_annotation_keys
 Usage   : $ac->get_all_annotation_keys()
 Function: gives back a list of annotation keys, which are simple text strings
 Returns : list of strings
 Args    : none

=cut

sub get_all_annotation_keys{
   my ($self) = @_;
   return keys %{$self->{'_annotation'}};
}

=head2 get_Annotations

 Title   : get_Annotations
 Usage   : my @annotations = $collection->get_Annotations('key')
 Function: Retrieves all the Bio::AnnotationI objects for a specific key
 Returns : list of Bio::AnnotationI - empty if no objects stored for a key
 Args    : string which is key for annotations

=cut

sub get_Annotations{
   my ($self,$key) = @_;

   if( !defined $self->{'_annotation'}->{$key} ) {
       return ();
   } else {
       return @{$self->{'_annotation'}->{$key}};
   }
}

=head2 get_num_of_annotations

 Title   : get_num_of_annotations
 Usage   : my $count = $collection->get_num_of_annotations()
 Function: Returns the count of all annotations stored in this collection 
 Returns : integer
 Args    : none


=cut

sub get_num_of_annotations{
   my ($self) = @_;
   my $count = 0;
   map { $count += scalar @$_ } values %{$self->{'_annotation'}};
   return $count;
}

=head2 Implementation specific functions - mainly for adding

=cut

=head2 add_Annotation

 Title   : add_Annotation
 Usage   : $self->add_Annotation('reference',$object);
           $self->add_Annotation('disease',$object,'Bio::MyInterface::DiseaseI');
 Function: Adds an annotation for a specific 
 Returns : none
 Args    : annotation key ('disease', 'dblink', ...)
           object to store (must be Bio::AnnotationI compliant)
           [optional] object archytype to map future storage of object 
                      of these types to

=cut

sub add_Annotation{
   my ($self,$key,$object,$archytype) = @_;
   
   if( !defined $object ) {
       $self->throw("Must have at least key and object in add_Annotation");
   }

   if( !ref $object ) {
       $self->throw("Must add an object. Use Bio::Annotation::Comment, SimpleValue or ControledVocabTerm for simple text additions");
   }

   if( !$object->isa("Bio::AnnotationI") ) {
       $self->throw("object must be AnnotationI compliant, otherwise we wont add it!");
   }

   # ok, now we are ready! If we don't have an archytype, set it
   # from the type of the object

   if( !defined $archytype ) {
       $archytype = ref $object;
   }

   # check typemap, storing if needed.
   my $stored_map = $self->_typemap->type_for_key($key);

   if( defined $stored_map ) {
       # check validity, irregardless of archytype. A little cheeky
       # this means isa stuff is executed correctly

       if( !$self->_typemap()->is_valid($key,$object) ) {
	   $self->throw("Object $object was not valid with key $key. If you were adding new keys in, perhaps you want to make use of the archytype method to allow registration to a more basic type");
       }
   } else {
       $self->_typemap->_add_type_map($key,$archytype);
   }

   # we are ok to store

   if( !defined $self->{'_annotation'}->{$key} ) {
       $self->{'_annotation'}->{$key} = [];
   }

   push(@{$self->{'_annotation'}->{$key}},$object);

   return 1;
}


=head2 Backward compatible functions

Functions put in for backward compatibility with old
Bio::Annotation.pm stuff

=cut

=head2 description

 Title   : description
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub description{
   my ($self,$value) = @_;

   $self->deprecated("Using old style annotation call on new Annotation::Collection object");

   if( defined $value ) {
       my $val = Bio::Annotation::SimpleValue->new();
       $val->value($value);
       $self->add_Annotation('description',$val);
   }

   my ($desc) = $self->get_Annotations('description');
   
   # If no description tag exists, do not attempt to call value on undef:
   return $desc ? $desc->value : undef;
}


=head2 add_gene_name

 Title   : add_gene_name
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_gene_name{
   my ($self,$value) = @_;

   $self->deprecated("Old style add_gene_name called on new style Annotation::Collection");

   my $val = Bio::Annotation::SimpleValue->new();
   $val->value($value);
   $self->add_Annotation('gene_name',$val);
}

=head2 each_gene_name

 Title   : each_gene_name
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_gene_name{
   my ($self) = @_;

   $self->deprecated("Old style each_gene_name called on new style Annotation::Collection");

   my @out;
   my @gene = $self->get_Annotations('gene_name');

   foreach my $g ( @gene ) {
       push(@out,$g->value);
   }

   return @out;
}

=head2 add_Reference

 Title   : add_Reference
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Reference{
   my ($self, @values) = @_;

   $self->deprecated("add_Reference (old style Annotation) on new style Annotation::Collection");
   
   # Allow multiple (or no) references to be passed, as per old method
   foreach my $value (@values) {
       $self->add_Annotation('reference',$value);
   }
}

=head2 each_Reference

 Title   : each_Reference
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_Reference{
   my ($self) = @_;

   $self->deprecated("each_Reference (old style Annotation) on new style Annotation::Collection");
   
   return $self->get_Annotations('reference');
}


=head2 add_Comment

 Title   : add_Comment
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Comment{
   my ($self,$value) = @_;

   $self->deprecated("add_Comment (old style Annotation) on new style Annotation::Collection");

   $self->add_Annotation('comment',$value);

}

=head2 each_Comment

 Title   : each_Comment
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_Comment{
   my ($self) = @_;

   $self->deprecated("each_Comment (old style Annotation) on new style Annotation::Collection");
   
   return $self->get_Annotations('comment');
}



=head2 add_DBLink

 Title   : add_DBLink
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_DBLink{
   my ($self,$value) = @_;

   $self->deprecated("add_DBLink (old style Annotation) on new style Annotation::Collection");

   $self->add_Annotation('dblink',$value);

}

=head2 each_DBLink

 Title   : each_DBLink
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_DBLink{
   my ($self) = @_;

   $self->deprecated("each_DBLink (old style Annotation) on new style Annotation::Collection");
   
   return $self->get_Annotations('dblink');
}



=head2 Implementation management functions

=cut

=head2 _typemap

 Title   : _typemap
 Usage   : $obj->_typemap($newval)
 Function: 
 Example : 
 Returns : value of _typemap
 Args    : newvalue (optional)


=cut

sub _typemap{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_typemap'} = $value;
    }
    return $self->{'_typemap'};

}

1;
