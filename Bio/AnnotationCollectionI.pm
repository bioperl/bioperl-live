# $Id$

#
# BioPerl module for Bio::AnnotationCollectionI
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::AnnotationCollectionI - Interface for annotation collections

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

Annotation Collections are a way of storing a series of "interesting
facts" about something. We call an "interesting fact" in Bioperl an
Annotation (this differs from a Sequence Feature, which is called
a Sequence Feature and may or may not have an Annotation Collection).

The trouble about this is we are not that sure what "interesting
facts" someone might want to store: the possibility is endless. 

Bioperl''s approach is that the "interesting facts" are represented by
Bio::AnnotationI objects. The interface Bio::AnnotationI guarantees
two methods

   $obj->as_text(); # string formated to display to users

and

   $obj->hash_tree(); # hash with defined rules for data-orientated discovery

The hash_tree method is designed to play well with XML output and
other "nested-tag-of-data-values" think BoulderIO and/or Ace stuff. For more
info read Bio::AnnotationI docs

Annotations are stored in AnnotationCollections, each Annotation under a
different "tag". The tags allow simple discovery of the available annotations,
and in some cases (like the tag "gene_name") indicate how to interpret the
data underneath the tag. The tag is only one tag deep and each tag can have an
array of values.

In addition, AnnotationCollections are guaranteed to maintain consistent
types of objects under each tag - at least that each object complies to one
interface. The "standard" AnnotationCollection insists the following rules
are set up

  Tag            Object
  ---            ------
  comment        Bio::Annotation::Comment
  dblink         Bio::Annotation::DBLink
  description    Bio::Annotation::SimpleValue
  gene_name      Bio::Annotation::SimpleValue
  ontology_term  Bio::Annotation::OntologyTerm
  reference      Bio::Annotation::Reference

These tags are the implict tags that the SeqIO system needs to round-trip
GenBank/EMBL/Swissprot.

However, you as a user and us collectively as a community can grow the
"standard" tag mapping over time and specifically for a particular
area.


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bio.perl.org

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::AnnotationCollectionI;
use vars qw(@ISA);
use strict;

# Interface preamble - inherits from Bio::Root::RootI

use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);

=head1 ACCESSOR METHODS

Use these for Bio::AnnotationI object access.

=cut

=head2 get_all_annotation_keys()

 Usage   : $ac->get_all_annotation_keys()
 Function: gives back a list of annotation keys, which are simple text strings
 Returns : list of strings
 Args    : none

=cut

sub get_all_annotation_keys{
    shift->throw_not_implemented();
}


=head2 get_Annotations()

 Usage   : my @annotations = $collection->get_Annotations('key')
 Function: Retrieves all the Bio::AnnotationI objects for a specific key
 Returns : list of Bio::AnnotationI - empty if no objects stored for a key
 Args    : string which is key for annotations

=cut

sub get_Annotations{
    shift->throw_not_implemented();
}

=head2 add_Annotation()

 Usage   : $self->add_Annotation('reference',$object);
           $self->add_Annotation($object,'Bio::MyInterface::DiseaseI');
           $self->add_Annotation($object);
           $self->add_Annotation('disease',$object,'Bio::MyInterface::DiseaseI');
 Function: Adds an annotation for a specific key.

           If the key is omitted, the object to be added must provide a value
           via its tagname().

           If the archetype is provided, this and future objects added under
           that tag have to comply with the archetype and will be rejected
           otherwise.

 Returns : none
 Args    : annotation key ('disease', 'dblink', ...)
           object to store (must be Bio::AnnotationI compliant)
           [optional] object archetype to map future storage of object
           of these types to

=cut

sub add_Annotation {
  shift->throw_not_implemented();
}

=head2 remove_Annotations()

 Usage   :
 Function: Remove the annotations for the specified key from this collection.
 Returns : an list of Bio::AnnotationI compliant objects which were stored
           under the given key(s)
 Args    : the key(s) (tag name(s), one or more strings) for which to
           remove annotations (optional; if none given, flushes all
           annotations)

=cut

sub remove_Annotations{
  shift->throw_not_implemented();
}

=head2 get_num_of_annotations()

 Usage   : my $count = $collection->get_num_of_annotations()
 Function: Returns the count of all annotations stored in this collection 
 Returns : integer
 Args    : none

=cut

sub get_num_of_annotations{
    shift->throw_not_implemented();
}

=head1 "*_tag_*" METHODS

The methods below allow mapping of the old "get_tag_values()"-style
annotation access to Bio::AnnotationCollectionI.  These need not be
implemented in a Bio::AnnotationCollectionI compliant class, as they
are built on top of the methods (see above L</ACCESSOR METHODS>).

 B<DEPRECATED>: DO NOT USE THESE FOR FUTURE DEVELOPMENT.

=cut

=head2 has_tag()

 Usage   : $count = $obj->has_tag($tag)
 Function: returns the number of annotations corresponding to $tag
 Returns : an integer
 Args    : tag name
 Note    : B<DEPRECATED>: this method is essentially scalar(L</get_Annotations()>).

=cut

sub has_tag {
  my ($self,$tag) = @_;
  $self->warn('has_tag() is deprecated.  use get_Annotations()');
  return scalar($self->get_Annotations($tag));
}

=head2 add_tag_value()

 Usage   : See L</add_Annotation()>.
 Function:
 Returns : 
 Args    : B<DEPRECATED>: this method is essentially L</add_Annotation()>.

=cut

sub add_tag_value {
  my ($self,@args) = @_;
  $self->warn('add_tag_value() is deprecated.  use add_Annotation()');
  return $self->add_Annotation(@args);
}

=head2 get_tag_values()

 Usage   : @annotations = $obj->get_tag_values($tag)
 Function: returns annotations corresponding to $tag
 Returns : a list of Bio::AnnotationI objects
 Args    : tag name
 Note    : B<DEPRECATED>: this method is essentially L</get_Annotations()>.

=cut

sub get_tag_values {
  my ($self,$tag) = @_;
  $self->warn('get_tag_values() is deprecated.  use get_Annotations()');
  return $self->get_Annotations($tag);
}

=head2 get_tagset_values()

 Usage   : @annotations = $obj->get_tagset_values($tag1,$tag2)
 Function: returns annotations corresponding to a list of tags.
           this is a convenience method equivalent to multiple calls
           to L</get_tag_values> with each tag in the list.
 Returns : a list of Bio::AnnotationI objects.
 Args    : a list of tag names
 Note    : B<DEPRECATED>: this method is essentially multiple calls to
           L</get_Annotations()>.

=cut

sub get_tagset_values {
  my ($self,@tags) = @_;
  $self->warn('get_tagset_values() is deprecated.  use get_Annotations()');
  my @r = ();
  foreach my $tag (@tags){
    push @r, $self->get_Annotations($tag);
  }
  return @r;
}

=head2 get_all_tags()

 Usage   : @tags = $obj->get_all_tags()
 Function: returns a list of annotation tag names.
 Returns : a list of tag names
 Args    : none
 Note    : B<DEPRECATED>: use L</get_all_annotation_keys()>.

=cut

sub get_all_tags {
  my ($self,@args) = @_;
  $self->warn('get_all_tags() is deprecated.  use get_all_annotation_keys()');
  return $self->get_all_annotation_keys(@args);
}

=head2 remove_tag()

 Usage   : See L</remove_Annotations()>.
 Function:
 Returns : 
 Args    : B<DEPRECATED>: use L</remove_Annotations()>.
 Note    : Contrary to what the name suggests, this method removes
           B<all> annotations corresponding to $tag, not just a
           single anntoation.

=cut

sub remove_tag {
  my ($self,@args) = @_;
  $self->warn('remove_tag() is deprecated.  use remove_Annotations()');
  return $self->remove_Annotations(@args);
}

1;
