

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

Bioperl\'s approach is that the "interesting facts" are represented by
Bio::AnnotationI objects. The interface Bio::AnnotationI guarentees
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

In addition, AnnotationCollectionI\'s are guarentee to maintain a consistent
set object values under each tag - at least that each object complies to one
interface. The "standard" AnnotationCollection insists the following rules
are set up

  Tag         Object
  ---         ------
  reference   Bio::Annotation::Reference
  comment     Bio::Annotation::Comment
  dblink      Bio::Annotation::DBLink
  gene_name   Bio::Annotation::SimpleValue
  description Bio::Annotation::SimpleValue

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
  http://bio.perl.org/bioperl-bugs/

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



=head2 get_all_annotation_keys

 Title   : get_all_annotation_keys
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_annotation_keys{
    my ($self) = @_;
    $self->throw("Abstract method get_all_annotation_keys implementing class did not provide method");
}


=head2 get_Annotations

 Title   : get_Annotations
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Annotations{
    my ($self) = @_;

    $self->throw("Abstract method get_Annotations implementing class did not provide method");
    
}

1;
