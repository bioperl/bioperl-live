#
# BioPerl module for Bio::AnnotationI
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

Bio::AnnotationI - Annotation interface

=head1 SYNOPSIS

  # generally you get AnnotationI's from AnnotationCollectionI's

   foreach $key ( $ac->get_all_annotation_keys() ) {
       @values = $ac->get_Annotations($key);
       foreach $value ( @values ) {
          # value is an Bio::AnnotationI, and defines a "as_text" method
          print "Annotation ",$key," stringified value ",$value->as_text,"\n";
          # you can also use a generic hash_tree method for getting
          # stuff out say into XML format
          $hash_tree = $value->hash_tree();
       }
   }

=head1 DESCRIPTION

Interface all annotations must support. There are two things that each
annotation has to support.

  $annotation->as_text()

Annotations have to support an "as_text" method. This should be a
single text string, without newlines representing the annotation,
mainly for human readability. It is not aimed at being able to
store/represent the annotation.

The second method allows annotations to at least attempt to represent
themselves as pure data for storage/display/whatever. The method
hash_tree should return an anonymous hash with "XML-like" formatting:

   $hash = $annotation->hash_tree();

The formatting is as follows.

  (1) For each key in the hash, if the value is a reference'd array -

  (2) For each element of the array if the value is a object -
          Assume the object has the method "hash_tree";
  (3) else if the value is a reference to a hash
          Recurse again from point (1)
  (4) else
          Assume the value is a scalar, and handle it directly as text
  (5) else (if not an array) apply rules 2,3 and 4 to value

The XML path in tags is represented by the keys taken in the
hashes. When arrays are encountered they are all present in the path
level of this tag

This is a pretty "natural" representation of an object tree in an XML
style, without forcing everything to inherit off some super-generic
interface for representing things in the hash.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...


package Bio::AnnotationI;
use strict;

# Object preamble - inherits from Bio::Root::Root

use base qw(Bio::Root::RootI);

=head2 as_text

 Title   : as_text
 Usage   :
 Function: single text string, without newlines representing the
           annotation, mainly for human readability. It is not aimed
           at being able to store/represent the annotation.
 Example :
 Returns : a string
 Args    : none


=cut

sub as_text{
    shift->throw_not_implemented();
}

=head2 display_text

 Title   : display_text
 Usage   : my $str = $ann->display_text();
 Function: returns a string. Unlike as_text(), this method returns a string
           formatted as would be expected for the specific implementation.

           Implementations should allow passing a callback as an argument which
           allows custom text generation; the callback will be passed the
           current implementation.

           Note that this is meant to be used as a simple representation
           of the annotation data but probably shouldn't be used in cases
           where more complex comparisons are needed or where data is
           stored.
 Example :
 Returns : a string
 Args    : [optional] callback

=cut

sub display_text {
    shift->throw_not_implemented();
}

=head2 hash_tree

 Title   : hash_tree
 Usage   :
 Function: should return an anonymous hash with "XML-like" formatting
 Example :
 Returns : a hash reference
 Args    : none

=cut

sub hash_tree{
    shift->throw_not_implemented();
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
    shift->throw_not_implemented();
}

1;
