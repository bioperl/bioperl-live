

#
# BioPerl module for Bio::AnnotationI
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

Interface all annotations must support. There are two things that each annotation
has to support.

  $annotation->as_text()

Annotations have to support an "as_text" method. This should be a
single text string, without newlines representing the annotation,
mainly for human readability. It is not aimed at being able to
store/represent the annotation

The second method allows annotations to at least attempt to represent
themselves as pure data for storage/display/whatever. The method
hash_tree

   $hash = $annotation->hash_tree();

should return an anonymous hash with "XML-like" formatting. The
formatting is as follows.

  (1) For each key in the hash, if the value is a reference'd array -

      (2) For each element of the array if the value is a object - 
          Assumme the object has the method "hash_tree";
      (3) else if the value is a referene to a hash
          Recurse again from point (1)
      (4) else 
          Assumme the value is a scalar, and handle it directly as text

   (5) else (if not an array) apply rules 2,3 and 4 to value

The XML path in tags is represented by the keys taken in the
hashes. When arrays are encountered they are all present in the path
level of this tag

This is a pretty "natural" representation of an object tree in an XML
style, without forcing everything to inheriet off some super-generic
interface for representing things in the hash.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bio.perl.org

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...


package Bio::AnnotationI;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;


@ISA = qw(Bio::Root::Root);
# new() can be inherited from Bio::Root::Root

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

   $self->throw("Implementing object did not provide as_text function");
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
   my ($self,@args) = @_;

   $self->throw("Implementing object did not provide hash_tree method");
}

1;
