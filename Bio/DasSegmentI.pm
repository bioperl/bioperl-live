# $Id$
#
# BioPerl module for Bio::DasI
#
# Cared for by Lincoln Stein <lstein@cshl.org>
#
# Copyright Lincoln Stein
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DasSegmentI - DAS-style access to a feature database

=head1 SYNOPSIS

  # Get a Bio::DasSegmentI object from a Bio::DasI database...


  @features = $segment->features(-type=>['type1','type2']);
  # each feature is a Bio::SeqFeatureI-compliant object

  $stream = $segment->get_feature_stream(-type=>['type1','type2','type3'];
  while (my $feature = $stream->next_seq) {
     # do something with feature
  }

  $count = $segment->features_callback(-type=>['type1','type2','type3'],
                                       -callback => sub { ... { }
                                       );

=head1 DESCRIPTION

Bio::DasSegmentI is a simplified alternative interface to sequence annotation
databases used by the distributed annotation system. In this scheme,
the genome is represented as a series of landmarks.  You fetch regions
of the genome by specifying a landmark, and optionally a start and end
position relative to the landmark.

Interface all
annotations must support. There are two things that each annotation
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
