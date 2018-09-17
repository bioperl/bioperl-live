#
# BioPerl module for PathI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Hilmar Lapp <hlapp at gmx.net>
#
# (c) Hilmar Lapp, hlapp at gmx.net, 2003.
# (c) GNF, Genomics Institute of the Novartis Research Foundation, 2003.
#
# You may distribute this module under the same terms as perl itself.
# Refer to the Perl Artistic License (see the license accompanying this
# software package, or see http://www.perl.com/language/misc/Artistic.html)
# for the terms under which you may use, modify, and redistribute this module.
#
# THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
# MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Ontology::PathI - Interface for a path between ontology terms

=head1 SYNOPSIS

    # see documentation of methods and an implementation, e.g.,
    # Bio::Ontology::Path

=head1 DESCRIPTION

This is the minimal interface for a path between two terms in
an ontology. Ontology engines may use this.

Essentially this is a very thin extension of the
L<Bio::Ontology::RelationshipI> interface. It basically adds an
attribute distance(). For a RelationshipI, you can think of distance as
equal to zero (subject == object) or 1 (subject != object).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Hilmar Lapp

Email hlapp at gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Ontology::PathI;
use strict;

use base qw(Bio::Ontology::RelationshipI);


=head2 distance

 Title   : distance
 Usage   : $obj->distance($newval)
 Function: Get (and set if the implementation allows it) the distance
           between the two terms connected by this path.

 Example : 
 Returns : value of distance (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub distance{
    return shift->throw_not_implemented();
}

=head1 Bio::Ontology::RelationshipI Methods

=cut

=head2 subject_term

 Title   : subject_term
 Usage   : $subj = $rel->subject_term();
 Function: Set/get for the subject term of this Relationship.

           The common convention for ontologies is to express
           relationships between terms as triples (subject, predicate,
           object).

 Returns : The subject term [Bio::Ontology::TermI].
 Args    : 

=cut

=head2 object_term

 Title   : object_term
 Usage   : $object = $rel->object_term();
 Function: Set/get for the object term of this Relationship.

           The common convention for ontologies is to express
           relationships between terms as triples (subject, predicate,
           object).

 Returns : The object term [Bio::Ontology::TermI].
 Args    : 

=cut

=head2 predicate_term

 Title   : predicate_term
 Usage   : $type = $rel->predicate_term();
 Function: Set/get for the predicate of this relationship.

           For a path the predicate (relationship type) is defined as
           the greatest common denominator of all predicates
           (relationship types) encountered along the path. I.e., if
           predicate A is-a predicate B, the greatest common
           denominator for a path containing both predicates A and B is B

 Returns : The predicate term [Bio::Ontology::TermI].
 Args    : 

=cut

=head2 ontology

 Title   : ontology
 Usage   : $ont = $obj->ontology()
 Function: Get the ontology that defined this relationship.
 Example : 
 Returns : an object implementing Bio::Ontology::OntologyI
 Args    : 

See L<Bio::Ontology::OntologyI>.

=cut

1;
