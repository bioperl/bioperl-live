#
# BioPerl module for RelationshipI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Peter Dimitrov <dimitrov@gnf.org>
#
# (c) Peter Dimitrov
# (c) GNF, Genomics Institute of the Novartis Research Foundation, 2002.
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

Bio::Ontology::RelationshipI - Interface for a relationship between ontology terms

=head1 SYNOPSIS

    # see documentation of methods and an implementation, e.g.,
    # Bio::Ontology::Relationship

=head1 DESCRIPTION

This is the minimal interface for a relationship between two terms in
an ontology. Ontology engines will use this.

The terminology we use here is the one commonly used for ontologies,
namely the triple of (subject, predicate, object), which in addition
is scoped in a namespace (ontology). It is called triple because it is
a tuple of three ontology terms.

There are other terminologies in use for expressing relationships. For
those who it helps to better understand the concept, the triple of
(child, relationship type, parent) would be equivalent to the
terminology chosen here, disregarding the question whether the notion
of parent and child is sensible in the context of the relationship
type or not. Especially in the case of ontologies with a wide variety
of predicates the parent/child terminology and similar ones can
quickly become ambiguous (e.g., A synthesises B), meaningless (e.g., A
binds B), or even conflicting (e.g., A is-parent-of B), and are
therefore strongly discouraged.

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

=head1 AUTHOR - Peter Dimitrov

Email dimitrov@gnf.org

=head1 CONTRIBUTORS

 Hilmar Lapp, email: hlapp at gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Ontology::RelationshipI;
use strict;

use base qw(Bio::Root::RootI);

=head2 identifier

 Title   : identifier
 Usage   : print $rel->identifier();
 Function: Set/get for the identifier of this Relationship.

           Note that this may not necessarily be used by a particular
           ontology.

 Returns : The identifier [scalar].
 Args    : 

=cut

sub identifier{
    shift->throw_not_implemented();
}

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

sub subject_term{
    shift->throw_not_implemented();
}

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

sub object_term{
    shift->throw_not_implemented();
}

=head2 predicate_term

 Title   : predicate_term
 Usage   : $type = $rel->predicate_term();
 Function: Set/get for the relationship type of this relationship.

           The common convention for ontologies is to express
           relationships between terms as triples (subject, predicate,
           object).

 Returns : The relationship type [Bio::Ontology::TermI].
 Args    : 

=cut

sub predicate_term{
    shift->throw_not_implemented();
}

=head2 ontology

 Title   : ontology
 Usage   : $ont = $obj->ontology()
 Function: Get the ontology that defined (is the scope for) this
           relationship.
 Example : 
 Returns : an object implementing Bio::Ontology::OntologyI
 Args    : 

See L<Bio::Ontology::OntologyI>.

=cut

sub ontology{
    shift->throw_not_implemented();
}

1;
