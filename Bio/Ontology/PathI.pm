# $Id$
#
# BioPerl module for PathI
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

PathI - Interface for a path between ontology terms

=head1 SYNOPSIS

    # see documentation of methods and an implementation, e.g.,
    # Bio::Ontology::Path

=head1 DESCRIPTION

This is the minimal interface for a path between two terms in
an ontology. Ontology engines may use this.

Essentially this is a very thin extension of the
L<Bio::Ontology::RelationshipI> interface. It basically adds an
attribute distance(). You can think of distance as equal to 0 (parent
== child) or 1 (parent != child) for RelationshipIs.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Hilmar Lapp

Email hlapp at gmx.net

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Ontology::PathI;
use vars qw(@ISA);
use strict;
use Bio::Ontology::RelationshipI;

@ISA = qw( Bio::Ontology::RelationshipI );


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

=head2 parent_term

 Title   : parent_term
 Usage   : $parent = $rel->parent_term();
 Function: Set/get for the parent term of this Path.
 Returns : The parent term [Bio::Ontology::TermI].
 Args    : 

=cut

=head2 child_term

 Title   : child_term
 Usage   : $child = $rel->child_term();
 Function: Set/get for the child term of this Path.
 Returns : The child term [Bio::Ontology::TermI].
 Args    : 

=cut

=head2 relationship_type

 Title   : relationship_type
 Usage   : $type = $rel->relationship_type();
 Function: Set/get for the relationship type of this relationship.

           For a path the relationship type is defined as the greatest
           common denominator of all relationship types encountered
           along the path. I.e., if relationship type A is-a
           relationship type B, the greatest common denominator for a
           path containing both types A and B is B

 Returns : The relationship type [Bio::Ontology::TermI].
 Args    : 

=cut

=head2 ontology

 Title   : ontology
 Usage   : $ont = $obj->ontology()
 Function: Get the ontology that defined this relationship.
 Example : 
 Returns : an object implementing L<Bio::Ontology::OntologyI>
 Args    : 


=cut

1;
