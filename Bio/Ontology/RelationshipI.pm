# $Id$
#
# BioPerl module for RelationshipI
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

RelationshipI - Interface for a relationship between ontology terms

=head1 SYNOPSIS

    # see documentation of methods and an implementation, e.g.,
    # Bio::Ontology::Relationship

=head1 DESCRIPTION

This is the minimal interface for a relationship between two terms in
an ontology. Ontology engines will use this.

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

=head1 AUTHOR - Peter Dimitrov

Email dimitrov@gnf.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Ontology::RelationshipI;
use vars qw(@ISA);
use strict;
use Bio::Root::RootI;

@ISA = qw( Bio::Root::RootI );

=head2 identifier

 Title   : identifier
 Usage   : print $rel->identifier();
 Function: Set/get for the identifier of this Relationship.
 Returns : The identifier [scalar].
 Args    : 

=cut

sub identifier{
    shift->throw_not_implemented();
}

=head2 parent_term

 Title   : parent_term
 Usage   : $parent = $rel->parent_term();
 Function: Set/get for the parent term of this Relationship.
 Returns : The parent term [Bio::Ontology::TermI].
 Args    : 

=cut

sub parent_term{
    shift->throw_not_implemented();
}

=head2 child_term

 Title   : child_term
 Usage   : $child = $rel->child_term();
 Function: Set/get for the child term of this Relationship.
 Returns : The child term [Bio::Ontology::TermI].
 Args    : 

=cut

sub child_term{
    shift->throw_not_implemented();
}

=head2 relationship_type

 Title   : relationship_type
 Usage   : $type = $rel->relationship_type();
 Function: Set/get for the relationship type of this relationship.
 Returns : The relationship type [Bio::Ontology::TermI].
 Args    : 

=cut

sub relationship_type{
    shift->throw_not_implemented();
}

=head2 ontology

 Title   : ontology
 Usage   : $ont = $obj->ontology()
 Function: Get the ontology that defined this relationship.
 Example : 
 Returns : an object implementing L<Bio::Ontology::OntologyI>
 Args    : 


=cut

sub ontology{
    shift->throw_not_implemented();
}

1;
