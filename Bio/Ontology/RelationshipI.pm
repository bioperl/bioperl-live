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

RelationshipI - DESCRIPTION of Interface

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the interface here

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
  http://bioperl.org/bioperl-bugs/

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
use Carp;
use Bio::Root::Root;

@ISA = qw( Bio::Root::Root );

=head2 identifier

 Title   : identifier
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub identifier{
  my ($self) = @_;

  $self->throw("Abstract method identifier implementing class did not provide method");
}

=head2 parent_term

 Title   : parent_term
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub parent_term{
  my ($self) = @_;

  $self->throw("Abstract method parent_term implementing class did not provide method");


}

=head2 child_term

 Title   : child_term
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub child_term{
  my ($self) = @_;

  $self->throw("Abstract method child_term implementing class did not provide method");
}

=head2 relationship_type

 Title   : relationship_type
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub relationship_type{
  my ($self) = @_;

  $self->throw("Abstract method relationship_type implementing class did not provide method");
}


1;
