# $Id$

#
# This module is licensed under the same terms as Perl itself. You use,
# modify, and redistribute it under the terms of the Perl Artistic License.
#

=head1 NAME

Bio::IdCollectionI - interface for objects with multiple identifiers

=head1 SYNOPSIS


    # to test this is an identifiable collection object

    $obj->isa("Bio::IdCollectionI") ||
      $obj->throw("$obj does not implement the Bio::IdCollectionI interface");

    # accessors
    @authorities = $obj->id_authorities();
    @ids         = $obj->ids();
    $id          = $obj->ids($authority);

=head1 DESCRIPTION

This interface describes methods expected on objects that have
multiple identifiers, each of which is controlled by a different
authority.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                       - General discussion
  http://bio.perl.org/MailList.html           - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Lincoln Stein

Email lstein@cshl.org

=cut

package Bio::IdCollectionI;
use vars qw(@ISA );
use strict;
use Bio::Root::RootI;


@ISA = qw(Bio::Root::RootI);

=head1 Implementation Specific Functions

These functions are the ones that a specific implementation must
define.

=head2 id_authorities

 Title   : id_authorities
 Usage   : @array    = $obj->id_authorities()
 Function: Return the authorities which have names for this object.
           The authorities can then be used to select ids.

 Returns : An array
 Status  : Virtual

=cut

sub id_authorities {
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 ids

 Title   : ids
 Usage   : @ids    = $obj->ids([$authority1,$authority2...])
 Function: return a list of Bio::LocallyIdentifiableI objects, optionally
           filtered by the list of authorities.

 Returns : A list of Bio::LocallyIdentifiableI objects.
 Status  : Virtual

=cut

sub ids {
   my ($self) = @_;
   my @authorities = @_;
   $self->throw_not_implemented();
}

1;
