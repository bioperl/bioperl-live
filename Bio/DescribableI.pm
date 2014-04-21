
#
# This module is licensed under the same terms as Perl itself. You use,
# modify, and redistribute it under the terms of the Perl Artistic License.
#

=head1 NAME

Bio::DescribableI - interface for objects with human readable names and descriptions

=head1 SYNOPSIS


    # to test this is a describable object

    $obj->isa("Bio::DescribableI") || 
      $obj->throw("$obj does not implement the Bio::DescribableI interface");

    # accessors

    $name = $obj->display_name();
    $desc = $obj->description();



=head1 DESCRIPTION

This interface describes methods expected on describable objects, ie
ones which have human displayable names and descriptions

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Ewan Birney

Email birney@sanger.ac.uk

=cut

package Bio::DescribableI;
use strict;


use base qw(Bio::Root::RootI);

=head1 Implementation Specific Functions

These functions are the ones that a specific implementation must
define.

=head2 display_name

 Title   : display_name
 Usage   : $string    = $obj->display_name()
 Function: A string which is what should be displayed to the user
           the string should have no spaces (ideally, though a cautious
           user of this interface would not assumme this) and should be
           less than thirty characters (though again, double checking 
           this is a good idea)
 Returns : A scalar
 Status  : Virtual

=cut

sub display_name {
   my ($self) = @_;
   $self->throw_not_implemented();
}


=head2 description

 Title   : description
 Usage   : $string    = $obj->description()
 Function: A text string suitable for displaying to the user a 
           description. This string is likely to have spaces, but
           should not have any newlines or formatting - just plain
           text. The string should not be greater than 255 characters
           and clients can feel justified at truncating strings at 255
           characters for the purposes of display
 Returns : A scalar
 Status  : Virtual

=cut

sub description {
   my ($self) = @_;
   $self->throw_not_implemented();
}

1;
