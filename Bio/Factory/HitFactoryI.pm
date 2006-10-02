#-----------------------------------------------------------------
# $Id$
#
# BioPerl module for Bio::Factory::HitFactoryI
#
# Cared for by Steve Chervitz <sac@bioperl.org>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

# POD documentation - main docs before the code

=head1 NAME

Bio::Factory::HitFactoryI - Interface for an object that builds Bio::Search::Hit::HitI objects

=head1 SYNOPSIS

To be completed.

=head1 DESCRIPTION

To be completed.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR 

Steve Chervitz E<lt>sac@bioperl.orgE<gt>

See L<the FEEDBACK section | FEEDBACK> for where to send bug reports and comments.

=head1 COPYRIGHT

Copyright (c) 2001 Steve Chervitz. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object methods. 

=cut

#'

package Bio::Factory::HitFactoryI;

use strict;
use Bio::Factory::ObjectFactoryI;


use base qw(Bio::Root::RootI);

=head2 create_hit

 Title   : create_hit
 Usage   : $hit = $factory->create_hit( %params );
 Function: Creates a new Bio::Search::Hit::HitI object.
 Returns : An object that implements the Bio::Search::Hit::HitI interface
 Args    : Named parameters (to be defined)

=cut

sub create_hit {
    my ($self, @args) = @_;
    $self->throw_not_implemented;
}


1;
