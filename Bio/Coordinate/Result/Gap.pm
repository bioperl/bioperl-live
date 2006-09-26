# $Id$
#
# BioPerl module for Bio::Coordinate::Result::Gap
#
# Cared for by Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
#
# Copywright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
#

=head1 NAME

Bio::Coordinate::Result::Gap - Another name for Bio::Location::Simple

=head1 SYNOPSIS

  $loc = new Bio::Coordinate::Result::Gap(-start=>10,
                                          -end=>30,
                                          -strand=>1);

=head1 DESCRIPTION

This is a location object for coordinate mapping results.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Heikki Lehvaslaiho

Email heikki-at-bioperl-dot-org

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Coordinate::Result::Gap;
use strict;


use base qw(Bio::Location::Simple Bio::Coordinate::ResultI);


1;
