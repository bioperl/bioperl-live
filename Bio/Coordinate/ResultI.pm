# $Id$
#
# bioperl module for Bio::Coordinate::ResultI
#
# Cared for by Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Coordinate::ResultI - Interface to identify coordinate mapper results

=head1 SYNOPSIS

  # not to be used directly

=head1 DESCRIPTION

ResultI identifies Bio::LocationIs returned by
Bio::Coordinate::MapperI implementing classes from other locations.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                        - General discussion
  http://bioperl.org/wiki/Mailing_lists             - About the mailing lists

=head2 Reporting Bugs

report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki-at-bioperl-dot-org

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Coordinate::ResultI;
use vars qw(@ISA );
use strict;

# Object preamble
use Bio::LocationI;

@ISA = qw(Bio::LocationI);


1;

