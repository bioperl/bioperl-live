# $Id $
#
# BioPerl module for Bio::Seq::LargeSeqI
#
# Cared for by Albert Vilella <avilella@ebi.ac.uk>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::LargeSeqI - Interface class for sequences that cache their
residues in a temporary file

=head1 SYNOPSIS

 #

=head1 DESCRIPTION

The interface class defines a group of sequence classes that do not
keep their sequence information in memory but store it in a file. This
makes it possible to work with very large files even with limited RAM.

The most important consequence of file caching for sequences is that
you do not want to inspect the sequence unless absolutely
necessary. These sequences typically override the length() method not
to check the sequence.

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

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Albert Vilella

Email avilella@ebi.ac.uk

=head1 CONTRIBUTORS

Heikki Lehvaslaiho, heikki@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::LargeSeqI;
use vars qw(@ISA);
use strict;

use Bio::Root::RootI;

@ISA = qw( Bio::Root::RootI );

1;
