# $Id$
#
# BioPerl module for Bio::Tree::Statistics
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::Statistics - Calculate certain statistics for a Tree

=head1 SYNOPSIS

  use Bio::Tree::Statistics;


=head1 DESCRIPTION

This should be where Tree statistics are calculated.  It was
previously where statistics from a Coalescent simulation.  Currently
it is empty because we have not added any Tree specific statistic
calculations to this module yet.  We welcome any contributions.


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
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 CONTRIBUTORS

none so far

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tree::Statistics;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;

@ISA = qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tree::Statistics();
 Function: Builds a new Bio::Tree::Statistics object 
 Returns : Bio::Tree::Statistics
 Args    :


=cut


1;
