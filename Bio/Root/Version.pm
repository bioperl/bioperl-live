# $Id$
#
# BioPerl module for Bio::Root::Version
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Root::Version - provide global, distribution-level versioning

=head1 SYNOPSIS

  package Bio::Tools::NiftyFeature;
  require Bio::Root::RootI;


  # later, in client code:
  package main;
  use Bio::Tools::NiftyFeature 3.14;


  ## alternative usage: NiftyFeature defines own $VERSION:
  package Bio::Tools::NiftyFeature;
  my $VERSION = 9.8;

  # later in client code:
  package main;

  # ensure we're using an up-to-date BioPerl distribution
  use Bio::Perl 3.14;

  # NiftyFeature has its own versioning scheme:
  use Bio::Tools::NiftyFeature 9.8;

=head1 DESCRIPTION

This module provides a mechanism by which all other BioPerl modules
can share the same $VERSION, without manually synchronizing each file.

Bio::Root::RootI itself uses this module, so any module that directly
(or indirectly) uses Bio::Root::RootI will get a global $VERSION
variable set if it's not already.


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Aaron Mackey

Email amackey@virginia.edu

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Root::Version;
use strict;
use vars qw($VERSION);

$VERSION = 1.3;

sub import {
    # try to handle multiple levels of inheritance:
    my $i = 0;
    my $pkg = caller($i);
    no strict 'refs';
    while ($pkg) {
	if ($pkg =~ m/^Bio::/o &&
	    not defined ${$pkg . "::VERSION"}) {
	    ${$pkg . "::VERSION"} = $VERSION;
	}
        $pkg = caller(++$i);
    }
}

1;
__END__
