# $Id$
#
# BioPerl module for Bio::LocationI
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::LocationI - Abstract interface of a Location on Sequence

=head1 SYNOPSIS

# get a LocationI somehow

=head1 DESCRIPTION

Descript to follow

=head1 FEEDBACK

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org             - General discussion
  http://bio.perl.org/MailList.html - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::LocationI;
use vars qw(@ISA);
use strict;

use Bio::RangeI;
use Carp;

@ISA = qw(Bio::RangeI);

# utility method Prints out a method like: 
# Abstract method stop defined in interface Bio::LocationI not
# implemented by package You::BadLocation

sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller)[1];
  
  my $msg = "Abstract method '$caller' defined in interface Bio::LocationI but not implemented by package $package";
  if( $self->can('throw') ) {
      $self->throw($msg);
  } else {
      confess($msg);
  }
}

1;

