# $Id$
#
# BioPerl module for Bio::Location::FuzzyLocationI
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::FuzzyLocationI - Abstract interface of a Location on a Sequence
which has unclear start/end location

=head1 SYNOPSIS

    # Get a FuzzyLocationI object somehow
    # methods have yet to be defined

=head1 DESCRIPTION

This interface encapsulates the necessary methods for representing a
Fuzzy Location, one that does not have clear start and/or end points.
This will initially serve to handle features from Genbank/EMBL feature
tables that are written as 1^100 meaning between bases 1 and 100 or
<100..300 meaning it starts somewhere before 100.  Advanced
implementations of this interface may be able to handle the necessary
logic of overlaps/intersection/contains/union, but initially this will
be just a holder for the Genbank/EMBL fuzzy location parsing and
producing.

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

=head1 AUTHOR - Jason Stajich

Email jason@chg.mc.duke.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::Location::FuzzyLocationI;
use vars qw(@ISA);
use strict;

use Bio::LocationI;
use Carp;

@ISA = qw(Bio::LocationI);

# utility method Prints out a method like: 
# Abstract method stop defined in interface Bio::LocationI not
# implemented by package You::BadLocation

sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller)[1];
  
  my $msg = "Abstract method '$caller' defined in interface Bio::ComplexLocationI but not implemented by package $package";
  if( $self->can('throw') ) {
      $self->throw($msg);
  } else {
      confess($msg);
  }
}

=head2

  Title   : start_fuzzy
  Usage   : $status = $fuzzy->start_fuzzy();
  Function: get/set if start point is fuzzy
  Returns : true if start point is fuzzy, false otherwise
  Args    : optionaly allows the status to be set
          : using $fuzzy->start_fuzzy($value)

=cut

sub start_fuzzy {
    my ($self, $value) = @_;
    $self->_abstractDeath();
}

=head2

  Title   : end_fuzzy
  Usage   : $status = $fuzzy->end_fuzzy();
  Function: get/set if end point is fuzzy
  Returns : true if end point is fuzzy, false otherwise
  Args    : optionaly allows the status to be set
          : using $fuzzy->end_fuzzy($value)

=cut

sub end_fuzzy {
    my ($self, $value) = @_;
    $self->_abstractDeath();
}

=head2

  Title   : range_fuzzy
  Usage   : $status = $fuzzy->range_fuzzy();
  Function: get/set if range is fuzzy (ie 10.20 )
  Returns : true if range is fuzzy, false otherwise
  Args    : optionaly allows the status to be set
          : using $fuzzy->range_fuzzy($value)

=cut

sub range_fuzzy {
    my ($self, $value) = @_;
    $self->_abstractDeath();
}

1;

