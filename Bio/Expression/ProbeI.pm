# $Id$
# BioPerl module for Bio::Expression::ProbeI
#
# Copyright Allen Day <allenday@ucla.edu>, Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expression::ProbeI - an interface class for DNA/RNA probes

=head1 SYNOPSIS

Do not use this module directly

=head1 DESCRIPTION

This provides a standard bioperl interface class for representing
DNA and RNA probes.  It cannot be instantiated directly, but serves
as an abstract base class for implementors.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org            - General discussion
  http://bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR

Allen Day E<lt>allenday@ucla.eduE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::Expression::ProbeI;

use strict;
use Bio::Root::Root;

use base qw(Bio::Root::Root);
use vars qw($DEBUG);

=head2 sequence()

  Title   : sequence
  Usage   : $seq = $probe->sequence()
  Function: get/set the sequence of the probe
  Returns : A Bio::Seq object or undef
  Args    : a new Bio::Seq object (optional)

=cut

sub sequence {
  shift->throw_not_implemented();
}

=head2 length()

  Title   : length
  Usage   : $len = $probe->length()
  Function: get/set the length of the probe sequence
  Returns : An integer
  Args    : a new integer (optional)
  Comments: this should probably defer to the
            Bio::Seq object returned by sequence(),
            but to allow ProbeI objects to be
            lightweight, we also have this method.

=cut

sub length {
  my($self,$arg) = @_;
  if($arg){
    $self->throw(__PACKAGE__ . "::value only accepts integer values") unless $arg =~ /^[\d]+$/;
    $self->{length} = $arg;
  }
  return $self->{length} || 0;
}

=head2 value()

  Title   : value
  Usage   : $val = $probe->value()
  Function: get/set the probe's observed value
  Returns : A numeric value
  Args    : a new numeric value (optional)

=cut

sub value {
  my($self,$arg) = @_;
  if($arg){
    $self->throw(__PACKAGE__ . "::value only accepts numeric values") unless $arg =~ /^[\d.]+$/;
    $self->{value} = $arg;
  }
  return $self->{value} || 0;
}

=head2 value_units()

  Title   : value_units
  Usage   : $units = $probe->units()
  Function: get/set the units of the probe's observed value
  Returns : A string or undef
  Args    : a new string (optional)

=cut

sub value_units {
  my($self,$arg) = @_;
  $self->{value_units} = $arg if defined $arg;
  return $self->{value_units};
}

=head2 standard_deviation()

  Title   : standard_deviation
  Usage   : $std_dev = $probe->standard_deviation()
  Function: get/set the probe's standard deviation of value()
  Returns : A numeric value
  Args    : a new numeric value (optional)
  Comments: no calculation is done here

=cut

sub standard_deviation {
  my($self,$arg) = @_;
  if($arg){
    $self->throw(__PACKAGE__ . "::standard_deviation only accepts numeric values") unless $arg =~ /^[\d.]+$/;
    $self->{standard_deviation} = $arg;
  }
  return $self->{deviation} || 0;
}

=head2 sample_count()

  Title   : sample_count
  Usage   : $sample_count = $probe->sample_count()
  Function: get/set the number of samples used to calculate
            value()
  Returns : An integer
  Args    : a new integer (optional)

=cut

sub sample_count {
  my($self,$arg) = @_;
  if($arg){
    $self->throw(__PACKAGE__ . "::sample_count only accepts integers") unless $arg =~ /^[\d]+$/;
    $self->{sample_count} = $arg;
  }
  return $self->{sample_count} || 0;
}

1;
