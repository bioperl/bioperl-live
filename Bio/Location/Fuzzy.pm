# $Id$
#
# BioPerl module for Bio::Location::Fuzzy
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Location::Fuzzy - Implementation of a Location on a Sequence
which has unclear start and/or end locations

=head1 SYNOPSIS

    my $fuzzylocation = new Bio::Location::Fuzzy(-start => '<30',
						 -end   => 90,
						 -rangefuzzy => '.');

    my $fuzzy_start = $loc->fuzzy_start();
    my $fuzzy_end = $loc->fuzzy_end();
    print $fuzzy_start, $loc->range_fuzzy ? "." : "..",
          $fuzzy_end, "\n";

=head1 DESCRIPTION

This module implements the necessary methods for representing a Fuzzy
Location, one that does not have clear start and/or end points.  This
will initially serve to handle features from Genbank/EMBL feature
tables that are written as 1^100 meaning between bases 1 and 100 or
<100..300 meaning the feature starts somewhere before 100.  Advanced
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


package Bio::Location::Fuzzy;
use vars qw(@ISA);
use strict;

use Bio::Location::FuzzyLocationI;
use Bio::Location::Simple;

@ISA = qw(Bio::Location::Simple Bio::Location::FuzzyLocationI );

sub new {
    my ($class, @args) = @_;

    my $self = $class->SUPER::new(@args);
    my ($rangefuzzy) = $self->_rearrange([qw(RANGEFUZZY)],
					 @args);
    
    $self->fuzzy_range( ($rangefuzzy) ? $rangefuzzy : '..');
    return $self;
}

=head2 length

  Title   : length
  Usage   : $length = $fuzzy->length();
  Function: get/set the length of this range
  Returns : the length of this range
  Args    : none

=cut

sub length {
    my($self) = @_;
    return $self->SUPER::length() if( !$self->fuzzy_start && !$self->fuzzy_end);
    $self->warn('Length is not valid for a FuzzyLocation'); 
    return 0;
}

=head2 start

  Title   : start
  Usage   : $start = $fuzzy->start();
  Function: get/set start of this range, handling fuzzy_starts
  Returns : and integer representing the start of the range
  Args    : start location (can be fuzzy string)

=cut

sub start {
    my($self,$value) = @_;

    if( defined $value ) {
	my ($encode, $start) = $self->_fuzzypoint($value);	
	if( defined $encode ) {
	    $self->fuzzy_start($value);
	    $value = $start;
	}
    }
    return $self->SUPER::start($value);
}

=head2 end

  Title   : end
  Usage   : $end = $fuzzy->end();
  Function: get/set end of this range, handling fuzzy_ends
  Returns : and integer representing the end of the range
  Args    : end location (can be fuzzy string)

=cut

sub end {
    my($self,$value) = @_;

    if( defined $value ) {
	my ($encode, $end) = $self->_fuzzypoint($value);	
	if( defined $encode ) {
	    $self->fuzzy_end($value);
	    $value = $end;
	}
    }
    return $self->SUPER::end($value);
}

=head2 fuzzy_start

  Title   : fuzzy_start
  Usage   : $status = $fuzzy->fuzzy_start();
  Function: get/set fuzzy startpoint
  Returns : fuzzy start string
  Args    : 

=cut

sub fuzzy_start {
    my ($self, $value) = @_;
    if( defined $value ) {	
	my ($encode, $end) = $self->_fuzzypoint($value);
	if( !defined $encode ) {
	    $self->throw("Trying to set fuzzy_end to an invalid value ($value)");
	}
	$self->{'_fuzzystart'} = $value;
    }
    return $self->{'_fuzzystart'} || $self->start;
}

=head2 fuzzy_end

  Title   : fuzzy_end
  Usage   : $status = $fuzzy->fuzzy_end();
  Function: get/set if end point is fuzzy
  Returns : true if end point is fuzzy, false otherwise
  Args    : optionaly allows the status to be set
          : using $fuzzy->fuzzy_end($value)

=cut

sub fuzzy_end {
    my ($self, $value) = @_;
    if( defined $value ) {
	my ($encode, $end) = $self->_fuzzypoint($value);
	if( !defined $encode ) {
	    $self->throw("Trying to set fuzzy_end to an invalid value ($value)");
	}
	$self->{'_fuzzyend'} = $value;
    }
    return $self->{'_fuzzyend'} || $self->end;
}

=head2 fuzzy_range

  Title   : fuzzy_range
  Usage   : $status = $fuzzy->fuzzy_range();
  Function: get/set range delimiter
  Returns : range delimiter
  Args    : optionaly allows the delimiter to be set
          : using $fuzzy->fuzzy_range($value)

=cut

sub fuzzy_range {
    my ($self, $value) = @_;
    if( defined $value ) {
	my ($encoded) = $self->_fuzzyrange($value);
	if( !defined $encoded ) {
	    $self->throw("Tried to set fuzzy_range to an invalid value ($value)");
	}
	$self->{'_fuzzyrange'} = $value;
    }
    return $self->{'_fuzzyrange'};
}

1;

