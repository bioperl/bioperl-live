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

    my $fuzzylocation = new Bio::Location::Fuzzy(-start => 30,
						 -end   => 90,
						 -startfuzzy => -1
						 -rangefuzzy => 1);

    my $fuzzy_start = $loc->fuzzy_string($loc->start, $loc->start_fuzzy);
    my $fuzzy_end = $loc->fuzzy_string($loc->end, $loc->end_fuzzy);
    print $fuzzy_start, $loc->range_fuzzy ? "." : "..",
          $fuzzy_end, "\n";


=head1 DESCRIPTION

This module implements the necessary methods for representing a
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


package Bio::Location::Fuzzy;
use vars qw(@ISA);
use strict;

use Bio::Location::FuzzyLocationI;
use Bio::Location::Simple;

@ISA = qw(Bio::Location::Simple Bio::Location::FuzzyLocationI );

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($startfuzzy, $endfuzzy,
	$rangefuzzy) = $self->_rearrange([qw(STARTFUZZY
					     ENDFUZZY 
					     RANGEFUZZY)],
						    @args);
    
    $startfuzzy = 0 unless( defined($startfuzzy)); 
    $endfuzzy = 0 unless( defined($endfuzzy)); 
    $rangefuzzy = 0 unless( defined($rangefuzzy)); 

    $self->start_fuzzy($startfuzzy);
    $self->end_fuzzy($endfuzzy);
    $self->range_fuzzy($rangefuzzy);

    return $self;
}

=head2

  Title   : length
  Usage   : $length = $fuzzy->length();
  Function: get/set the length of this range
  Returns : the length of this range
  Args    : none

=cut

sub length {
    my($self) = @_;
    return $self->SUPER::length() if( !$self->start_fuzzy && !$self->end_fuzzy);
    $self->warn('Length is not valid for a FuzzyLocation'); 
    return 0;
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
    if( defined $value ) {
	$self->{'_startfuzzy'} = $value;
    }
    return $self->{'_startfuzzy'};
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
    if( defined $value ) {
	$self->{'_endfuzzy'} = $value;
    }
    return $self->{'_endfuzzy'};
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
    if( defined $value ) {
	$self->{'_rangefuzzy'} = $value;
    }
    return $self->{'_rangefuzzy'};
}

1;

