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
    my $fuzzy_start = $loc->fuzzy_start();
    my $fuzzy_end = $loc->fuzzy_end();
    print $fuzzy_start, $loc->range_fuzzy ? "." : "..",
          $fuzzy_end, "\n";

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
use vars qw(@ISA %FUZZYPOINTENCODE %FUZZYRANGEENCODE);
use strict;

use Bio::LocationI;
use Carp;

@ISA = qw(Bio::LocationI);

BEGIN { 
    %FUZZYPOINTENCODE = ( 
		     '\>(\d+)' => -2,
		     '\<(\d+)' => -1,
		     '(\d+)'  => 0,
		     '(\d+)\>' => 1,
		     '(\d+)\<' => 2,
		     );
    
    %FUZZYRANGEENCODE  = ( '.' => -1,
			   '..' => 0,
			   '^' => 1 );
}

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

=head2 fuzzy_start

  Title   : fuzzy_start
  Usage   : $fuzzystr = $fuzzy->fuzzy_start();
  Function: get/set if start point as a fuzzystring
  Returns : fuzzy startpoint string
  Args    : [optional] fuzzy startpoint string to set

=cut

sub fuzzy_start {
    my ($self, $value) = @_;
    $self->_abstractDeath();
}

=head2 fuzzy_end

  Title   : fuzzy_end
  Usage   : $fuzzystr = $fuzzy->fuzzy_end();
  Function: get/set fuzzy endpoint
  Returns : fuzzy endpoint string
  Args    : [optional] fuzzy endpoint string to set

=cut

sub fuzzy_end {
    my ($self, $value) = @_;
    $self->_abstractDeath();
}

=head2 fuzzy_range

  Title   : fuzzy_range
  Usage   : $status = $fuzzy->fuzzy_range();
  Function: get/set if range is fuzzy (ie 10.20 )
  Returns : true if range is fuzzy, false otherwise
  Args    : optionaly allows the status to be set
          : using $fuzzy->fuzzy_range($value)
=cut

sub fuzzy_range {
    my ($self, $value) = @_;
    $self->_abstractDeath();
}

=head2 _fuzzypointencode

  Title   : _fuzzypointencode
  Usage   : $fuzzy->_fuzzypointencode('5>');
  Function: encode a fuzzy string
  Returns : array of fuzzy encoding and the integer value of the point
          : empty array on fail
  Args    : fuzzypoint string

=cut

sub _fuzzypointencode {
    my ($self, $string) = @_;
    return () if( !defined $string);
    foreach my $pattern ( keys %FUZZYPOINTENCODE ) {
	if( $string =~ /^\s*$pattern\s*$/ ) {
	    return ($FUZZYPOINTENCODE{$pattern}, $1);
	}
    }
    if( $self->verbose > 1 ) {
	$self->warn("could not find a valid fuzzy encoding for $string");
    }
    return ();
}

=head2 _fuzzyrangeencode

  Title   : _fuzzyrange
  Usage   : $fuzzy->_fuzzyrange('.');
  Function: encode a fuzzy range
  Returns : fuzzy range encoding string or undef on fail 
  Args    : fuzzy range string [ '.', '..', '^' ]

=cut

sub _fuzzyrange {
    my ($self, $string) = @_;
    my $encode = $FUZZYRANGEENCODE{$string};
    if( !defined $encode ) {
	$self->warn("could not find a valid fuzzy encoding for $string") if( $self->verbose > 1);
    }
    return $encode;
}

1;

