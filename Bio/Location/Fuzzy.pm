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
						 -loc_type => '.');

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
use vars qw(@ISA %FUZZYCODES  %FUZZYPOINTENCODE %FUZZYRANGEENCODE);
use strict;

use Bio::Location::FuzzyLocationI;
use Bio::Location::Simple;

@ISA = qw(Bio::Location::Simple Bio::Location::FuzzyLocationI );

BEGIN {
    %FUZZYCODES = ( 'EXACT' => '..', # Position is 'exact
   # Exact position is unknown, but is within the range specified, ((1.2)..100)
		    'WITHIN' => '.', 
		    # 1^2
		    'BETWEEN' => '^',
		    # <100
		    'BEFORE'  => '<',
		    # >10
		    'AFTER'   => '>');   
   
    %FUZZYPOINTENCODE = ( 
			  '\>(\d+)' => 'AFTER',
			  '\<(\d+)' => 'BEFORE',
			  qr/^(\d+)$/  => 'EXACT',
			  '(\d+)\>' => 'AFTER',
			  '(\d+)\<' => 'BEFORE',
			  '(\d+)\.(\d+)' => 'WITHIN',
			  '(\d+)\^(\d+)' => 'BETWEEN',
		     );
    
    %FUZZYRANGEENCODE  = ( '\.' => 'WITHIN',
			   '\.\.' => 'EXACT',
			   '\^' => 'BETWEEN' );

}

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($loc_type) = $self->_rearrange([qw(LOC_TYPE)], @args);

    $loc_type && $self->loc_type($loc_type);

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
    return $self->SUPER::length() if( !$self->start || !$self->end);
    $self->warn('Length is not valid for a FuzzyLocation'); 
    return 0;
}

=head2 start

  Title   : start
  Usage   : $start = $fuzzy->start();
  Function: get/set start of this range, handling fuzzy_starts
  Returns : an integer representing the start of the range or undef if
            location is fuzzy
  Args    : start location (can be fuzzy string)

=cut

sub start {
    my($self,$value) = @_;
    if( defined $value ) {
	my ($encode, $min,$max) = $self->_fuzzypointencode($value);	
	$self->{'_start_pos_type'} = $encode;
	$self->{'_max_start'} = $self->{'_min_start'} = $self->{'_start'} = undef;
	if( $encode eq 'EXACT' ) {
	    # max and min should be equal to the start value 
	    # if we have an exact point

	    $self->{'_start'} = $value;
	    $self->{'_max_start'} = $value;
	    $self->{'_min_start'} = $value;
	} else { 
	    # so only if we are talking about a WITHIN or BETWEEN 
	    # for both max and min get set to values

	    $self->{'_max_start'} = $max if( $encode ne 'BEFORE' );
	    $self->{'_min_start'} = $min if( $encode ne 'AFTER' );
	}
    }
    # start is undef if we don't know where it is exactly

    return $self->{'_start'};
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
	my ($encode, $min,$max) = $self->_fuzzypointencode($value);		
	$self->{'_end_pos_type'} = $encode;
	# reset all the store values
	$self->{'_max_end'} = $self->{'_min_end'} = $self->{'_end'} = undef;
	if( $encode eq 'EXACT' ) {
	    # max and min should be equal to the end value 
	    # if we have an exact point

	    $self->{'_end'} = $value;
	    $self->{'_max_end'} = $value;
	    $self->{'_min_end'} = $value;
	} else { 
	    # so only if we are talking about a WITHIN or BETWEEN 
	    # for both max and min get set to values
	    $self->{'_max_end'} = $max if( $encode ne 'BEFORE');
	    $self->{'_min_end'} = $min if( $encode ne 'AFTER' );
	}
    }
    # end is undef if we don't know where it is exactly

    return $self->{'_end'};
}

=head2 LocationI methods

=head2 min_start

  Title   : min_start
  Usage   : $min_start = $fuzzy->min_start();
  Function: get the minimum starting point
  Returns : the minimum starting point from the contained sublocations
  Args    : none

=cut

sub min_start {
    my ($self) = @_;    
    return $self->{'_min_start'};
}

=head2 max_start

  Title   : max_start
  Usage   : my $maxstart = $location->max_start();
  Function: Get maximum starting location of feature startpoint  
  Returns : integer or undef if no maximum starting point.
  Args    : none

=cut

sub max_start {
    my ($self) = @_;
    return $self->{'_max_start'};
}

=head2 start_pos_type

  Title   : start_pos_type
  Usage   : my $start_pos_type = $location->start_pos_type();
  Function: Get start position type (ie <,>, ^) 
  Returns : type of position coded as text 
            ('BEFORE', 'AFTER', 'EXACT','WITHIN', 'BETWEEN')
  Args    : none

=cut

sub start_pos_type {
    my ($self) = @_;
    return $self->{'_start_pos_type'};
}

=head2 min_end

  Title   : min_end
  Usage   : my $minend = $location->min_end();
  Function: Get minimum ending location of feature endpoint 
  Returns : integer or undef if no minimum ending point.
  Args    : none

=cut

sub min_end {
    my ($self) = @_;
    return $self->{'_min_end'};
}

=head2 max_end

  Title   : max_end
  Usage   : my $maxend = $location->max_end();
  Function: Get maximum ending location of feature endpoint 
  Returns : integer or undef if no maximum ending point.
  Args    : none

=cut

sub max_end {
    my ($self) = @_;
    return $self->{'_max_end'};
}

=head2 end_pos_type

  Title   : end_pos_type
  Usage   : my $end_pos_type = $location->end_pos_type();
  Function: Get end position type (ie <,>, ^) 
  Returns : type of position coded as text 
            ('BEFORE', 'AFTER', 'EXACT','WITHIN', 'BETWEEN')
  Args    : none

=cut

sub end_pos_type {
    my ($self) = @_;
    return $self->{'_end_pos_type'};
}

=head2 loc_type

  Title   : loc_type
  Usage   : my $location_type = $location->loc_type();
  Function: Get location type encoded as text
  Returns : string ('EXACT', 'WITHIN', 'BETWEEN')
  Args    : none

=cut

sub loc_type {
    my ($self,$value) = @_;
    if( defined $value || ! defined $self->{'_location_type'} ) {
	$value = 'EXACT' unless defined $value;
	$value = uc($value);
	if( $value =~ /\.\./ ) {
	    $value = 'EXACT';
	} elsif( $value =~ /^\.$/ ) {
	    $value = 'WITHIN';
	} elsif( $value =~ /\^/ ) {
	    $value = 'BETWEEN';
	} elsif( $value ne 'EXACT' && $value ne 'WITHIN' && 
		 $value ne 'BETWEEN' ) {
	    $self->throw("Did not specify a valid location type");
	}
	$self->{'_location_type'} = $value;
    }
    return $self->{'_location_type'};
}

=head2 to_FTstring

  Title   : to_FTstring
  Usage   : my $locstr = $location->to_FTstring()
  Function: Get/Set seq_id that location refers to
  Returns : seq_id
  Args    : [optional] seq_id value to set

=cut

sub to_FTstring {
    my ($self) = @_;
    my (%vals) = ( 'start' => $self->start,
		   'min_start' => $self->min_start,
		   'max_start' => $self->max_start,
		   'start_code' => $self->start_pos_type,
		   'end' => $self->end,
		   'min_end' => $self->min_end,
		   'max_end' => $self->max_end,
		   'end_code' => $self->end_pos_type );
    
    my (%strs) = ( 'start' => '',
		   'end'   => '');
    my ($delimiter) = $FUZZYCODES{$self->loc_type};
    # I'm lazy, lets do this in a loop since behaviour will be the same for 
    # start and end
    foreach my $point ( qw(start end) ) {
	if( !defined $vals{$point} ) {	
	    if( (!defined $vals{"min_$point"} || !defined $vals{"max_$point"})
		&& ( $vals{"$point\_code"} eq 'WITHIN' || 
		     $vals{"$point\_code"} eq 'BETWEEN')
		     ) {
		$vals{"min_$point"} = '' unless defined $vals{"min_$point"};
		$vals{"max_$point"} = '' unless defined $vals{"max_$point"};
		
		$self->warn("Fuzzy codes for start are in a strange state, (".
			    join(",", ($vals{"min_$point"}, 
				       $vals{"max_$point"},
				       $vals{"$point\_code"})). ")");
		return '';
	    }
	    
	    if( defined $vals{"min_$point"} ) {
		$strs{$point} .= &_create_point_location($vals{"min_$point"}, 
							 $vals{"$point\_code"});
	    }
	    if( defined $vals{"$point\_code"} ) {
		$strs{$point} .= $FUZZYCODES{$vals{"$point\_code"}}
	    }
	    if( defined $vals{"max_$point"} ) {
		$strs{$point} .= &_create_point_location($vals{"max_$point"}, 
							 $vals{"$point\_code"});
	    }
	    
	} else { 
	    $strs{$point} = $vals{$point};
	}
    }
    my $str = $strs{'start'} . $delimiter . $strs{'end'};
    return $str;
}

sub _create_point_location {
    my ($point,$code) = @_;
    my $str = '';
    return $str if( !defined $point );
    $str = $point;
    return $str if( ! defined $code );
    if( $code eq 'BEFORE') { $str = "<" . $str; }
    elsif( $code eq 'AFTER' ) { $str .= ">"; }    
    return $str;
}

=head2 seq_id

  Title   : seq_id
  Usage   : my $seqid = $location->seq_id();
  Function: Get/Set seq_id that location refers to
  Returns : seq_id
  Args    : [optional] seq_id value to set

=cut

=head2 _fuzzypointencode

  Title   : _fuzzypointencode
  Usage   : $fuzzy->_fuzzypointencode('5>');
  Function: Decode a fuzzy string
  Returns : A two-element array consisting of a integer code of the fuzzy 
            encoding being used, and the integer value of the point.
            A fuzzy code of 0 means 'non-fuzzy', any other code indicates a
            fuzzy location.
          : Returns empty array on fail.
  Args    : fuzzypoint string

=cut

sub _fuzzypointencode {
    my ($self, $string) = @_;
    return () if( !defined $string);
    foreach my $pattern ( keys %FUZZYPOINTENCODE ) {
	if( $string =~ /^\s*$pattern\s*$/ ) {
	    my ($min,$max) = ($1,$2);
	    if( ! defined $max ) {
		$max = $min;
	    }
	    return ($FUZZYPOINTENCODE{$pattern}, $min,$max);
	}
    }
    if( $self->verbose > 1 ) {
	$self->warn("could not find a valid fuzzy encoding for $string");
    }
    return ();
}

1;

