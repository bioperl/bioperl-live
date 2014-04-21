#
# BioPerl module for Bio::Map::CytoPosition
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::CytoPosition - Marker class with cytogenetic band storing attributes

=head1 SYNOPSIS

  $m1 = Bio::Map::CytoPosition->new ( '-id' => 'A1',
				       '-value' => '2q1-3'
					     );
  $m2 = Bio::Map::CytoPosition->new ( '-id' => 'A2',
				       '-value' => '2q2'
					     );

  if ($m1->cytorange->overlaps($m2->cytorange)) {
      print "Makers overlap\n";
  }


=head1 DESCRIPTION

CytoPosition is marker (Bio::Map::MarkerI compliant) with a location in a
cytogenetic map. See L<Bio::Map::MarkerI> for more information.

Cytogenetic locations are names of bands visible in stained mitotic
eucaryotic chromosomes. The naming follows strict rules which are
consistant at least in higher vertebates, e.g. mammals. The chromosome
name preceds the band names.

The shorter arm of the chromosome is called 'p' ('petit') and usually
drawn pointing up. The lower arm is called 'q' ('queue'). The bands
are named from the region separting these, a centromere (cen), towards
the tips or telomeric regions (ter) counting from 1 upwards. Depending
of the resolution used the bands are identified with one or more
digit. The first digit determines the major band and subsequent digits
sub bands: p1 band can be divided into subbands p11, p12 and 13 and
p11 can furter be divided into subbands p11.1 and p11.2. The dot after
second digit makes it easier to read the values. A region between ands
is given from the centromere outwards towards the telomere (e.g. 2p2-5
or 3p21-35) or from a band in the p arm to a band in the q arm.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki-at-bioperl-dot-org

=head1 CONTRIBUTORS

Sendu Bala  bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::CytoPosition;

use strict;
use integer;
use Bio::Range;

use base qw(Bio::Map::Position);

=head2 cytorange

 Title   : cytorange
 Usage   : my $range = $obj->cytorange();
 Function:
            Converts cytogenetic location set by value method into
            an integer range. The chromosome number determines the
            "millions" in the values.  Human X and Y chromosome
            symbols are represented by values 100 and 101.

            The localization within chromosomes are converted into
            values between the range of 0 and 200,000:

            pter                    cen                       qter
            |------------------------|-------------------------|
            0                     100,000                   200,000

            The values between -100,000 through 0 for centromere to
            100,000 would have reflected the band numbering better but
            use of positive integers was chosen since the
            transformation is very easy. These values are not metric.

            Each band defines a range in a chromosome. A band string
            is converted into a range by padding it with lower and and
            higher end digits (for q arm: '0' and '9') to the length
            of five. The arm and chromosome values are added to these:
            e.g. 21000 & 21999 (band 21) + 100,000 (q arm) + 2,000,000
            (chromosome 2) => 2q21 : 2,121,000 .. 2,121,999. Note that
            this notation breaks down if there is a band or a subband
            using digit 9 in its name!  This is not the case in human
            karyotype.

            The full algorithm used for bands:

            if arm is 'q' then
               pad char for start is '0', for end '9'
               range is chromosome + 100,000 + padded range start or end
            elsif arm is 'p' then
               pad char for start is '9', for end '0'
               range is chromosome + 100,000 - padded range start or end

 Returns : Bio::Range object or undef
 Args    : none

=cut

sub cytorange {
    my ($self) = @_;
    my ($chr, $r, $band, $band2, $arm, $arm2, $lc, $uc, $lcchar, $ucchar);

    return $r if not defined $self->{_value}; # returns undef
    $self->{_value} =~
	#  -----1-----  --------2---------   -----3-----     -------4-------   ---6---
	m/([XY]|[0-9]+)(cen|qcen|pcen|[pq])?(ter|[.0-9]+)?-?([pq]?(cen|ter)?)?([.0-9]+)?/;
    $self->warn("Not a valid value: ". $self->{_value}), return $r
	if not defined $1 ; # returns undef

    $chr = uc $1;
    $self->chr($chr);

    $chr = 100 if $chr eq 'X';
    $chr = 101 if $chr eq 'Y';
    $chr *= 1000000;

    $r = Bio::Range->new();

    $band = '';
    if (defined $3 ) {
	$2 || $self->throw("$& does not make sense: 'arm' or 'cen' missing");
	$band = $3;
	$band =~ tr/\.//d;
    }
    if (defined $6 ) {
	$arm2 = $4;
	$arm2 = $2 if $4 eq ''; # it is not necessary to repeat the arm [p|q]
	$band2 = $6;
	$band2 =~ tr/\.//d;
    
	#find the correct order
    #print STDERR "-|$&|----2|$2|-----3|$band|---4|$4|--------arm2|$arm2|-------------\n";
	if ($band ne '' and (defined $arm2 and $2 ne $arm2 and $arm2 eq 'q') ) {
	    $lc = 'start'; $lcchar = '9';
	    $uc = 'end'; $ucchar = '9';
	}
	elsif ($band ne 'ter' and $2 ne $arm2 and $arm2 eq 'p') {
	    $lc = 'end'; $lcchar = '9';
	    $uc = 'start'; $ucchar = '9';
	}
	elsif ($band eq 'ter' and  $arm2 eq 'p') {
	    $uc = 'start'; $ucchar = '9';
	} # $2 eq $arm2
	elsif ($arm2 eq 'q') {
	    if (_pad($band, 5, '0') < _pad($band2, 5, '0')) {
		$lc = 'start'; $lcchar = '0';
		$uc = 'end'; $ucchar = '9';
	    } else {
		$lc = 'end'; $lcchar = '9';
		$uc = 'start'; $ucchar = '0';		
	    }
	}
	elsif ($arm2 eq 'p') {
	    if (_pad($band, 5, '0') < _pad($band2, 5, '0')) {
		$lc = 'end'; $lcchar = '0';
		$uc = 'start'; $ucchar = '9';
	    } else {
		$lc = 'start'; $lcchar = '9';
		$uc = 'end'; $ucchar = '0';		
	    }
	}
	else {
	    $self->throw("How did you end up here? $&");
	}

	#print STDERR "-------$arm2--------$band2---------$ucchar--------------\n";
	if ( (defined $arm2 and $arm2 eq 'p') or (defined $arm2 and $arm2 eq 'p') ) {
	    $r->$uc(-(_pad($band2, 5, $ucchar)) + 100000 + $chr );
	    if (defined $3 and $3 eq 'ter') {
		$r->end(200000 + $chr);
	    }
	    elsif ($2 eq 'cen' or $2 eq 'qcen' or $2 eq 'pcen'){
		$r->$lc(100000 + $chr);
	    } 
	    elsif ($2 eq 'q') {
		$r->$lc(_pad($band, 5, $lcchar) + 100000 + $chr );
	    } else {
		$r->$lc(-(_pad($band, 5, $lcchar)) + 100000 + $chr );
	    }
	} else { #if:$arm2=q e.g. 9p22-q32
	    #print STDERR "-------$arm2--------$band2---------$ucchar--------------\n";
	    $r->$uc(_pad($band2, 5, $ucchar) +  100000 + $chr);
	    if ($2 eq 'cen' or $2 eq 'pcen') {
		$r->$lc(100000 + $chr);
	    }
	    elsif ($2 eq 'p') {
		if ($3 eq 'ter') {
		    $r->$lc(200000 + $chr);
		} else {
		    $r->$lc(-(_pad($band, 5, $lcchar)) + 100000 + $chr);
		}
	    } else { #$2.==q
		$r->$lc(_pad($band, 5, $lcchar) + 100000 + $chr);
	    }
	}
    }
    #
    # e.g. 10p22.1-cen
    #
    elsif (defined $4 and $4 ne '') {
	#print STDERR "$4-----$&----\n";
	if ($4 eq 'cen' || $4 eq 'qcen' || $4 eq 'pcen') { # e.g. 10p22.1-cen;
	    # '10pcen-qter' does not really make sense but lets have it in anyway
	    $r->end(100000 + $chr);
	    if ($2 eq 'p') {
		if ($3 eq 'ter') {
		    $r->start($chr);
		} else {
		    $r->start(_pad($band, 5, '9') + $chr);
		}
	    }
	    elsif ($2 eq 'cen') {
		$self->throw("'cen-cen' does not make sense: $&");
	    } else {
		$self->throw("Only order p-cen is valid: $&");
	    }
	}
	elsif ($4 eq 'qter' || $4 eq 'ter') { # e.g. 10p22.1-qter, 1p21-qter, 10pcen-qter, 7q34-qter
	    $r->end(200000 + $chr);
	    if ($2 eq 'p'){
		$r->start(-(_pad($band, 5, '9')) + 100000 + $chr); #??? OK?
	    }
	    elsif ($2 eq 'q') {
		$r->start(_pad($band, 5, '0') + 100000 + $chr);
	    }
	    elsif ($2 eq 'cen' || $2 eq 'qcen' || $2 eq 'pcen' ) {
		$r->start(100000 + $chr);
	    }
	}
	elsif ($4 eq 'pter' ) {
	    #print STDERR "$2,$3--$4-----$&----\n";
	    $r->start( $chr);
	     if ($2 eq 'p'){
		$r->end(-(_pad($band, 5, '0')) + 100000 + $chr);
	    }
	    elsif ($2 eq 'q') {
		$r->end(_pad($band, 5, '9') + 100000 + $chr);
	    }
	    elsif ($2 eq 'cen' || $2 eq 'qcen' || $2 eq 'pcen' ) {
		$r->end(100000 + $chr);
	    }
	} else { # -p or -q at the end of the range
	    $self->throw("lone '$4' in $& does not make sense");
	}
    }
    #
    #  e.g 10p22.1, 10pter
    #
    elsif (defined $3 ) {
	if ($2 eq 'p') {
	    if ($3 eq 'ter') { # e.g. 10pter
		$r = Bio::Range->new('-start' => $chr,
				    '-end' => $chr,
				    );
	    } else { # e.g 10p22.1
		$r = Bio::Range->new('-start' => -(_pad($band, 5, '9')) + 100000 + $chr,
				    '-end' => -(_pad($band, 5, '0')) + 100000 + $chr,
				    );
	    }
	} elsif ($2 eq 'q') {
	    if ($3 eq 'ter') { # e.g. 10qter
		$r = Bio::Range->new('-start' => 200000 + $chr,
				    '-end' => 200000 + $chr,
				    );
	    } else { # e.g 10q22.1
		$r = Bio::Range->new('-start' => _pad($band, 5, '0') + 100000 + $chr,
				    '-end' => _pad($band, 5, '9') + 100000 + $chr,
				    );
	    }
	} else { # e.g. 10qcen1.1 !
	    $self->throw("'cen' in $& does not make sense");
	}
    }
    #
    # e.g. 10p
    #
    elsif (defined $2 ) { # e.g. 10p
	if ($2 eq 'p' ) {
	    $r = Bio::Range->new('-start' => $chr,
				'-end' => 100000  + $chr
				);
	}
	elsif ($2 eq 'q' )  {
	    $r = Bio::Range->new('-start' => 100000 + $chr,
				'-end' => 200000 + $chr
				);
	} else { # $2 eq 'cen' || 'qcen'
	    $r = Bio::Range->new('-start' => 100000 + $chr,
				'-end' => 100000 + $chr
				);
	}
    }
    #
    # chr only, e.g. X
    #
    else {
	$r = Bio::Range->new('-start' => $chr,
			    '-end' => 200000 + $chr
			    );
    }
    
    if ($r) {
        $self->start($r->start);
        $self->end($r->end);
    }
    return $r;
}


sub _pad {
    my ($string, $len, $pad_char) = @_;
    __PACKAGE__->throw("function _pad needs a positive integer length, not [$len]") 
	unless $len =~ /^\+?\d+$/;
    __PACKAGE__->throw("function _pad needs a single character pad_char, not [$pad_char]") 
	unless length $pad_char == 1;
    $string ||= '';
    return $string . $pad_char x ( $len - length( $string ) );
}

=head2 range2value

 Title   : range2value
 Usage   : my $value = $obj->range2value($range);
 Function: Sets and returns the value string based on start and end values of
           the Bio::Range object passes as an argument.
 Returns : string or false
 Args    : Bio::Range object

=cut

sub range2value {
    my ($self,$value) = @_;
    if( defined $value) {
	if( ! $value->isa('Bio::Range') ) {
	    $self->throw("Is not a Bio::Range object but a [$value]");
	    return;
	}
	if( ! $value->start ) {
	    $self->throw("Start is not defined in [$value]");
	    return;
	}
	if( ! $value->end ) {
	    $self->throw("End is not defined in [$value]");
	    return;
	}
	if( $value->start < 100000 ) {
	    $self->throw("Start value has to be in millions, not ". $value->start);
	    return;
	}
	if( $value->end < 100000 ) {
	    $self->throw("End value has to be in millions, not ". $value->end);
	    return;
	}

	my ($chr, $arm, $band) = $value->start =~ /(\d+)(\d)(\d{5})/;	
	my ($chr2, $arm2, $band2) = $value->end =~ /(\d+)(\d)(\d{5})/;	

	my ($chrS, $armS, $bandS, $arm2S, $band2S, $sep) = ('', '', '', '', '', '' );
      LOC: {
	  #
	  # chromosome
	  #
	  if ($chr == 100) {
	      $chrS = 'X';
	  }
	  elsif ($chr == 100) {
	      $chrS = 'Y';
	  } else {
	      $chrS = $chr;
	  }
	  last LOC if  $arm == 0 and $arm2 == 2 and $band == 0 and $band2 == 0 ;
	  #
	  # arm
	  #
	  if ($arm == $arm2 ) {
	      if ($arm == 0) {
		  $armS = 'p';
		  #$armS = 'pter' if $band == 0 and $band2 == 0;
		  $bandS = 'ter' if $band == 0;
		  #$arm2S = 'p'; #?
	      }
	      elsif ($arm == 2) {
		  $armS = 'q';
		  $bandS = 'ter' if $band == 0;
	      } else {
		  $armS = 'q';
		  #$arm2S = 'q'; #?
		  $armS = 'cen',  if $band == 0;# and $band2 == 0;
	      }
	  } else {
	      if ($arm == 0) {
		  $armS = 'p';
		  $arm2S = 'q';
		  $arm2S = '' if $band == 0 and $band2 == 0;
	      } else {
		  $armS = 'q';
		  $arm2S = 'p';
		  $arm2S = '' if $arm2 == 2 and $band == 0 and $band2 == 0;
	      }
	  }
	  last LOC if $band == $band2 ;
	  my $c;
	  #
	  # first band (ter is hadled with the arm)
	  #
	  if ($bandS ne 'ter') {
	      if ($armS eq 'p') {
		  $band = 100000 - $band;
		  $c = '9';
	      } else {
		  $c = '0';
	      }
	      $band =~ s/$c+$//; 
	      $bandS = $band;
	      $bandS = substr($band, 0, 2). '.'. substr($band, 2) if length $band > 2;
	  }
	  last LOC unless $band2;
	  #
	  # second band
	  #
	  if ($arm2 == 0) {
	      $arm2S = 'p';
	      $band2 = 100000 - $band2;
	      $c = '0';
	  } else { # 1 or 2
	      $arm2S = 'q';
	      $c = '9';
	  }
	  if ($band2 == 0) {
	      if ($arm2 == 1) {
		  $arm2S = 'p';
		  $band2S = 'cen';
	      } else {
		  $band2S = 'ter';
	      }
	      last LOC;
	  }
	  last LOC if $band eq $band2 and $arm == $arm2;

	  $band2 =~ s/$c+$//; 
	  $band2S = $band2;
	  $band2S = substr($band2, 0, 2). '.'. substr($band2, 2) if length $band2 > 2;

      } # end of LOC:

	if ($armS eq 'p' and $arm2S eq 'p') {
	    my $tmp = $band2S;
	    $band2S = $bandS;
	    $bandS = $tmp;
	}
	$band2S = '' if $bandS eq $band2S ;
	$armS = '' if $bandS eq 'cen';
	$arm2S = '' if $armS eq $arm2S and $band2S ne 'ter';
	$sep = '-' if $arm2S || $band2S;
	$self->value( $chrS. $armS. $bandS. $sep. $arm2S. $band2S);
    }
   return $self->value;
}

=head2 value

 Title   : value
 Usage   : my $pos = $position->value;
 Function: Get/Set the value for this postion
 Returns : scalar, value
 Args    : none to get, OR scalar to set

=cut

sub value {
   my ($self,$value) = @_;
   if( defined $value ) {
       $self->{'_value'} = $value;
       $self->cytorange;
   }
   return $self->{'_value'};
}

=head2 numeric

 Title   : numeric
 Usage   : my $num = $position->numeric;
 Function: Read-only method that is guarantied to return a numeric 
           representation of the start of this position.
 Returns : int (the start of the range)
 Args    : optional Bio::RangeI object 

=cut

sub numeric {
   my $self = shift;
   return $self->start(@_);
}

=head2 chr

 Title   : chr
 Usage   : my $mychr = $position->chr();
 Function: Get/Set method for the chromosome string of the location.
 Returns : chromosome value
 Args    : none to get, OR scalar to set

=cut

sub chr {
   my ($self,$chr) = @_;
   if( defined $chr ) {
       $self->{'_chr'} = $chr;
   }
   return $self->{'_chr'};
}


1;
