# $Id$
#
# BioPerl module for Bio::Location::SplitLocation
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Location::Split - Implementation of a Location on a Sequence
which has multiple locations (start/end points)

=head1 SYNOPSIS

    my $splitlocation = new Bio::Location::Split();
    $splitlocation->add_sub_Location(new Bio::Location::Simple(-start=>1,
							       -end=>30,
							       -strand=>1));
    $splitlocation->add_sub_Location(new Bio::Location::Simple(-start=>50,
							       -end=>61,
							       -strand=>1));   
    my @sublocs = $splitlocation->sub_Location();

    my $count = 1;
    # print the start/end points of the sub locations
    foreach my $location ( sort { $a->start <=> $b->start } 
			   @sublocs ) {
	printf "sub feature %d [%d..%d]\n", $location->start,$location->end;
        $count++;
    }

=head1 DESCRIPTION

This implementation handles locations which span more than one
start/end location.

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


package Bio::Location::Split;
use vars qw(@ISA);
use strict;

use Bio::Root::RootI;
use Bio::Location::SplitLocationI;
use Bio::Location::Simple;

@ISA = qw(Bio::Location::Simple Bio::Location::SplitLocationI );

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    # initialize
    $self->{'_sublocations'} = [];
    my ( $type, $locations ) = $self->_rearrange([qw(SPLITTYPE 
						     LOCATIONS)], @args);
    if( defined $locations && ref($locations) =~ /array/i ) {
	$self->add_subLocation(@$locations);
    }
    $type = lc ($type);    
    $self->splittype($type || 'join');
    return $self;
}

=head2 sub_Location

 Title   : sub_Location
 Usage   : @locations = $feat->sub_Location();
 Function: Returns an array of LocationI objects
 Returns : An array
 Args    : none

=cut

sub sub_Location {
    my ($self) = @_;
    return @{$self->{'_sublocations'}};
}

=head2 add_sub_Location

 Title   : add_sub_Location
 Usage   : $feat->add_sub_Location(@locationIobjs);
 Function: add an additional sublocation
 Returns : number of current sub locations
 Args    : list of LocationI object(s) to add

=cut

sub add_sub_Location {
    my ($self,@args) = @_;
    my @locs;    
    foreach my $loc ( @args ) {
	if( !ref($loc) || ! $loc->isa('Bio::LocationI') ) {
	    $self->warn("Trying to add $loc as a sub Location but it is not a Bio::LocationI object!");
	    next;
	}
	push @locs, $loc;
    }
    # insert in sorted order, somewhat inefficient
    $self->{'_sublocations'} = [ sort { return 1 unless defined $a->start;
					return -1 unless defined $b->start;
					
					$a->start<=> $b->start          ||
					$a->min_start <=> $b->min_start ||
					$a->max_start <=> $b->max_start ||
					$a->end       <=> $b->end       ||
					$a->min_end   <=> $b->min_end;
				      } 
				 (@locs, @{$self->{'_sublocations'}}) ];
    return scalar @{$self->{'_sublocations'}};
}

=head2 splittype

  Title   : splittype
  Usage   : $splittype = $fuzzy->splittype();
  Function: get/set the split splittype
  Returns : the splittype of split feature (join, order)
  Args    : splittype to set

=cut

sub splittype {
    my ($self, $value) = @_;
    if( defined $value || ! defined $self->{'_splittype'} ) {
	$value = 'join' unless( defined $value );
	$self->{'_splittype'} = lc ($value);
    }
    return $self->{'_splittype'};
}

=head2 LocationI methods

=head2 start

  Title   : start
  Usage   : $start = $location->start();
  Function: get the starting point of the first (sorted) item 
  Returns : integer
  Args    : none

=cut

sub start {
    my ($self,$value) = @_;    
    if( defined $value ) {
	$self->warn("Trying to set the starting point of a split location, that is not possible, try manipulating the sub Locations");
    }
    my @locs = $self->sub_Location();
    return ( @locs ) ? $locs[0]->start : undef;
}

=head2 end

  Title   : end
  Usage   : $end = $location->end();
  Function: get the ending point of the last (sorted) item 
  Returns : integer
  Args    : none

=cut

sub end {
    my ($self,$value) = @_;    
    if( defined $value ) {
	$self->warn("Trying to set the ending point of a split location, that is not possible, try manipulating the sub Locations");
    }
    my @locs = $self->sub_Location();
    return ( @locs ) ? $locs[-1]->end : undef;
}

=head2 min_start

  Title   : min_start
  Usage   : $min_start = $location->min_start();
  Function: get the minimum starting point
  Returns : the minimum starting point from the contained sublocations
  Args    : none

=cut

sub min_start {
    my ($self, $value) = @_;    
    if( defined $value ) {
	$self->warn("Trying to set the minimum starting point of a split location, that is not possible, try manipulating the sub Locations");
    }

    my @locs = $self->sub_Location();
    if( @locs ) {
	return ( defined $locs[0]->min_start ? $locs[0]->min_start :
		 defined $locs[0]->start ? $locs[0]->start : 
		 $locs[0]->max_start); 
    } 
    return undef;
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
    my @locs = $self->sub_Location();
    if( @locs ) {
	return ( defined $locs[0]->max_start ? $locs[0]->max_start :
		 defined $locs[0]->start ? $locs[0]->start : 
		 $locs[0]->min_start); 
    } 
    return undef;
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
    my @locs = $self->sub_Location();
    return ( @locs ) ? $locs[0]->start_pos_type : undef;    
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
    my @locs = $self->sub_Location();
    if( @locs ) {
	return ( defined $locs[-1]->min_end ? $locs[-1]->min_end :
		 defined $locs[-1]->start ? $locs[-1]->start : 
		 $locs[-1]->max_end); 
    } 
    return undef;
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
    my @locs = $self->sub_Location();
    if( @locs ) {
	return ( defined $locs[-1]->max_end ? $locs[-1]->max_end :
		 defined $locs[-1]->end ? $locs[-1]->end : 
		 $locs[-1]->min_end); 
    } 
    return undef;
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
    my @locs = $self->sub_Location();
    return ( @locs ) ? $locs[0]->end_pos_type : undef;    
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
    my @strs;
    foreach my $loc ( $self->sub_Location() ) {
	push @strs, $loc->to_FTstring();
    }
    my $str = sprintf("%s(%s)",$self->splittype, join(",", @strs));
    return $str;
}

=head2 seq_id

  Title   : seq_id
  Usage   : my $seqid = $location->seq_id();
  Function: Get/Set seq_id that location refers to
  Returns : seq_id
  Args    : [optional] seq_id value to set

=cut

# we'll probably need to override the RangeI methods since our locations will
# not be contiguous.

1;
