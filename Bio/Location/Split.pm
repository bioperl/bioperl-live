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

    use Bio::Location::Split;

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
	printf "sub feature %d [%d..%d]\n", 
	       $count, $location->start,$location->end, "\n";
        $count++;
    }

=head1 DESCRIPTION

This implementation handles locations which span more than one
start/end location, or and/or lie on different sequences.

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
 Usage   : @sublocs = $splitloc->sub_Location();
 Function: Returns the array of sublocations making up this compound (split)
           location. Those sublocations referring to the same sequence as
           the root split location will be sorted by start position (forward
           sort) or end position (reverse sort) and come first (before
           those on other sequences).

           The sort order can be optionally specified or suppressed by the
           value of the first argument. The default is a forward sort.

 Returns : an array of Bio::LocationI implementing objects
 Args    : Optionally 1, 0, or -1 for specifying a forward, no, or reverse
           sort order

=cut

sub sub_Location {
    my ($self, $order) = @_;

    $order = 0 unless defined $order;
    if( defined($order) && ($order !~ /^-?\d+$/) ) {
	$self->throw("value $order passed in to sub_Location is $order, an invalid value");
    } 
    $order = 1 if($order > 1);
    $order = -1 if($order < -1);

    my @sublocs = defined $self->{'_sublocations'} ? 
	          @{$self->{'_sublocations'}} : ();

    # return the array if no ordering requested
    return @sublocs if( ($order == 0) || (! @sublocs) );
    
    # sort those locations that are on the sequence as the top (`master')
    # if the top seq is undefined, we take the first defined in a sublocation
    my $seqid = $self->seq_id();
    my $i = 0;
    while((! defined($seqid)) && ($i <= $#sublocs)) {
	$seqid = $sublocs[$i++]->seq_id();
    }
    if((! $self->seq_id()) && $seqid) {
	$self->warn("sorted sublocation array requested but ".
		    "root location doesn't define seq_id ".
		    "(at least one sublocation does!)");
    }
    my @locs = ($seqid ?
		grep { $_->seq_id() eq $seqid; } @sublocs :
		@sublocs);
    if(@locs) {
	if($order == 1) {
	    @locs = sort { $a->start() <=> $b->start() } @locs;
	} else { # $order == -1
	    @locs = sort { $b->end() <=> $a->end() } @locs;
	}
    }
    # push the rest unsorted
    if($seqid) {
	push(@locs, grep { $_->seq_id() ne $seqid; } @sublocs);
    }
    # done!
    return @locs;
}

=head2 add_sub_Location

 Title   : add_sub_Location
 Usage   : $splitloc->add_sub_Location(@locationIobjs);
 Function: add an additional sublocation
 Returns : number of current sub locations
 Args    : list of Bio::LocationI implementing object(s) to add

=cut

sub add_sub_Location {
    my ($self,@args) = @_;
    my @locs;    
    foreach my $loc ( @args ) {
	if( !ref($loc) || ! $loc->isa('Bio::LocationI') ) {
	    $self->throw("Trying to add $loc as a sub Location but it doesn't implement Bio::LocationI!");
	    next;
	}	
	push @{$self->{'_sublocations'}}, $loc;
    }

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

=head2 is_single_sequence

  Title   : is_single_sequence
  Usage   : if($splitloc->is_single_sequence()) {
                print "Location object $splitloc is split ".
                      "but only across a single sequence\n";
	    }
  Function: Determine whether this location is split across a single or
            multiple sequences.

            This implementation ignores (sub-)locations that do not define
            seq_id(). The same holds true for the root location.

  Returns : TRUE if all sublocations lie on the same sequence as the root
            location (feature), and FALSE otherwise.
  Args    : none

=cut

sub is_single_sequence {
    my ($self) = @_;

    my $seqid = $self->seq_id();
    foreach my $loc ($self->sub_Location(0)) {
	$seqid = $loc->seq_id() if(! $seqid);
	if(defined($loc->seq_id()) && ($loc->seq_id() ne $seqid)) {
	    return 0;
	}
    }
    return 1;
}

=head1 LocationI methods

=head2 start

  Title   : start
  Usage   : $start = $location->start();
  Function: get the starting point of the first (sorted) sublocation
  Returns : integer
  Args    : none

=cut

sub start {
    my ($self,$value) = @_;    
    if( defined $value ) {
	$self->throw("Trying to set the starting point of a split location, that is not possible, try manipulating the sub Locations");
    }
    return $self->SUPER::start();
}

=head2 end

  Title   : end
  Usage   : $end = $location->end();
  Function: get the ending point of the last (sorted) sublocation
  Returns : integer
  Args    : none

=cut

sub end {
    my ($self,$value) = @_;    
    if( defined $value ) {
	$self->throw("Trying to set the ending point of a split location, that is not possible, try manipulating the sub Locations");
    }
    return $self->SUPER::end();
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
	$self->throw("Trying to set the minimum starting point of a split location, that is not possible, try manipulating the sub Locations");
    }
    my @locs = $self->sub_Location(1);
    return $locs[0]->min_start() if @locs; 
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
    my ($self,$value) = @_;

    if( defined $value ) {
	$self->throw("Trying to set the maximum starting point of a split location, that is not possible, try manipulating the sub Locations");
    }
    my @locs = $self->sub_Location(1);
    return $locs[0]->max_start() if @locs; 
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
    my ($self,$value) = @_;

    if( defined $value ) {
	$self->throw("Trying to set the start_pos_type of a split location, that is not possible, try manipulating the sub Locations");
    }
    my @locs = $self->sub_Location();
    return ( @locs ) ? $locs[0]->start_pos_type() : undef;    
}

=head2 min_end

  Title   : min_end
  Usage   : my $minend = $location->min_end();
  Function: Get minimum ending location of feature endpoint 
  Returns : integer or undef if no minimum ending point.
  Args    : none

=cut

sub min_end {
    my ($self,$value) = @_;

    if( defined $value ) {
	$self->throw("Trying to set the minimum end point of a split location, that is not possible, try manipulating the sub Locations");
    }
    # reverse sort locations by largest ending to smallest ending
    my @locs = $self->sub_Location(-1);
    return $locs[0]->min_end() if @locs; 
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
    my ($self,$value) = @_;

    if( defined $value ) {
	$self->throw("Trying to set the maximum end point of a split location, that is not possible, try manipulating the sub Locations");
    }
    # reverse sort locations by largest ending to smallest ending
    my @locs = $self->sub_Location(-1);
    return $locs[0]->max_end() if @locs; 
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
    my ($self,$value) = @_;

    if( defined $value ) {
	$self->throw("Trying to set end_pos_type of a split location, that is not possible, try manipulating the sub Locations");
    }
    my @locs = $self->sub_Location();
    return ( @locs ) ? $locs[0]->end_pos_type() : undef;    
}


=head2 seq_id

  Title   : seq_id
  Usage   : my $seqid = $location->seq_id();
  Function: Get/Set seq_id that location refers to
  Returns : seq_id
  Args    : [optional] seq_id value to set

=head2 coordinate_policy

  Title   : coordinate_policy
  Usage   : $policy = $location->coordinate_policy();
            $location->coordinate_policy($mypolicy); # set may not be possible
  Function: Get the coordinate computing policy employed by this object.

            See Bio::Location::CoordinatePolicyI for documentation about
            the policy object and its use.

            The interface *does not* require implementing classes to accept
            setting of a different policy. The implementation provided here
            does, however, allow to do so.

            Implementors of this interface are expected to initialize every
            new instance with a CoordinatePolicyI object. The implementation
            provided here will return a default policy object if none has
            been set yet. To change this default policy object call this
            method as a class method with an appropriate argument. Note that
            in this case only subsequently created Location objects will be
            affected.

  Returns : A Bio::Location::CoordinatePolicyI implementing object.
  Args    : On set, a Bio::Location::CoordinatePolicyI implementing object.

=head2 to_FTstring

  Title   : to_FTstring
  Usage   : my $locstr = $location->to_FTstring()
  Function: returns the FeatureTable string of this location
  Returns : string
  Args    : none

=cut

sub to_FTstring {
    my ($self) = @_;
    my @strs;
    foreach my $loc ( $self->sub_Location() ) {	
	my $str = $loc->to_FTstring();
	if( defined $self->seq_id && 
	    defined $loc->seq_id && 
	    $loc->seq_id ne $self->seq_id ) {
	    $str = sprintf("%s:%s", $loc->seq_id, $str);
	} 
	push @strs, $str;
    }    

    my $str = sprintf("%s(%s)",$self->splittype, join(",", @strs));
    if( $self->strand == -1 ) {
	$str = sprintf("complement(%s)",$str);
    }
    return $str;
}

# we'll probably need to override the RangeI methods since our locations will
# not be contiguous.

1;
