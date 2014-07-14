#
# BioPerl module for Bio::Location::Split
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
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

    my $splitlocation = Bio::Location::Split->new();
    $splitlocation->add_sub_Location(Bio::Location::Simple->new(-start=>1,
							       -end=>30,
							       -strand=>1));
    $splitlocation->add_sub_Location(Bio::Location::Simple->new(-start=>50,
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
start/end location, or and/or lie on different sequences, and can
work with split locations that depend on the specific order of the
sublocations ('join') or don't have a specific order but represent
a feature spanning noncontiguous sublocations ('order', 'bond').

Note that the order in which sublocations are added may be very important,
depending on the specific split location type.  For instance, a 'join'
must have the sublocations added in the order that one expects to
join the sublocations, whereas all other types are sorted based on the
sequence location.

=head1 FEEDBACK

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-AT-bioperl_DOT_org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Location::Split;

# as defined by BSANE 0.03
our @CORBALOCATIONOPERATOR = ('NONE','JOIN', undef, 'ORDER');;

use Bio::Root::Root;

use base qw(Bio::Location::Atomic Bio::Location::SplitLocationI);

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    # initialize
    $self->{'_sublocations'} = [];
    my ( $type, $seqid, $locations ) = 
	$self->_rearrange([qw(SPLITTYPE
                              SEQ_ID
			      LOCATIONS
                              )], @args);
    if( defined $locations && ref($locations) =~ /array/i ) {
	$self->add_sub_Location(@$locations);
    }
    $seqid  && $self->seq_id($seqid);
    $type ||= 'JOIN';
    $type = lc ($type);
    $self->splittype($type);
    return $self;
}

=head2 each_Location

 Title   : each_Location
 Usage   : @locations = $locObject->each_Location($order);
 Function: Conserved function call across Location:: modules - will
           return an array containing the component Location(s) in
           that object, regardless if the calling object is itself a
           single location or one containing sublocations.
 Returns : an array of Bio::LocationI implementing objects
 Args    : Optional sort order to be passed to sub_Location()

=cut

sub each_Location {
    my ($self, $order) = @_;
    my @locs = ();
    foreach my $subloc ($self->sub_Location($order)) {
	# Recursively check to get hierarchical split locations:
	push @locs, $subloc->each_Location($order);
    }
    return @locs;
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
           value of the first argument. The default is no sort.

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
    my @sublocs = defined $self->{'_sublocations'} ? @{$self->{'_sublocations'}} : ();

    # return the array if no ordering requested
    return @sublocs if( ($order == 0) || (! @sublocs) );
    	
    # sort those locations that are on the same sequence as the top (`master')
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
			# Schwartzian transforms for performance boost	  
			@locs = map { $_->[0] }
			sort {
				(defined $a && defined $b) ? $a->[1] <=> $b->[1] :
                $a                         ?  -1                 : 1
				}
			map {
				[$_, (defined $_->start ? $_->start : $_->end)]
				} @locs;;
		} else { # $order == -1
			@locs = map { $_->[0]}
			sort { 
				(defined $a && defined $b) ? $b->[1] <=> $a->[1] :
				$a                         ? -1                  : 1
				}
			map {
				[$_, (defined $_->end ? $_->end : $_->start)]
				} @locs;
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
  Usage   : $splittype = $location->splittype();
  Function: get/set the split splittype
  Returns : the splittype of split feature (join, order)
  Args    : splittype to set

=cut

sub splittype {
    my ($self, $value) = @_;
    if( defined $value || ! defined $self->{'_splittype'} ) {
	$value = 'JOIN' unless( defined $value );
	$self->{'_splittype'} = uc ($value);
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

=head2 guide_strand

  Title   : guide_strand
  Usage   : $str = $loc->guide_strand();
  Function: Get/Set the guide strand.  Of use only if the split type is
            a 'join' (this helps determine the order of sublocation
			retrieval)
  Returns : value of guide strand (1, -1, or undef)
  Args    : new value (-1 or 1, optional)

=cut

sub guide_strand {
    my $self = shift;
    return $self->{'strand'} = shift if @_;

    # Sublocations strand values consistency check to set Guide Strand
    my @subloc_strands;
    foreach my $loc ($self->sub_Location(0)) {
        push @subloc_strands, $loc->strand || 1;
    }
    if ($self->isa('Bio::Location::SplitLocationI')) {
        my $identical   = 0;
        my $first_value = $subloc_strands[0];
        foreach my $strand (@subloc_strands) {
            $identical++ if ($strand == $first_value);
        }

        if ($identical == scalar @subloc_strands) {
            $self->{'strand'} = $first_value;
        }
        else {
            $self->{'strand'} = undef;
        }
    }
    return $self->{'strand'};
}

=head1 LocationI methods

=head2 strand

 Title   : strand
 Usage   : $obj->strand($newval)
 Function: For SplitLocations, setting the strand of the container
           (this object) is a short-cut for setting the strand of all
           sublocations.

           In get-mode, checks if no sub-location is remote, and if
           all have the same strand. If so, it returns that shared
           strand value. Otherwise it returns undef.

 Example : 
 Returns : on get, value of strand if identical between sublocations 
           (-1, 1, or undef)
 Args    : new value (-1 or 1, optional)


=cut

sub strand{
    my ($self,$value) = @_;
    if( defined $value) {
		$self->{'strand'} = $value;
		# propagate to all sublocs
		foreach my $loc ($self->sub_Location(0)) {
			$loc->strand($value);
		}
    } else {
		my ($strand, $lstrand);
		foreach my $loc ($self->sub_Location(0)) {
			# we give up upon any location that's remote or doesn't have
			# the strand specified, or has a differing one set than 
			# previously seen.
			# calling strand() is potentially expensive if the subloc is also
			# a split location, so we cache it
			$lstrand = $loc->strand();
			if((! $lstrand) ||
			   ($strand && ($strand != $lstrand)) ||
			   $loc->is_remote()) {
			$strand = undef;
			last;
			} elsif(! $strand) {
			$strand = $lstrand;
			}
		}
		return $strand;
    }
}

=head2 flip_strand

  Title   : flip_strand
  Usage   : $location->flip_strand();
  Function: Flip-flop a strand to the opposite.  Also sets Split strand
            to be consistent with the sublocation strands
            (1, -1 or undef for mixed strand values)
  Returns : None
  Args    : None

=cut

sub flip_strand {
    my $self = shift;
    my @sublocs;
    my @subloc_strands;

    for my $loc ( $self->sub_Location(0) ) {
        # Atomic "flip_strand" now initialize strand if necessary
        my $new_strand = $loc->flip_strand;

        # Store strand values for later consistency check
        push @sublocs, $loc;
        push @subloc_strands, $new_strand;
    }

    # Sublocations strand values consistency check to set Guide Strand
    if ($self->isa('Bio::Location::SplitLocationI')) {
        my $identical   = 0;
        my $first_value = $subloc_strands[0];
        foreach my $strand (@subloc_strands) {
            $identical++ if ($strand == $first_value);
        }

        if ($identical == scalar @subloc_strands) {
            $self->guide_strand($first_value);
        }
        else {
            # Mixed strand values, must reverse the sublocations order
            $self->guide_strand(undef);
            @{ $self->{_sublocations} } = reverse @sublocs;
        }
    }
}

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
	$self->throw("Trying to set the starting point of a split location, ".
				 "that is not possible, try manipulating the sub Locations");
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
	$self->throw("Trying to set the ending point of a split location, ".
				 "that is not possible, try manipulating the sub Locations");
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
	$self->throw("Trying to set the minimum starting point of a split ".
				 "location, that is not possible, try manipulating the sub Locations");
    }
    my @locs = $self->sub_Location(1);
    return $locs[0]->min_start() if @locs; 
    return;
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
	$self->throw("Trying to set the maximum starting point of a split ".
				 "location, that is not possible, try manipulating the sub Locations");
    }
    my @locs = $self->sub_Location(1);
    return $locs[0]->max_start() if @locs; 
    return;
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
	$self->throw("Trying to set the start_pos_type of a split location, ".
				 "that is not possible, try manipulating the sub Locations");
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
	$self->throw("Trying to set the minimum end point of a split location, ".
				 "that is not possible, try manipulating the sub Locations");
    }
    # reverse sort locations by largest ending to smallest ending
    my @locs = $self->sub_Location(-1);
    return $locs[0]->min_end() if @locs; 
    return;
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
	$self->throw("Trying to set the maximum end point of a split location, ".
				 "that is not possible, try manipulating the sub Locations");
    }
    # reverse sort locations by largest ending to smallest ending
    my @locs = $self->sub_Location(-1);
    return $locs[0]->max_end() if @locs; 
    return;
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
	$self->throw("Trying to set end_pos_type of a split location, ".
				 "that is not possible, try manipulating the sub Locations");
    }
    my @locs = $self->sub_Location();
    return ( @locs ) ? $locs[0]->end_pos_type() : undef;    
}

=head2 length

 Title   : length
 Usage   : $len = $loc->length();
 Function: get the length in the coordinate space this location spans
 Example :
 Returns : an integer
 Args    : none

=cut

sub length {
    my ($self) = @_;
    my $length = 0;
    # Mixed strand values means transplicing (where exons can even
    # be in different chromosomes), so in that case only give the sum
    # of the lengths of the individual segments
    if (! defined $self->guide_strand) {
        for my $loc ( $self->sub_Location(0) ) {
            $length += abs($loc->end - $loc->start) + 1
        }
    }
    else {
        my @sublocs = $self->sub_Location(0);
        my $start   = $sublocs[0]->start;
        my $end     = $sublocs[-1]->end;

        # If Start > ·End, its a possible case of cut by origin
        # location in circular sequences (e.g "join(16..20,1..2)")
        if ($start > $end) {
            # Figure out which segments are located before
            # and which are located after coordinate 1
            # (END_SEQ - 1 - START_SEQ)
            my @end_seq_segments;
            my @start_seq_segments;
            my $switch = 0;
            foreach my $subloc (@sublocs) {
                if ($switch == 0) {
                    if ($subloc->start == 1) {
                        $switch = 1;
                        push @start_seq_segments, $subloc;
                    }
                    else {
                        push @end_seq_segments, $subloc;
                    }
                }
                else {
                    push @start_seq_segments, $subloc;
                }
            }

            # If its a cut by origin location, sum the whole length of each group
            if (scalar @end_seq_segments > 0 and @start_seq_segments > 0) {
                my $end_segments_length   = abs(  $end_seq_segments[0]->start
                                                - $end_seq_segments[-1]->end)
                                                + 1;
                my $start_segments_length = abs(  $start_seq_segments[0]->start
                                                - $start_seq_segments[-1]->end)
                                                + 1;
                $length = $end_segments_length + $start_segments_length;
            }
        }
        else {
            $length = $end - $start + 1;
        }
    }

    # If for some reason nothing worked, fall back to previous behaviour
    if ($length == 0) {
        $length = abs($self->end - $self->start) + 1
    }

    return $length;
}

=head2 seq_id

  Title   : seq_id
  Usage   : my $seqid = $location->seq_id();
  Function: Get/Set seq_id that location refers to

            We override this here in order to propagate to all sublocations
            which are not remote (provided this root is not remote either)
  Returns : seq_id
  Args    : [optional] seq_id value to set


=cut

sub seq_id {
    my $self = shift;

    if(@_ && !$self->is_remote()) {
	foreach my $subloc ($self->sub_Location(0)) {
	    $subloc->seq_id(@_) if !$subloc->is_remote();
	}
    }
    return $self->SUPER::seq_id(@_);
}

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
	my $strand = $self->strand() || 0;
	my $stype = lc($self->splittype());

    if( $strand < 0 ) {
		$self->flip_strand; # this will recursively set the strand
							# to +1 for all the sub locations
    }
    
    foreach my $loc ( $self->sub_Location(0) ) {
		$loc->verbose($self->verbose);
		my $str = $loc->to_FTstring();
		# we only append the remote seq_id if it hasn't been done already
		# by the sub-location (which it should if it knows it's remote)
		# (and of course only if it's necessary)
		if( (! $loc->is_remote) &&
			defined($self->seq_id) && defined($loc->seq_id) &&
			($loc->seq_id ne $self->seq_id) ) {
			$str = sprintf("%s:%s", $loc->seq_id, $str);
		} 
		push @strs, $str;
	}
	$self->flip_strand if $strand < 0;
	my $str;
	if( @strs == 1 ) {
		($str) = @strs;
	} elsif( @strs == 0 ) {
		$self->warn("no Sublocations for this splitloc, so not returning anything\n");
	} else { 
		$str = sprintf("%s(%s)",lc $self->splittype, join(",", @strs));
	}
	if( $strand < 0 ) {  # wrap this in a complement if it was unrolled
		$str = sprintf("%s(%s)",'complement',$str);
	}

    return $str;
}

=head2 valid_Location

 Title   : valid_Location
 Usage   : if ($location->valid_location) {...};
 Function: boolean method to determine whether location is considered valid
           (has minimum requirements for Simple implementation)
 Returns : Boolean value: true if location is valid, false otherwise
 Args    : none

=cut

# we'll probably need to override the RangeI methods since our locations will
# not be contiguous.

1;
