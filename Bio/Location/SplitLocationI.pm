#
# BioPerl module for Bio::Location::SplitLocationI
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Location::SplitLocationI - Abstract interface of a Location on a Sequence
which has multiple locations (start/end points)

=head1 SYNOPSIS

  # get a SplitLocationI somehow
    print $splitlocation->start, "..", $splitlocation->end, "\n";
    my @sublocs = $splitlocation->sub_Location();

    my $count = 1;
    # print the start/end points of the sub locations
    foreach my $location ( sort { $a->start <=> $b->start }  @sublocs ) {
		printf "sub feature %d [%d..%d]\n", $location->start,$location->end;
        $count++;
    }

=head1 DESCRIPTION

This interface encapsulates the necessary methods for representing the
location of a sequence feature that has more that just a single
start/end pair.  Some examples of this are the annotated exons in a
gene or the annotated CDS in a sequence file.

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

Email jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::Location::SplitLocationI;
use strict;

use Carp;

use base qw(Bio::LocationI);


=head2 sub_Location

 Title   : sub_Location
 Usage   : @locations = $feat->sub_Location();
 Function: Returns an array of LocationI objects
 Returns : An array
 Args    : none

=cut

sub sub_Location {
    my ($self,@args) = @_;
    $self->throw_not_implemented();
}

=head2 splittype

  Title   : splittype
  Usage   : $splittype = $fuzzy->splittype();
  Function: get/set the split splittype
  Returns : the splittype of split feature (join, order)
  Args    : splittype to set

=cut

sub splittype {
    my($self) = @_;
    $self->throw_not_implemented();
}


=head2 is_single_sequence

  Title   : is_single_sequence
  Usage   : if($splitloc->is_single_sequence()) {
                print "Location object $splitloc is split ".
                      "but only across a single sequence\n";
	    }
  Function: Determine whether this location is split across a single or
            multiple sequences.
  Returns : TRUE if all sublocations lie on the same sequence as the root
            location (feature), and FALSE otherwise.
  Args    : none

=cut

sub is_single_sequence {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head1 Bio::LocationI methods

Bio::LocationI inherited methods follow

=head2 min_start

  Title   : min_start
  Usage   : my $minstart = $location->min_start();
  Function: Get minimum starting location of feature startpoint   
  Returns : integer or undef if no maximum starting point.
  Args    : none

=cut

=head2 max_start

  Title   : max_start
  Usage   : my $maxstart = $location->max_start();
  Function: Get maximum starting location of feature startpoint  
  Returns : integer or undef if no maximum starting point.
  Args    : none

=cut

=head2 start_pos_type

  Title   : start_pos_type
  Usage   : my $start_pos_type = $location->start_pos_type();
  Function: Get start position type (ie <,>, ^) 
  Returns : type of position coded as text 
            ('BEFORE', 'AFTER', 'EXACT','WITHIN', 'BETWEEN')
  Args    : none

=cut

=head2 min_end

  Title   : min_end
  Usage   : my $minend = $location->min_end();
  Function: Get minimum ending location of feature endpoint 
  Returns : integer or undef if no minimum ending point.
  Args    : none

=cut

=head2 max_end

  Title   : max_end
  Usage   : my $maxend = $location->max_end();
  Function: Get maximum ending location of feature endpoint 
  Returns : integer or undef if no maximum ending point.
  Args    : none

=cut

=head2 end_pos_type

  Title   : end_pos_type
  Usage   : my $end_pos_type = $location->end_pos_type();
  Function: Get end position type (ie <,>, ^) 
  Returns : type of position coded as text 
            ('BEFORE', 'AFTER', 'EXACT','WITHIN', 'BETWEEN')
  Args    : none

=cut

=head2 seq_id

  Title   : seq_id
  Usage   : my $seqid = $location->seq_id();
  Function: Get/Set seq_id that location refers to
  Returns : seq_id
  Args    : [optional] seq_id value to set

=cut

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

=cut

=head2 to_FTstring

  Title   : to_FTstring
  Usage   : my $locstr = $location->to_FTstring()
  Function: returns the FeatureTable string of this location
  Returns : string
  Args    : none

=cut

1;

