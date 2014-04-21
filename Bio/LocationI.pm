#
# BioPerl module for Bio::LocationI
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::LocationI - Abstract interface of a Location on a Sequence

=head1 SYNOPSIS

    # get a LocationI somehow
    printf( "start = %d, end = %d, strand = %s, seq_id = %s\n", 
	    $location->start, $location->end, $location->strand,
	    $location->seq_id);
    print "location str is ", $location->to_FTstring(), "\n"; 


=head1 DESCRIPTION

This Interface defines the methods for a Bio::LocationI, an object
which encapsulates a location on a biological sequence.  Locations
need not be attached to actual sequences as they are stand alone
objects.  LocationI objects are used by L<Bio::SeqFeatureI> objects to
manage and represent locations for a Sequence Feature.

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

package Bio::LocationI;
use strict;

use Carp;

use base qw(Bio::RangeI);

=head2 location_type

  Title   : location_type
  Usage   : my $location_type = $location->location_type();
  Function: Get location type encoded as text
  Returns : string ('EXACT', 'WITHIN', 'IN-BETWEEN')
  Args    : none

=cut

sub location_type { 
    my ($self,@args) = @_;
    $self->throw_not_implemented();
}

=head2 start

  Title   : start
  Usage   : $start = $location->start();
  Function: Get the start coordinate of this location as defined by
            the currently active coordinate computation policy. In
            simple cases, this will return the same number as
            min_start() and max_start(), in more ambiguous cases like
            fuzzy locations the number may be equal to one or neither
            of both.

            We override this here from RangeI in order to delegate
            'get' to a L<Bio::Location::CoordinatePolicy> implementing
            object.  Implementing classes may also wish to provide
            'set' functionality, in which case they *must* override
            this method. The implementation provided here will throw
            an exception if called with arguments.

  Returns : A positive integer value.
  Args    : none

See L<Bio::Location::CoordinatePolicy> for more information

=cut

sub start {
    my ($self,@args) = @_;

    # throw if @args means that we don't support updating information
    # in the interface but will delegate to the coordinate policy object
    # for interpreting the 'start' value

    $self->throw_not_implemented if @args;
    return $self->coordinate_policy()->start($self);
}

=head2 end

  Title   : end
  Usage   : $end = $location->end();
  Function: Get the end coordinate of this location as defined by the
            currently active coordinate computation policy. In simple
            cases, this will return the same number as min_end() and
            max_end(), in more ambiguous cases like fuzzy locations
            the number may be equal to one or neither of both.

            We override this here from Bio::RangeI in order to delegate
            'get' to a L<Bio::Location::CoordinatePolicy> implementing
            object. Implementing classes may also wish to provide
            'set' functionality, in which case they *must* override
            this method. The implementation provided here will throw
            an exception if called with arguments.

  Returns : A positive integer value.
  Args    : none

See L<Bio::Location::CoordinatePolicy> and L<Bio::RangeI> for more
information

=cut

sub end {
    my ($self,@args) = @_;

    # throw if @args means that we don't support updating information
    # in the interface but will delegate to the coordinate policy object
    # for interpreting the 'end' value
    $self->throw_not_implemented if @args;
    return $self->coordinate_policy()->end($self);
}

=head2 min_start

  Title   : min_start
  Usage   : my $minstart = $location->min_start();
  Function: Get minimum starting point of feature.

            Note that an implementation must not call start() in this method.

  Returns : integer or undef if no minimum starting point.
  Args    : none

=cut

sub min_start {
    my($self) = @_;
    $self->throw_not_implemented();
}

=head2 max_start

  Title   : max_start
  Usage   : my $maxstart = $location->max_start();
  Function: Get maximum starting point of feature.

            Note that an implementation must not call start() in this method
            unless start() is overridden such as not to delegate to the
            coordinate computation policy object.

  Returns : integer or undef if no maximum starting point.
  Args    : none

=cut

sub max_start {
    my($self) = @_;
    $self->throw_not_implemented();
}

=head2 start_pos_type

  Title   : start_pos_type
  Usage   : my $start_pos_type = $location->start_pos_type();
  Function: Get start position type encoded as text

            Known valid values are 'BEFORE' (<5..100), 'AFTER' (>5..100), 
            'EXACT' (5..100), 'WITHIN' ((5.10)..100), 'BETWEEN', (5^6), with
            their meaning best explained by their GenBank/EMBL location string
            encoding in brackets.

  Returns : string ('BEFORE', 'AFTER', 'EXACT','WITHIN', 'BETWEEN')
  Args    : none

=cut

sub start_pos_type {
    my($self) = @_;
    $self->throw_not_implemented();
}


=head2 flip_strand

  Title   : flip_strand
  Usage   : $location->flip_strand();
  Function: Flip-flop a strand to the opposite
  Returns : None
  Args    : None

=cut


sub flip_strand {
    my $self= shift;
    $self->strand($self->strand * -1);
}

=head2 min_end

  Title   : min_end
  Usage   : my $minend = $location->min_end();
  Function: Get minimum ending point of feature. 

            Note that an implementation must not call end() in this method
            unless end() is overridden such as not to delegate to the
            coordinate computation policy object.

  Returns : integer or undef if no minimum ending point.
  Args    : none

=cut

sub min_end {
    my($self) = @_;
    $self->throw_not_implemented();
}

=head2 max_end

  Title   : max_end
  Usage   : my $maxend = $location->max_end();
  Function: Get maximum ending point of feature.

            Note that an implementation must not call end() in this method
            unless end() is overridden such as not to delegate to the
            coordinate computation policy object.

  Returns : integer or undef if no maximum ending point.
  Args    : none

=cut

sub max_end {
    my($self) = @_;
    $self->throw_not_implemented();
}

=head2 end_pos_type

  Title   : end_pos_type
  Usage   : my $end_pos_type = $location->end_pos_type();
  Function: Get end position encoded as text.

            Known valid values are 'BEFORE' (5..<100), 'AFTER' (5..>100), 
            'EXACT' (5..100), 'WITHIN' (5..(90.100)), 'BETWEEN', (5^6), with
            their meaning best explained by their GenBank/EMBL location string
            encoding in brackets.

  Returns : string ('BEFORE', 'AFTER', 'EXACT','WITHIN', 'BETWEEN')
  Args    : none

=cut

sub end_pos_type {
    my($self) = @_;
    $self->throw_not_implemented();
}

=head2 seq_id

  Title   : seq_id
  Usage   : my $seqid = $location->seq_id();
  Function: Get/Set seq_id that location refers to
  Returns : seq_id (a string)
  Args    : [optional] seq_id value to set

=cut

sub seq_id {
    my($self) = @_;
    $self->throw_not_implemented();
}

=head2 is_remote

 Title   : is_remote
 Usage   : $is_remote_loc = $loc->is_remote()
 Function: Whether or not a location is a remote location.

           A location is said to be remote if it is on a different
           'object' than the object which 'has' this
           location. Typically, features on a sequence will sometimes
           have a remote location, which means that the location of
           the feature is on a different sequence than the one that is
           attached to the feature. In such a case, $loc->seq_id will
           be different from $feat->seq_id (usually they will be the
           same).

           While this may sound weird, it reflects the location of the
           kind of AL445212.9:83662..166657 which can be found in GenBank/EMBL
           feature tables.

 Example : 
 Returns : TRUE if the location is a remote location, and FALSE otherwise
 Args    : Value to set to


=cut

sub is_remote{
    shift->throw_not_implemented();
}

=head2 coordinate_policy

  Title   : coordinate_policy
  Usage   : $policy = $location->coordinate_policy();
            $location->coordinate_policy($mypolicy); # set may not be possible
  Function: Get the coordinate computing policy employed by this object.

            See L<Bio::Location::CoordinatePolicyI> for documentation
            about the policy object and its use.

            The interface *does not* require implementing classes to
            accept setting of a different policy. The implementation
            provided here does, however, allow to do so.

            Implementors of this interface are expected to initialize
            every new instance with a
            L<Bio::Location::CoordinatePolicyI> object. The
            implementation provided here will return a default policy
            object if none has been set yet. To change this default
            policy object call this method as a class method with an
            appropriate argument. Note that in this case only
            subsequently created Location objects will be affected.

  Returns : A L<Bio::Location::CoordinatePolicyI> implementing object.
  Args    : On set, a L<Bio::Location::CoordinatePolicyI> implementing object.

See L<Bio::Location::CoordinatePolicyI> for more information


=cut

sub coordinate_policy {
    shift->throw_not_implemented();
}

=head2 to_FTstring

  Title   : to_FTstring
  Usage   : my $locstr = $location->to_FTstring()
  Function: returns the FeatureTable string of this location
  Returns : string
  Args    : none

=cut

sub to_FTstring { 
    my($self) = @_;
    $self->throw_not_implemented();
}

=head2 each_Location

 Title   : each_Location
 Usage   : @locations = $locObject->each_Location($order);
 Function: Conserved function call across Location:: modules - will
           return an array containing the component Location(s) in
           that object, regardless if the calling object is itself a
           single location or one containing sublocations.
 Returns : an array of Bio::LocationI implementing objects
 Args    : Optional sort order to be passed to sub_Location() for Splits

=cut

sub each_Location {
    my ($self,@args) = @_;
    $self->throw_not_implemented();
}


=head2 valid_Location

 Title   : valid_Location
 Usage   : if ($location->valid_location) {...};
 Function: boolean method to determine whether location is considered valid
           (has minimum requirements for a specific LocationI implementation)
 Returns : Boolean value: true if location is valid, false otherwise
 Args    : none

=cut

sub valid_Location {
    my ($self,@args) = @_;
    $self->throw_not_implemented();
}

1;

