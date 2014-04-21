#
# BioPerl module for Bio::Location::Simple
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Location::Simple - Implementation of a Simple Location on a Sequence

=head1 SYNOPSIS

    use Bio::Location::Simple;

    my $location = Bio::Location::Simple->new(
        -start  => 1,
        -end    => 100,
        -strand => 1,
    );

    if( $location->strand == -1 ) {
        printf "complement(%d..%d)\n", $location->start, $location->end;
    } else {
        printf "%d..%d\n", $location->start, $location->end;
    }

=head1 DESCRIPTION

This is an implementation of Bio::LocationI to manage exact location
information on a Sequence: '22' or '12..15' or '16^17'.

You can test the type of the location using length() function () or
directly location_type() which can one of two values: 'EXACT' or
'IN-BETWEEN'.


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

=head1 AUTHOR - Heikki Lehvaslaiho

Email heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::Location::Simple;
use strict;

use base qw(Bio::Location::Atomic);

our %RANGEENCODE  = ('\.\.' => 'EXACT',
                     '\^'   => 'IN-BETWEEN' );

our %RANGEDECODE  = ('EXACT'      => '..',
                     'IN-BETWEEN' => '^' );

sub new { 
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($locationtype) = $self->_rearrange([qw(LOCATION_TYPE)],@args);

    $locationtype && $self->location_type($locationtype);

    return $self;
}

=head2 start

  Title   : start
  Usage   : $start = $loc->start();
  Function: get/set the start of this range
  Returns : the start of this range
  Args    : optionaly allows the start to be set
            using $loc->start($start)

=cut

sub start {
    my ($self, $value) = @_;
    $self->{'_start'} = $value if defined $value ;
  
    $self->throw("Only adjacent residues when location type ".
             "is IN-BETWEEN. Not [". $self->{'_start'}. "] and [".
             $self->{'_end'}. "]" )
        if defined $self->{'_start'} && defined $self->{'_end'} && 
        $self->location_type eq 'IN-BETWEEN' &&
        ($self->{'_end'} - 1 != $self->{'_start'});
    return $self->{'_start'};
}


=head2 end

  Title   : end
  Usage   : $end = $loc->end();
  Function: get/set the end of this range
  Returns : the end of this range
  Args    : optionaly allows the end to be set
          : using $loc->end($start)
  Note    : If start is set but end is undefined, this now assumes that start
            is the same as end but throws a warning (i.e. it assumes this is
            a possible error). If start is undefined, this now throws an
            exception.

=cut

sub end {
    my ($self, $value) = @_;
  
    $self->{'_end'} = $value if defined $value ;
    
    # Assume end is the same as start if not defined
    if (!defined $self->{'_end'}) {
        if (!defined $self->{'_start'}) {
            $self->warn('Can not set Bio::Location::Simple::end() equal to start; start not set');
            return;
        }
        $self->warn('Setting end to equal start['. $self->{'_start'}. ']');
        $self->{'_end'} = $self->{'_start'};
    }
    $self->throw("Only adjacent residues when location type ".
            "is IN-BETWEEN. Not [". $self->{'_start'}. "] and [".
             $self->{'_end'}. "]" )
        if defined $self->{'_start'} && defined $self->{'_end'} && 
        $self->location_type eq 'IN-BETWEEN' &&
        ($self->{'_end'} - 1 != $self->{'_start'});
  
    return $self->{'_end'};
}

=head2 strand

  Title   : strand
  Usage   : $strand = $loc->strand();
  Function: get/set the strand of this range
  Returns : the strandedness (-1, 0, +1)
  Args    : optionaly allows the strand to be set
          : using $loc->strand($strand)

=cut

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
    if ($self->location_type eq 'IN-BETWEEN' ) {
        return 0;
    } else {
        return abs($self->end - $self->start) + 1;
    }
}


=head2 min_start

  Title   : min_start
  Usage   : my $minstart = $location->min_start();
  Function: Get minimum starting location of feature startpoint
  Returns : integer or undef if no minimum starting point.
  Args    : none

=cut

=head2 max_start

  Title   : max_start
  Usage   : my $maxstart = $location->max_start();
  Function: Get maximum starting location of feature startpoint.

            In this implementation this is exactly the same as min_start().

  Returns : integer or undef if no maximum starting point.
  Args    : none

=cut

=head2 start_pos_type

  Title   : start_pos_type
  Usage   : my $start_pos_type = $location->start_pos_type();
  Function: Get start position type (ie <,>, ^).

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

            In this implementation this is exactly the same as min_end().

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

=head2 location_type

  Title   : location_type
  Usage   : my $location_type = $location->location_type();
  Function: Get location type encoded as text
  Returns : string ('EXACT' or 'IN-BETWEEN')
  Args    : 'EXACT' or '..' or 'IN-BETWEEN' or '^'

=cut

sub location_type {
    my ($self, $value) = @_;

    if( defined $value || ! defined $self->{'_location_type'} ) {
        $value = 'EXACT' unless defined $value;
        $value = uc $value;
        if (! defined $RANGEDECODE{$value}) {
            $value = '\^' if $value eq '^';
            $value = '\.\.' if $value eq '..';
            $value = $RANGEENCODE{$value};
        }
        $self->throw("Did not specify a valid location type. [$value] is no good")
            unless defined $value;
        $self->{'_location_type'} = $value;
    }
    $self->throw("Only adjacent residues when location type ".
            "is IN-BETWEEN. Not [". $self->{'_start'}. "] and [".
            $self->{'_end'}. "]" )
        if $self->{'_location_type'} eq 'IN-BETWEEN' &&
        defined $self->{'_start'} &&
        defined $self->{'_end'} &&
            ($self->{'_end'} - 1 != $self->{'_start'});

    return $self->{'_location_type'};
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

=head2 to_FTstring

  Title   : to_FTstring
  Usage   : my $locstr = $location->to_FTstring()
  Function: returns the FeatureTable string of this location
  Returns : string
  Args    : none

=cut

sub to_FTstring { 
    my($self) = @_;

    my $str;
    if( $self->start == $self->end ) {
        $str =  $self->start;
    } else {
        $str = $self->start . $RANGEDECODE{$self->location_type} . $self->end;
    }
    if($self->is_remote() && $self->seq_id()) {
        $str = $self->seq_id() . ":" . $str;
    }
    if( defined $self->strand &&
        $self->strand == -1 ) {
        $str = "complement(".$str.")";
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

# comments, not function added by jason 
#
# trunc is untested, and as of now unannounced method for truncating a
# location.  This is to eventually be part of the procedure to
# truncate a sequence with annotation and properly remap the location
# of all the features contained within the truncated segment.

# presumably this might do things a little differently for the case 
# where the truncation splits the location in half
# 
# in short- you probably don't want to use  this method.

sub trunc {
    my ($self,$start,$end,$relative_ori) = @_;
    my $newstart  = $self->start - $start+1;
    my $newend    = $self->end   - $start+1;
    my $newstrand = $relative_ori * $self->strand;

    my $out;
    if( $newstart < 1 || $newend > ($end-$start+1) ) {
        $out = Bio::Location::Simple->new();
        $out->start($self->start);
        $out->end($self->end);
        $out->strand($self->strand);
        $out->seq_id($self->seqid);
        $out->is_remote(1);
    } else {
        $out = Bio::Location::Simple->new();
        $out->start($newstart);
        $out->end($newend);
        $out->strand($newstrand);
       $out->seq_id();
    }

    return $out;
}

1;

