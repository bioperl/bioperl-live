# $Id$
#
# BioPerl module for Bio::SeqFeature::PositionProxy
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::PositionProxy - handle features when truncation/revcom sequences span a feature

=head1 SYNOPSIS

   $proxy = new Bio::SeqFeature::PositionProxy ( -loc => $loc,
                                                 -parent => $basefeature);

   $seq->add_SeqFeature($feat);



=head1 DESCRIPTION

PositionProxy is a Proxy Sequence Feature to handle truncation
and revcomp without duplicating all the data within the sequence features.
It holds a new location for a sequence feature and the original feature
it came from to provide the additional annotation information.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Ewan Birney

Ewan Birney E<lt>birney@sanger.ac.ukE<gt>

=head1 DEVELOPERS

This class has been written with an eye out of inheritence. The fields
the actual object hash (in addition to those defined in superclasses) are:

  _gsf_peer      = the seqfeature peer object
  _gsf_seq
  _gsf_location
  _gsf_alternative_locations

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::SeqFeature::PositionProxy;
use vars qw( @ISA $AUTOLOAD );
use strict;

use Bio::RelRange;
use Bio::SeqFeature::SegmentI;
use Bio::SeqFeatureI;
use Bio::AnnotatableI;
use Bio::AlternativeLocationHolderI;
@ISA = qw( Bio::RelRange
           Bio::SeqFeature::SegmentI
           Bio::SeqFeatureI
           Bio::AnnotatableI
           Bio::AlternativeLocationHolderI );

use Bio::DB::GFF::Util::Rearrange; # for 'rearrange'

sub new {
  my ( $caller, @args ) = @_;
  my $self = $caller->SUPER::new( @args );

  my ( $feature, $location );
  if( scalar( @args ) && $args[ 0 ] =~ /^-/ ) {
    ( $feature, $location ) =
      rearrange( [ [ qw( PARENT PEER ) ],
                   [ qw( LOCATION LOC ) ] ], @args );
  }

  if( !defined( $feature ) ||
      !ref( $feature ) ||
      !$feature->isa( 'Bio::SeqFeatureI' ) ) {
    $self->throw( "Must have a Bio::SeqFeatureI peer, not $feature, a ".ref( $feature ) );
  }

  if( $feature->isa( "Bio::SeqFeature::PositionProxy" ) ) {
    $feature = $feature->peer();
  }
  $self->peer( $feature );

  if( defined( $location ) &&
      ( !ref( $location ) || !$location->isa( 'Bio::LocationI' ) ) ) {
    $self->throw( "Must have a location, not a [$location]" );
  }
  $self->location( $location ) if defined( $location );

  return $self;
} # new(..)

=head2 peer

 Title   : peer
 Usage   : my $peer_feature = $proxy->peer()
 Function: Get/Set the L<Bio::SeqFeatureI> peer of this proxy
 Returns : The current (or former, if used as a set method) L<Bio::SeqFeatureI>
           peer of this PositionProxy.
 Args    : [optional] a new peer.

=cut

sub peer {
  my $self = shift;
  my $new_value = shift;

  my $old_value = $self->{ '_gsf_peer' };
  if( defined( $new_value ) ) {
    unless( ref( $new_value ) && $new_value->isa( 'Bio::SeqFeatureI' ) ) {
      $self->throw( "Unable to set the peer to $new_value.  It is a ".
                    ref( $new_value ).
                    ", not a Bio::SeqFeatureI." );
    }
    $self->{ '_gsf_peer' } = $new_value;
  }
  return $old_value;
} # peer(..)

=head2 location

 Title   : location
 Usage   : my $location = $seqfeature->location()
 Function: returns a location object suitable for identifying location
	   of feature on sequence or parent feature
 Returns : Bio::LocationI object
 Args    : [optional] Bio::LocationI object to set the value to.


=cut

sub location {
  my $self = shift;
  my $new_location = shift;

  my $old_location = $self->{ '_location' } ||
    Bio::Location::Simple->new(
      '-seq_id' => $self->seq_id(),
      '-start' => $self->start(),
      '-end' => $self->end(),
      '-strand' => $self->strand()
    );
  if( defined( $new_location ) ) {
    unless( ref( $new_location ) and $new_location->isa( 'Bio::LocationI' ) ) {
      $self->throw( "object $new_location pretends to be a location but ".
                    "does not implement Bio::LocationI" );
    }
    $self->{ '_location' } = $new_location;
  }
  return $old_location;
} # location(..)

sub add_alternative_locations {
    my $self = shift;
    foreach (@_) {
      $self->throw("object $_ pretends to be a location but ".
			 "does not implement Bio::LocationI")
          unless ref($_) and $_->isa('Bio::LocationI');
      push @{$self->{'_gsf_alternative_locations'}},$_;
    }
}

*add_alternative_location = \&add_alternative_locations;

=head2 alternative_locations

 Title   : alternative_locations
 Usage   : @locations = $seqfeature->alternative_locations([$seq_id])
 Function: returns alternative locations
 Returns : list of alternative locations
 Args    : optionally, a seq_id to filter on

=cut

sub alternative_locations {
    my $self = shift;
    my $seqid_filter = shift;
    return unless $self->{'_gsf_alternative_locations'};
    if ( $seqid_filter ) {
       return grep {$seqid_filter eq $_->seq_id} @{$self->{'_gsf_alternative_locations'}};
    } else {
       return @{$self->{'_gsf_alternative_locations'}};
    }
}

=head2 clear_alternative_locations

 Title   : clear_alternative_locations
 Usage   : $seqfeature->clear_alternative_locations([$seqid])
 Function: clears all alternative locations
 Returns : void
 Args    : optionally, a seq_id to clear locations on

=cut

sub clear_alternative_locations {
    my $self = shift;
    my $seqid_filter = shift;
    return unless $self->{'_gsf_alternative_locations'};
    if ( $seqid_filter ) {
       my @locations = grep {$seqid_filter ne $_->seq_id} @{$self->{'_gsf_alternative_locations'}};
       return $self->{'_gsf_alternative_locations'} = \@locations;
    } else {
       $self->{'_gsf_alternative_locations'} = [];
   }
}

=head2 attach_seq

 Title   : attach_seq
 Usage   : $sf->attach_seq($seq)
 Function: Attaches a Bio::Seq object to this feature. This
           Bio::Seq object is for the *entire* sequence: ie
           from 1 to 10000
 Example :
 Returns : TRUE on success
 Args    : a Bio::PrimarySeqI compliant object

=cut

sub attach_seq {
  my ( $self, $seq ) = @_;

  if( !( $seq && ref( $seq ) && $seq->isa( "Bio::PrimarySeqI" ) ) ) {
    $self->throw( "Must attach Bio::PrimarySeqI objects to SeqFeatures" );
  }

  $self->{ '_gsf_seq' } = $seq;

  # If we share the same abs_seq_id with our peer then we can attatch
  # it to the peer as well.
  if( ( $self->abs_seq_id() eq $self->peer()->abs_seq_id() ) &&
      not defined( $self->peer()->seq() ) ) {
    $self->peer()->attach_seq( $seq );
  }

  # Attatch to sub features if they want it
  foreach ( $self->sub_SeqFeature() ) {
    $_->attach_seq($seq);
  }

  return 1;
} # attach_seq(..)

=head2 seq

 Title   : seq
 Usage   : $tseq = $sf->seq()
 Function: returns the truncated sequence (if there) for this
 Example :
 Returns : sub seq (a Bio::PrimarySeqI compliant object) on attached sequence
           bounded by start & end, or undef if there is no sequence attached
 Args    : none


=cut

sub seq {
  my ( $self, $arg ) = @_;

  if( defined $arg ) {
    $self->throw("Calling SeqFeature::Generic->seq with an argument. You probably want attach_seq");
  }

  if ( !exists $self->{ '_gsf_seq' } ) {
    if( $self->abs_seq_id() eq $self->peer()->abs_seq_id() ) {
      my $peer_seq = $self->peer()->entire_seq();
      if( defined( $peer_seq ) ) {
        $self->attach_seq( $peer_seq );
      }
    } else {
      return undef;
    }
  }

  # assumming our seq object is sensible, it should not have to yank
  # the entire sequence out here.
  my $seq = $self->{'_gsf_seq'}->trunc( $self->start(), $self->end() );

  if ( $self->strand == -1 ) {
    # ok. this does not work well (?)
    #print STDERR "Before revcom", $seq->str, "\n";
    $seq = $seq->revcom;
    #print STDERR "After  revcom", $seq->str, "\n";
  }

  return $seq;
} # seq(..)

=head2 entire_seq

 Title   : entire_seq
 Usage   : $whole_seq = $sf->entire_seq()
 Function: gives the entire sequence that this seqfeature is attached to
 Example :
 Returns : a Bio::PrimarySeqI compliant object, or undef if there is no
           sequence attached
 Args    :


=cut

sub entire_seq {
  my ( $self ) = @_;

  if ( !exists $self->{ '_gsf_seq' } ) {
    if( $self->abs_seq_id() eq $self->peer()->abs_seq_id() ) {
      my $peer_seq = $self->peer()->seq();
      if( defined( $peer_seq ) ) {
        $self->attach_seq( $peer_seq );
      }
    } else {
      return undef;
    }
  }
  return $self->{'_gsf_seq'};
} # entire_seq(..)

=head2 orientation_policy

  Title   : orientation_policy
  Usage   : my $orientation_policy =
              $range->orientation_policy( [new_policy] );
  Function: Get/Set the oriention policy that this RelRangeI uses.
  Returns : The current (or former, if used as a set method) value of
            the orientation policy.
  Args    : [optional] A new (string) orientation_policy value.

  The BioPerl community has various opinions about how coordinates
  should be returned when the strand is negative.  Some folks like the
  start to be the lesser-valued position in all circumstances
  ('independent' of the strand value).  Others like the start to be
  the lesser-valued position when the strand is 0 or 1 and the
  greater-valued position when the strand is -1 ('dependent' on the
  strand value).  Others expect that the start and end values are
  whatever they were set to ('ignorant' of the strand value).

  Legal values of orientation_policy are:
      Value          Assertion
   ------------   -------------------
   'independent'  ( start() <= end() )
   'dependent'    (( strand() < 0 )?( end() <= start() ):( start() <= end() ))
   'ignorant'     # No assertion.

  See also ensure_orientation().  Note that changing the
  orientation_policy will not automatically ensure that the
  orientation policy assertion holds, so you should call
  ensure_orientation() also.

=cut

sub orientation_policy {
  my $self = shift;
  unless( $self->peer() ) {
    return 'ignorant';
  }
  $self->peer()->orientation_policy( @_ );
} # orientation_policy(..)

# Internal overridable getter/setter for the actual stored value of
# seq_id.  Delegates to the location object at $self->{ '_location' }
# if there is one; otherwise delegates to the superclass _seq_id(..).
sub _seq_id {
  my $self = shift;
  if( $self->{ '_location' } ) {
    my ( $new_val ) = @_;
    my $old_val = $self->{ '_location' }->seq_id();
    if( defined $new_val ) {
      $self->{ '_location' }->seq_id( $new_val );
    }
    return $old_val;
  } else {
    return $self->SUPER::_seq_id( @_ );
  }
} # _seq_id(..)

# Internal overridable getter/setter for the actual stored value of
# start.  Delegates to the location object at $self->{ '_location' }
# if there is one; otherwise delegates to the superclass _start(..).
sub _start {
  my $self = shift;
  if( $self->{ '_location' } ) {
    my ( $new_val ) = @_;
    my $old_val = $self->{ '_location' }->start();
    if( defined $new_val ) {
      $self->{ '_location' }->start( $new_val );
    }
    return $old_val;
  } else {
    return $self->SUPER::_start( @_ );
  }
} # _start(..)

# Internal overridable getter/setter for the actual stored value of
# end.  Delegates to the location object at $self->{ '_location' }
# if there is one; otherwise delegates to the superclass _end(..).
sub _end {
  my $self = shift;
  if( $self->{ '_location' } ) {
    my ( $new_val ) = @_;
    my $old_val = $self->{ '_location' }->end();
    if( defined $new_val ) {
      $self->{ '_location' }->end( $new_val );
    }
    return $old_val;
  } else {
    return $self->SUPER::_end( @_ );
  }
} # _end(..)

# Internal overridable getter/setter for the actual stored value of
# strand.  Delegates to the location object at $self->{ '_location' }
# if there is one; otherwise delegates to the superclass _strand(..).
sub _strand {
  my $self = shift;
  if( $self->{ '_location' } ) {
    my ( $new_val ) = @_;
    my $old_val = $self->{ '_location' }->strand();
    if( defined $new_val ) {
      $self->{ '_location' }->strand( $new_val );
    }
    return $old_val;
  } else {
    return $self->SUPER::_strand( @_ );
  }
} # _strand(..)

=head2 add_alternative_locations

 Title   : add_alternative_locations
 Usage   : $seqfeature->add_alternative_locations(@locationi)
 Function: adds new alternative location to the object
 Returns : void
 Args    : one or more LocationI-implementing object

 This adds one or more alternative locations to the feature.  These are
 to be viewed as alternative coordinate systems, such as
 assembly-to-assembly alignments, and not as alternative locations in
 the same coordinate space.

=cut

=head2 Autogenerated Methods

Any method that does not pertain to location/range/position will be delegated
to AUTOLOAD and treated as a call to the peer's method of the same name.
For instance, this call:

  @subfeatures = $feature->sub_SeqFeature();

is equivalent to this call:

  @subfeatures = $feature->peer()->sub_SeqFeature();

=cut

sub AUTOLOAD {
  my( $pack, $func_name ) = ( $AUTOLOAD =~ /(.+)::([^:]+)$/ );
  my $self = shift;

  # ignore DESTROY calls
  return if $func_name eq 'DESTROY';

  if( $self->{ '_gsf_peer' }->can( $func_name ) ) {
    return $self->{ '_gsf_peer' }->$func_name( @_ );
  } else {
    # error message of last resort
    $self->throw( "Can't locate object method \"$func_name\" via package \"$pack\"" );
  }
} # AUTOLOAD(..)

1;

__END__
