package Bio::Graphics::Feature;

=head1 NAME

Bio::Graphics::Feature - A simple feature object for use with Bio::Graphics::Panel

=head1 SYNOPSIS

 use Bio::Graphics::Feature;

 # create a simple feature with no internal structure
 $f = Bio::Graphics::Feature->new(-start => 1000,
                                  -stop  => 2000,
                                  -type  => 'transcript',
                                  -name  => 'alpha-1 antitrypsin',
				  -desc  => 'an enzyme inhibitor',
                                 );

 # create a feature composed of multiple segments, all of type "similarity"
 $f = Bio::Graphics::Feature->new(-segments => [[1000,1100],[1500,1550],[1800,2000]],
                                  -name     => 'ABC-3',
                                  -type     => 'gapped_alignment',
                                  -subtype  => 'similarity');

 # build up a gene exon by exon
 $e1 = Bio::Graphics::Feature->new(-start=>1,-stop=>100,-type=>'exon');
 $e2 = Bio::Graphics::Feature->new(-start=>150,-stop=>200,-type=>'exon');
 $e3 = Bio::Graphics::Feature->new(-start=>300,-stop=>500,-type=>'exon');
 $f  = Bio::Graphics::Feature->new(-segments=>[$e1,$e2,$e3],-type=>'gene');

=head1 DESCRIPTION

This is a simple Bio::SeqFeatureI-compliant object that is compatible
with Bio::Graphics::Panel.  With it you can create lightweight feature
objects for drawing.

All methods are as described in L<Bio::SeqFeatureI> with the following additions:

=head2 The new() Constructor

 $feature = Bio::Graphics::Feature->new(@args);

This method creates a new feature object.  You can create a simple
feature that contains no subfeatures, or a hierarchically nested object.

Arguments are as follows:

*  -start       the start position of the feature
*  -end         the stop position of the feature
*  -stop        an alias for end
*  -name        the feature name (returned by seqname())
*  -type        the feature type (returned by primary_tag())
*  -source      the source tag
  -desc        a description of the feature
*  -segments    a list of subfeatures (see below)
  -subtype     the type to use when creating subfeatures
*  -strand      the strand of the feature (one of -1, 0 or +1)
*  -id          an alias for -name (but now for -seq_id)
*  -seqname     an alias for -name
*  -primary_id  an alias for -name (but now for unique_id)
*  -display_id  an alias for -name (but now for display_name)
*  -display_name an alias for -name  (no longer an alias)
*  -attributes  a hashref of tag value attributes, in which the key is the tag
               and the value is an array reference of values

The subfeatures passed in -segments may be an array of
Bio::Graphics::Feature objects, or an array of [$start,$stop]
pairs. Each pair should be a two-element array reference.  In the
latter case, the feature type passed in -subtype will be used when
creating the subfeatures.

If no feature type is passed, then it defaults to "feature".

=head2 Non-SeqFeatureI methods

A number of new methods are provided for compatibility with
Ace::Sequence, which has a slightly different API from SeqFeatureI:

=over 4

=item add_segment(@segments)

Add one or more segments (a subfeature).  Segments can either be
Feature objects, or [start,stop] arrays, as in the -segments argument
to new().  The feature endpoints are automatically adjusted.

=item segments()

An alias for sub_SeqFeature().

=item merged_segments()

Another alias for sub_SeqFeature().

=item stop()

An alias for end().

=item name()

An alias for seqname().

=item exons()

An alias for sub_SeqFeature() (you don't want to know why!)

=back

=cut

use strict;
use Bio::SeqFeature::Generic;

use vars '$VERSION','@ISA';
$VERSION = '1.41';
@ISA  = qw( Bio::SeqFeature::Generic );

# usage:
# Bio::Graphics::Feature->new(
#                         -start => 1,
#                         -end   => 100,
#                         -name  => 'fred feature',
#                         -strand => +1);
#
# Alternatively, use -segments => [ [start,stop],[start,stop]...]
# to create a multisegmented feature.
sub new {
  my $caller = shift;
  my $self = $caller->SUPER::new( @_ ); 

  my ( $desc, $subtype, $url ) =
    $self->_rearrange( [
      'DESC',
      'SUBTYPE',
      'URL'
    ], @_ );

  $self->desc( $desc ) if defined( $desc );
  $self->subtype( $subtype ) if defined( $subtype );
  $self->url( $url ) if defined( $url );

  unless( $self->type() ) {
    $self->type( 'feature' );
  }

  return $self;
} # new(..)

=head2 subtype

 Title   : subtype
 Usage   : $feature->subtype( [$new_subtype] )
 Function: Getter/setter for the type of new subfeatures of this feature.
 Returns : the current (or former, if used as a set method) subtype (either a
           string or a Bio::SeqFeature::TypeI)
 Args    : (optional) A new subtype (either a string or, preferably, a
           Bio::SeqFeature::TypeI)
 Status  : Public

=cut

sub subtype {
  my $self = shift;
  my ( $new_value ) = @_;

  my $old_value = $self->{ '_subtype' };
  if( defined( $new_value ) ) {
    $self->{ '_subtype' } = $new_value;
  }
  return $old_value;
} # subtype(..)

=head2 _create_feature

 Title   : create_feature
 Usage   : my $new_seq_feature = $feature->_create_feature( $feature_data );
 Function: Factory method for instantiating a SeqFeatureI object.
 Returns : A new L<Bio::Graphics::Feature> object, or $feature_data.
 Args    : A single argument of any data type (see below).
 Status  : Protected

 The single argument may be of any type.  _create_feature(..) will be
 called by _insert_feature whenever the feature argument to that
 method is not a SeqFeatureI object.  If _create_feature(..) is unable
 to make a feature from the given data it must return the argument it
 was given.

=cut

## This one makes features when given [ $start, $end ] pairs.  New
## features will have the type returned by the $self->subtype() method
## or, if that is undef, the result of $self->type().  They will have
## $self as their $seq_id.  The given start and end positions are
## interpreted relative to $self.
sub _create_feature {
  my $self = shift;
  my $feature_data = shift;
  if( ref( $feature_data ) eq 'ARRAY' ) {
    # This is for feature data like [ [ start0, end0 ], [ start1, end1 ] ]
    my ( $start, $end ) = @{ $feature_data };
    
    # The following line should be unnecessary.
    #next unless defined $start && defined $end;
    
    # Strandedness defaults to our own, but start > end forces -1.
    my $strand = $self->strand();
    if( $start > $end ) {
      ( $start, $end ) = ( $end, $start );
      $strand = -1;
    }
    return $self->new( -start  => $start,
                       -end    => $end,
                       -strand => $strand,
                       -type   => $self->subtype() || $self->type(),
                       -seq_id => $self,
                       -parent => $self ) || $feature_data;
  } elsif( !ref( $feature_data ) ) {
    # This is for feature data like [ 'ref1:start1-end1', 'ref2:start2-end2' ].
    my ( $seq_id, $start, $end ) =
      ( $feature_data =~ /(.+):(-?\d+)(?:-|\.\.)(-?\d+)/ );

    ## TODO: REMOVE
    #print STDERR "_create_feature(..): $seq_id:$start-$end\n";
    
    # Strandedness defaults to our own, but start > end forces -1.
    my $strand = $self->strand();
    if( $start > $end ) {
      ( $start, $end ) = ( $end, $start );
      $strand = -1;
    }
    ## TODO: ERE I AM.  Need to make the coords relative to $self if possible.
    return $self->new( -start  => $start,
                       -end    => $end,
                       -strand => $strand,
                       -type   => $self->subtype() || $self->type(),
                       -seq_id => $seq_id,
                       -parent => $self ) || $feature_data;
  }
  # If we're at this point then we've been unable to make a new feature.
  return $feature_data;
} # _create_feature(..)

sub add_segment {
  shift->add_feature( @_ );
}

sub name {
  return shift->display_name( @_ );
}

sub info {
  shift->name( @_ );
}
sub exons {
  shift->features( @_ );
}
sub merged_segments {
  shift->features( @_ );
}
sub segments {
  shift->features( @_ );
}
sub method {
  shift->type( @_ );
}
sub source {
  shift->source_tag( @_ );
}
sub stop {
  shift->end( @_ );
}
sub target { return; }
sub hit    { return; }
sub ref {
  shift->seq_id( @_ );
}

sub phase {
  shift->frame( @_ );
}

sub dna {
  shift->seq( @_ );
}

=head2 accession_number

 Title   : accession_number
 Usage   : $unique_biological_key = $obj->accession_number;
 Function: Returns the unique biological id for a sequence, commonly
           called the accession_number. For sequences from established
           databases, the implementors should try to use the correct
           accession number. Notice that primary_id() provides the
           unique id for the implemetation, allowing multiple objects
           to have the same accession number in a particular implementation.

           For sequences with no accession number, this method should return
           "unknown".
 Returns : A string
 Args    : None


=cut

sub accession_number {
  return 'unknown';
}

=head2 alphabet

 Title   : alphabet
 Usage   : if( $obj->alphabet eq 'dna' ) { /Do Something/ }
 Function: Returns the type of sequence being one of
           'dna', 'rna' or 'protein'. This is case sensitive.

           This is not called <type> because this would cause
           upgrade problems from the 0.5 and earlier Seq objects.

 Returns : a string either 'dna','rna','protein'. NB - the object must
           make a call of the type - if there is no type specified it
           has to guess.
 Args    : none
 Status  : Virtual


=cut

sub alphabet{
  return 'dna'; # no way this will be anything other than dna!
}

=head2 desc

 Title   : desc
 Usage   : $seqobj->desc($string) or $seqobj->desc()
 Function: Sets or gets the description of the sequence
 Example :
 Returns : The description
 Args    : The description or none


=cut

sub desc {
  my $self = shift;
  my $d    = $self->{ '_desc' };
  $self->{ '_desc' } = shift if @_;
  $d;
}

## TODO: ?
#=head2 location
#
# Title   : location
# Usage   : my $location = $seqfeature->location()
# Function: returns a location object suitable for identifying location
#	   of feature on sequence or parent feature
# Returns : Bio::LocationI object
# Args    : none
#
#=cut
#
#sub location {
#   my $self = shift;
#   require Bio::Location::Split unless Bio::Location::Split->can('new');
#   my $location;
#   if (my @segments = $self->segments) {
#       $location = Bio::Location::Split->new();
#       foreach (@segments) {
#	 $location->add_sub_Location($_);
#       }
#   } else {
#       $location = $self;
#   }
#   $location;
#}
#
#sub coordinate_policy {
#   require Bio::Location::WidestCoordPolicy unless Bio::Location::WidestCoordPolicy->can('new');
#   return Bio::Location::WidestCoordPolicy->new();
#}

sub min_start { shift->low }
sub max_start { shift->low }
sub min_end   { shift->high }
sub max_end   { shift->high}
sub start_pos_type { 'EXACT' }
sub end_pos_type   { 'EXACT' }
sub to_FTstring {
  my $self = shift;
  my $low  = $self->min_start;
  my $high = $self->max_end;
  return "$low..$high";
}
sub db { return }


# This probably should be deleted.  Not sure why it's here, but might
# have been added for Ace::Sequence::Feature-compliance.
sub introns { return; }

# get/set the configurator (Bio::Graphics::FeatureFile) for this feature
sub configurator {
  my $self = shift;
  my $d = $self->{ '_configurator' };
  $self->{ '_configurator' } = shift if @_;
  $d;
}

# get/set the url for this feature
sub url {
  my $self = shift;
  my $d = $self->{ '_url' };
  $self->{ '_url' } = shift if @_;
  $d;
}

# make a link
sub make_link {
  my $self = shift;
  if (my $url = $self->url) {
    return $url;
  }

  elsif (my $configurator = $self->configurator) {
    return $configurator->make_link($self);
  }

  else {
    return;
  }
} # make_link(..)

## TODO: REMOVE?  Why is this here?
sub DESTROY { }

1;

__END__

=head1 SEE ALSO

L<Bio::Graphics::Panel>,L<Bio::Graphics::Glyph>,
L<GD>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
