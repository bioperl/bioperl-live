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

  -seq_id      the reference sequence
  -start       the start position of the feature
  -end         the stop position of the feature
  -stop        an alias for end
  -name        the feature name (returned by seqname())
  -type        the feature type (returned by primary_tag())
  -primary_tag the same as -type
  -source      the source tag
  -score       the feature score (for GFF compatibility)
  -desc        a description of the feature
  -segments    a list of subfeatures (see below)
  -subtype     the type to use when creating subfeatures
  -strand      the strand of the feature (one of -1, 0 or +1)
  -phase       the phase of the feature (0..2)
  -id          an alias for -name
  -seqname     an alias for -name
  -display_id  an alias for -name
  -display_name an alias for -name  (do you get the idea the API has changed?)
  -primary_id  unique database ID
  -url         a URL to link to when rendered with Bio::Graphics
  -configurator an object (like a Bio::Graphics::FeatureFile) that knows how 
               to configure the graphical representation of the object based
               on its type.
  -attributes  a hashref of tag value attributes, in which the key is the tag
               and the value is an array reference of values
  -factory     a reference to a feature factory, used for compatibility with
               more obscure parts of Bio::DB::GFF

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

=item url()

Get/set the URL that the graphical rendering of this feature will link to.

=item add_segment(@segments)

Add one or more segments (a subfeature).  Segments can either be
Feature objects, or [start,stop] arrays, as in the -segments argument
to new().  The feature endpoints are automatically adjusted.

=item segments()

An alias for sub_SeqFeature().

=item get_SeqFeatures()

Alias for sub_SeqFeature()

=item get_all_SeqFeatures()

Alias for sub_SeqFeature()

=item merged_segments()

Another alias for sub_SeqFeature().

=item stop()

An alias for end().

=item name()

An alias for seqname().

=item exons()

An alias for sub_SeqFeature() (you don't want to know why!)

=item configurator()

Get/set the configurator that knows how to adjust the graphical
representation of this feature based on its type.  Currently the only
configurator that will work is Bio::Graphics::FeatureFile.

=back

=cut

use strict;
use base 'Bio::Graphics::FeatureBase';

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
  my $self = shift->SUPER::new(@_);

  my %arg = @_;
  for my $option (qw(factory configurator)) {
    $self->{$option} = $arg{"-$option"} if exists $arg{"-$option"};
  }
  $self;
}

=head2 factory

 Title   : factory
 Usage   : $factory = $obj->factory([$new_factory])
 Function: Returns the feature factory from which this feature was generated.
           Mostly for compatibility with weird dependencies in gbrowse.
 Returns : A feature factory
 Args    : None

=cut

sub factory {
  my $self = shift;
  my $d = $self->{factory};
  $self->{factory} = shift if @_;
  $d;
}

=head2 display_name

 Title   : display_name
 Usage   : $id = $obj->display_name or $obj->display_name($newid);
 Function: Gets or sets the display id, also known as the common name of
           the Seq object.

           The semantics of this is that it is the most likely string
           to be used as an identifier of the sequence, and likely to
           have "human" readability.  The id is equivalent to the LOCUS
           field of the GenBank/EMBL databanks and the ID field of the
           Swissprot/sptrembl database. In fasta format, the >(\S+) is
           presumed to be the id, though some people overload the id
           to embed other information. Bioperl does not use any
           embedded information in the ID field, and people are
           encouraged to use other mechanisms (accession field for
           example, or extending the sequence object) to solve this.

           Notice that $seq->id() maps to this function, mainly for
           legacy/convenience issues.
 Returns : A string
 Args    : None or a new id


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

=head2 desc

 Title   : desc
 Usage   : $seqobj->desc($string) or $seqobj->desc()
 Function: Sets or gets the description of the sequence
 Example :
 Returns : The description
 Args    : The description or none

=head2 location

 Title   : location
 Usage   : my $location = $seqfeature->location()
 Function: returns a location object suitable for identifying location
	   of feature on sequence or parent feature
 Returns : Bio::LocationI object
 Args    : none


=head2 location_string

 Title   : location_string
 Usage   : my $string = $seqfeature->location_string()
 Function: Returns a location string in a format recognized by gbrowse
 Returns : a string
 Args    : none

This is a convenience function used by the generic genome browser. It
returns the location of the feature and its subfeatures in the compact
form "start1..end1,start2..end2,...".  Use
$seqfeature-E<gt>location()-E<gt>toFTString() to obtain a standard
GenBank/EMBL location representation.

=head2 configurator

 Title   : configurator
 Usage   : my $configurator = $seqfeature->configurator([$new_configurator])
 Function: Get/set an object that provides configuration information for this feature
 Returns : configurator object
 Args    : new configurator object (optional)

A configurator object provides hints to the Bio::Graphics::Feature as
to how to display itself on a canvas. Currently this stores the
Bio::Graphics::FeatureFile and descendents.

=cut


# get/set the configurator (Bio::Graphics::FeatureFile) for this feature
sub configurator {
  my $self = shift;
  my $d = $self->{configurator};
  $self->{configurator} = shift if @_;
  $d;
}

=head2 url

 Title   : url
 Usage   : my $url = $seqfeature->url([$new_url])
 Function: Get/set the URL associated with this feature
 Returns : a URL string
 Args    : new URL (optional)

Features link to URLs when displayed as a clickable image map. This
field holds that information.

=cut


# get/set the url for this feature
sub url {
  my $self = shift;
  my $d = $self->{url};
  $self->{url} = shift if @_;
  $d;
}

=head2 make_link

 Title   : make_link
 Usage   : my $url = $seqfeature->make_link()
 Function: Create a URL for the feature
 Returns : a URL string
 Args    : none

This method will invoke the configurator in order to turn the feature
into a link. Used by Bio::Graphics::Panel to create imagemaps.

=cut

# make a link
sub make_link {
  my $self = shift;

  if (my $url = $self->url) {
    return $url;
  }

  elsif (my $configurator = $self->configurator) {
    return $configurator->make_link($self) if $configurator->can('make_link');
  }

  else {
    return;
  }
}

1;

__END__

=head1 SEE ALSO

L<Bio::Graphics::Panel>,L<Bio::Graphics::Glyph>, L<Bio::Graphics::FeatureBase>
L<GD>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
