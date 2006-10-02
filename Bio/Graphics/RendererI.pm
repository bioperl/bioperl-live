# $Id$

=head1 NAME

Bio::Graphics::RendererI - A renderer for the Bio::Graphics class that
renders Bio::SeqFeature::CollectionI objects onto
Bio::Graphics::Panels using configuration information provided by a
Bio::Graphics::ConfiguratorI.

=head1 SYNOPSIS

 # Get a renderer somehow, called $renderer

 # create a new panel and render contents a feature collection onto it
 my $config = new ConfigIO( $config_file )->getConfig();
 my $features = $data_provider->getCollection();
 my ( $tracks_rendered, $panel ) = $renderer->render( $features, $prefs );

=head1 DESCRIPTION

Renderer of Bio::SeqFeature::CollectionIs (collections of features)
onto a Bio::Graphics::Panel using a Bio::Graphics::ConfiguratorI for
general and track-specific rendering options.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR

Paul Edlefsen E<lt>paul@systemsbiology.orgE<gt>.

Copyright (c) 2003 Institute for Systems Biology

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Graphics::RendererI;
use strict;

use base qw(Bio::Root::RootI);

=head2 render

 Title   : render
 Usage   : ( $rendered, $panel ) =
               $renderer->render( $collection, $configurator [, $panel ] );
 Function: Renders the SeqFeatures in the given collection onto a
           Bio::Graphics::Panel (if no panel is given, one will be
           created), using the given Bio::Graphics::ConfiguratorI for general
           and track-specific rendering options.
 Returns : In a scalar context returns the number of tracks rendered.
           In a list context, returns a two-element list containing
           the number of features rendered and the panel.
 Args    : A Bio::SeqFeature::CollectionI and a Bio::Graphics::ConfiguratorI
           and optionally a Bio::Graphics::Panel.
 Status  : Public

=cut

sub render {
  shift->throw_not_implemented();
}

1;
