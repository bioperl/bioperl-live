# $Id$

=head1 NAME

Bio::Graphics::DasAutoAggregatorFactory - A factory that creates
L<Bio::SeqFeature::AggregatorI> objects automatically when given a
L<Bio::Graphics::ConfiguratorI> and a section, using the keys
DAS_id, DAS_category, and DAS_method from that section.

=head1 SYNOPSIS

# TODO

=head1 DESCRIPTION

# TODO

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR

Paul Edlefsen E<lt>paul@systemsbiology.orgE<gt>.

Copyright (c) 2003 Institute for Systems Biology

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=head1 CONTRIBUTORS

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Graphics::DasAutoAggregatorFactory;

use strict;
use vars qw( $VERSION @ISA );

use Bio::Graphics::AutoAggregatorFactoryI;
use Bio::Root::Root;

$VERSION = '0.01';
@ISA = qw( Bio::Root::Root Bio::Graphics::AutoAggregatorFactoryI );

use Bio::DB::GFF::Aggregator;
use Digest::MD5 'md5_hex';

sub new {
  my $pack = shift;
  return $pack->SUPER::new( @_ );
}

=head2 create_aggregator

 Title   : create_aggregator
 Usage   : my $aggregator = $factory->create_aggregator( $config, $section )
 Function: Creates and returns a new L<Bio::SeqFeature::AggregatorI> to
           aggregate features for the given section of the given Configurator.
 Returns : a new L<Bio::SeqFeature::AggregatorI> object.
 Args    : a L<Bio::Graphics::ConfiguratorI> and a string section label.
 Status  : Public

  Note that this also changes the feature name in the config file to match the aggregator\'s aggregated feature type.

=cut

sub create_aggregator {
  my $self = shift;
  my ( $config, $section ) = @_;
  my $DAS_category = $config->get( $section, 'das_category' );
  unless( $DAS_category ) {
    return;
  }
  ## TODO: It seems that the DAS_method just screws it up.
  my $DAS_method = #$config->get( $section, 'das_method' ) ||
    $config->get( $section, 'das_id' );
  my $meta_method = 'meta_'.( $DAS_method || md5_hex( rand() ) );

  ## TODO: REMOVE
  #warn "Creating aggregator with method $meta_method, main_method $DAS_category, and parts [ $DAS_method ].";
  my $aggregator = Bio::DB::GFF::Aggregator->new( '-method' => $meta_method,
                                             '-main_method' => $meta_method,
                                             '-sub_parts' => [ $DAS_method ]
                                           );

  $config->set( $section, 'feature', $meta_method.':'.$DAS_category );

  return $aggregator;
}

sub toString {
  return "DasAutoAggregatorFactory";
}

1;
