# $Id$

=head1 NAME

Bio::Graphics::AutoAggregatorFactoryI - A factory that creates
L<Bio::SeqFeature::AggregatorI> objects automatically when given a
L<Bio::Graphics::ConfiguratorI> and a section.

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

package Bio::Graphics::AutoAggregatorFactoryI;

use strict;
use vars qw( $VERSION @ISA );
use overload 
  '""' => 'toString';

use Bio::Root::RootI;

$VERSION = '0.01';
@ISA = qw( Bio::Root::RootI );

=head2 create_aggregator

 Title   : create_aggregator
 Usage   : my $aggregator = $factory->create_aggregator( $config, $section )
 Function: Creates and returns a new L<Bio::SeqFeature::AggregatorI> to
           aggregate features for the given section of the given Configurator.
 Returns : a new L<Bio::SeqFeature::AggregatorI> object.
 Args    : a L<Bio::Graphics::ConfiguratorI> and a string section label.
 Status  : Public

=cut

sub create_aggregator {
  shift->throw_not_implemented( @_ );
}

sub toString {
  shift->throw_not_implemented( @_ );
}

1;
