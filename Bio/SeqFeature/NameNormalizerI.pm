package Bio::SeqFeature::NameNormalizerI;

# $Id$

=head1 NAME

Bio::SeqFeature::NameNormalizerI - A name normalizer for
L<Bio::SeqFeatureI> objects.  There are so many ways to name a feature,
and so many ontologies, etc.  L<NameNormalizerI>s can take a SeqFeatureI
object and set the display_name to a standard, normalized name for
that feature.

=head1 SYNOPSIS

  my $normalizer =
    Bio::DB::LocusLinkHugoNormalizer->new( @args );
  my $iterator = $seq_feature_collection->features( '-iterator' => 1 );
  while( $iterator->has_more_elements() ) {
    my $seq_feature = $iterator->next_feature();
    $normalizer->normalize( $seq_feature );
    print "The feature's normalized name is ".$seq_feature->display_name();
  }

=head1 DESCRIPTION

  A name normalizer for L<Bio::SeqFeatureI> objects.  There are so
  many ways to name a feature, and so many ontologies, etc.
  NameNormalizers can take a SeqFeatureI object and set the
  display_name to a standard, normalized name for that feature.  Of
  course, to do this, the SeqFeatureI object must already have its
  unique_id() or display_name() set to some value known by the
  normalizer to be an alias for some standard name.

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

use strict;
use vars qw( $VERSION @ISA );

$VERSION = '0.01';
use Bio::Root::RootI;
@ISA = qw( Bio::Root::RootI );

=head2 normalize

 Title   : normalize
 Usage   : $normalizer->normalize( $seq_feature );
 Function: Changes the given feature\'s display_name to something
           standard if a standard name can be found that corresponds
           to the feature\'s unique_id() or display_name() or string value,
           tested in that order.
 Returns : the given L<Bio::SeqFeatureI>
 Args    : a L<Bio::SeqFeatureI> object
 Status  : Public

=cut

sub normalize {
  shift->throw_not_implemented( @_ );
} # normalize()

1;

__END__
