package Bio::DB::FeatureProviderI;

# $Id $
# A Bio::SeqFeature::CollectionProviderI with additional methods
# suiting an implementation with an external or non-simple backing
# store.

=head1 NAME

Bio::DB::FeatureProviderI -- A provider of collections of sequence
features from a database or other non-trivial backing store.

=head1 SYNOPSIS

 use Bio::DB::SimpleFeatureProvider;

 use Bio::SeqFeature::Generic;
 use Bio::SeqFeature::SimpleCollection;

 my $data_provider =
   Bio::SeqFeature::SimpleFeatureProvider->new();

 # Add some features

 $data_provider->insert_collection(
   new Bio::SeqFeature::SimpleCollection(
     new Bio::SeqFeature::Generic(
       -id => 'foo',
       -start => 10,
       -end => 100,
       -strand => -1,
       -primary => 'repeat',
       -source_tag => 'repeatmasker',
       -score  => 1000
     ),
     new Bio::SeqFeature::Generic(
       -id => 'bar',
       -start => 100,
       -end => 200,
       -strand => -1
     )
   );
 );

 # Add another feature
 my $baz =
   new Bio::SeqFeature::Generic(
     -id => 'baz',
     -start => 1,
     -end => 200
   );
 $data_provider->insert_collection(
   new Bio::SeqFeature::SimpleCollection( $baz )
 );

 # Update one that we'd previously inserted.

 $baz->strand( -1 );
 $data_provider->update_collection(
   new Bio::SeqFeature::SimpleCollection( $baz );
 );

=head1 DESCRIPTION

The Bio::DB::FeatureProviderI interface provides access to
Bio::SeqFeature::CollectionIs stored in a database or other (generally
external) backing store.  It is a Bio::SeqFeature::CollectionProviderI
with additional methods for managing the connection, caching, etc.

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
use vars qw( @ISA );

use Bio::SeqFeature::CollectionProviderI;
@ISA = qw( Bio::SeqFeature::CollectionProviderI );

use vars '$VERSION';
$VERSION = '1.00';

1;

__END__
