package Bio::SeqFeature::SegmentI.pm

# $Id $
# A Bio::SeqFeature::CollectionI that is also a Bio::LocationI.

=head1 NAME

Bio::SeqFeature::SegmentI -- A segment of sequence that contains features.

=head1 SYNOPSIS

# get a Bio::SeqFeature::SegmentI somehow
# perhaps a Bio::SeqFeature::SimpleSegment

 use Bio::SeqFeature::SimpleSegment;
 my $segment =
   new Bio::SeqFeature::SimpleSegment(
     '-start' => 10,
     '-end' => 3000,
     '-strand' => 1,
     '-seq_id' => $seq_id
   );
 $segment->add_features( @featurelist );

 # Now the -location passed to features(..) will be relative to the
 # location [10,3000] on the forward strand.
 $segment->features(
    -location => new Bio::Location::Simple( '-start' => 1, '-end' => 300 ),
    -rangetype =>'overlaps'
 );

=head1 DESCRIPTION

The Bio::SeqFeature::SegmentI interface represents a segment of a
sequence that can contain features.  As such it implements the
Bio::LocationI interface and the Bio::SeqFeature::CollectionI
interface.

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

use Bio::LocationI;
use Bio::SeqFeature::CollectionI;
@ISA = qw( Bio::LocationI
           Bio::SeqFeature::CollectionI );

use vars '$VERSION';
$VERSION = '1.00';

1;

__END__
