# $Id$
#
# BioPerl module for Bio::GFFI
#
# Cared for by Lincoln Stein <lstein@cshl.org>
#
# Copyright Lincoln Stein
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::GFFI - GFF-style access to a feature database

=head1 SYNOPSIS

  # Open up a feature database somehow...
  $db = Bio::DB::GFFI->new(@args);

  $segment = $db->segment(-name  => 'NT_29921.4',
                          -start => 1,
		          -end   => 1000000);

  # segment is a Bio::DB::GFF::Segment - compliant object

  # fetch a list of features
  @features = $segment->features(-type=>['type1','type2','type3']);

  # invoke a callback over features
  $segment->features(-type=>['type1','type2','type3'],
                     -callback => sub { ... }
		    );

  $stream   = $segment->get_seq_stream(-type=>['type1','type2','type3']);
  while (my $feature = $stream->next_seq) {
     # each feature is a Bio::SeqFeatureI-compliant object
  }

  # get all feature types
  @types   = $db->types;

  # count types
  %types   = $db->types(-enumerate=>1);

  @feature = $segment->get_feature_by_name($class=>$name);
  @feature = $segment->get_feature_by_target($target_name);
  @feature = $segment->get_feature_by_attribute($att1=>$value1,$att2=>$value2);
  $feature = $segment->get_feature_by_id($id);

  $error = $db->error;

=head1 DESCRIPTION

Bio::DB::GFFI is a simplified alternative interface to sequence annotation
databases used by the GFF-format classes (see
L<Bio::DB::GFF>). In this scheme, the genome is represented as a series of
features, a subset of which are named.  Named features can be used as
reference points for retrieving "segments" (see L<Bio::DB::GFF::SegmentI>),
and these can, in turn, be used as the basis for exploring the genome
further.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bio.perl.org

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Lincoln Stein

Email lstein@cshl.org

=head1 CONTRIBUTORS

Paul Edlefsen E<lt>paul@systemsbiology.orgE<gt>.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...

package Bio::DB::GFFI;
use strict;

use vars qw( @ISA );
use Bio::DB::SegmentProviderI;
use Bio::DB::AggregatingFeatureProviderI;
@ISA = qw( Bio::DB::SegmentProviderI Bio::DB::AggregatingFeatureProviderI );

=head2 search_notes

 Title   : search_notes
 Usage   : $db->search_notes($search_term,$max_results)
 Function: full-text search on features, ENSEMBL-style
 Returns : an array of [$name,$description,$score]
 Args    : see below
 Status  : public

This routine performs a full-text search on feature attributes (which
attributes depend on implementation) and returns a list of
[$name,$description,$score], where $name is the feature ID,
$description is a human-readable description such as a locus line, and
$score is the match strength.

Since this is a decidedly non-standard thing to do (but the generic
genome browser uses it), the default method returns an empty list.
You do not have to implement it.

=cut

sub search_notes { return }

=head2 refclass

 Title   : refclass
 Usage   : $class = $db->refclass
 Function: returns the default class to use for segment() calls
 Returns : a string
 Args    : none
 Status  : public

For data sources which use namespaces to distinguish reference
sequence accessions, this returns the default namespace (or "class")
to use.  This interface defines a default of "Accession".

=cut

sub refclass { "Accession" }

1;
