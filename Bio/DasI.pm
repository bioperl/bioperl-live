# $Id$
#
# BioPerl module for Bio::DasI
#
# Cared for by Lincoln Stein <lstein@cshl.org>
#
# Copyright Lincoln Stein
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DasI - DAS-style access to a feature database

=head1 SYNOPSIS

  # Open up a feature database somehow...
  $db = Bio::DasI->new(@args);

  @segments = $db->segment(-name  => 'NT_29921.4',
                           -start => 1,
			   -end   => 1000000);

  # segments are Bio::Das::SegmentI - compliant objects

  # fetch a list of features
  @features = $db->features(-type=>['type1','type2','type3']);

  # invoke a callback over features
  $db->features(-type=>['type1','type2','type3'],
                -callback => sub { ... }
		);

  $stream   = $db->get_seq_stream(-type=>['type1','type2','type3']);
  while (my $feature = $stream->next_seq) {
     # each feature is a Bio::SeqFeatureI-compliant object
  }

  # get all feature types
  @types   = $db->types;

  # count types
  %types   = $db->types(-enumerate=>1);

  @feature = $db->get_feature_by_name($class=>$name);
  @feature = $db->get_feature_by_target($target_name);
  @feature = $db->get_feature_by_attribute($att1=>$value1,$att2=>$value2);
  $feature = $db->get_feature_by_id($id);

  $error = $db->error;

=head1 DESCRIPTION

Bio::DasI is a simplified alternative interface to sequence annotation
databases used by the distributed annotation system (see
L<Bio::Das>). In this scheme, the genome is represented as a series of
features, a subset of which are named.  Named features can be used as
reference points for retrieving "segments" (see L<Bio::Das::SegmentI>),
and these can, in turn, be used as the basis for exploring the genome
further.

In addition to a name, each feature has a "class", which is
essentially a namespace qualifier and a "type", which describes what
type of feature it is.  Das uses the GO consortium's ontology of
feature types, and so the type is actually an object of class
Bio::Das::FeatureTypeI (see L<Bio::Das::FeatureTypeI>). Bio::DasI
provides methods for interrogating the database for the types it
contains and the counts of each type.

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

package Bio::DasI;
use strict;

use vars qw( @ISA );
use Bio::DB::SegmentProviderI;
use Bio::DB::AggregatingFeatureProviderI;
@ISA = qw( Bio::DB::SegmentProviderI Bio::DB::AggregatingFeatureProviderI );

=head2 new

 Title   : new
 Usage   : Bio::DasI->new(@args)
 Function: Create new Bio::DasI object
 Returns : a Bio::DasI object
 Args    : see below

The new() method creates a new object.  The argument list is either a
single argument consisting of a connection string, or the following
list of -name=E<gt>value arguments:

   Argument        Description
   --------        -----------

   -dsn            Connection string for database
   -adaptor        Name of an adaptor class to use when connecting
   -aggregator     Array ref containing list of aggregators
                     "semantic mappers" to apply to database
   -user           Authentication username
   -pass           Authentication password

Implementors of DasI may add other arguments.

=cut

sub new {shift->throw_not_implemented}

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
