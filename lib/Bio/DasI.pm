#
# BioPerl module for Bio::DasI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

  bioperl-l@bioperl.org

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Lincoln Stein

Email lstein@cshl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...

package Bio::DasI;
use strict;

use Bio::Das::SegmentI;
# Object preamble - inherits from Bio::Root::Root;
use base qw(Bio::Root::RootI Bio::SeqFeature::CollectionI);

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

=head2 types

 Title   : types
 Usage   : $db->types(@args)
 Function: return list of feature types in database
 Returns : a list of Bio::Das::FeatureTypeI objects
 Args    : see below

This routine returns a list of feature types known to the database. It
is also possible to find out how many times each feature occurs.

Arguments are -option=E<gt>value pairs as follows:

  -enumerate  if true, count the features

The returned value will be a list of Bio::Das::FeatureTypeI objects
(see L<Bio::Das::FeatureTypeI>.

If -enumerate is true, then the function returns a hash (not a hash
reference) in which the keys are the stringified versions of
Bio::Das::FeatureTypeI and the values are the number of times each
feature appears in the database.

=cut

sub types {  shift->throw_not_implemented; }

=head2 parse_types

 Title   : parse_types
 Usage   : $db->parse_types(@args)
 Function: parses list of types
 Returns : an array ref containing ['method','source'] pairs
 Args    : a list of types in 'method:source' form
 Status  : internal

This method takes an array of type names in the format "method:source"
and returns an array reference of ['method','source'] pairs.  It will
also accept a single argument consisting of an array reference with
the list of type names.

=cut

# turn feature types in the format "method:source" into a list of [method,source] refs
sub parse_types {
  my $self  = shift;
  return []   if !@_ or !defined($_[0]);
  return $_[0] if ref $_[0] eq 'ARRAY' && ref $_[0][0];
  my @types = ref($_[0]) ? @{$_[0]} : @_;
  my @type_list = map { [split(':',$_,2)] } @types;
  return \@type_list;
}

=head2 segment

 Title   : segment
 Usage   : $db->segment(@args);
 Function: create a segment object
 Returns : segment object(s)
 Args    : see below

This method generates a Bio::Das::SegmentI object (see
L<Bio::Das::SegmentI>).  The segment can be used to find overlapping
features and the raw sequence.

When making the segment() call, you specify the ID of a sequence
landmark (e.g. an accession number, a clone or contig), and a
positional range relative to the landmark.  If no range is specified,
then the entire region spanned by the landmark is used to generate the
segment.

Arguments are -option=E<gt>value pairs as follows:

 -name         ID of the landmark sequence.

 -class        A namespace qualifier.  It is not necessary for the
               database to honor namespace qualifiers, but if it
               does, this is where the qualifier is indicated.

 -version      Version number of the landmark.  It is not necessary for
               the database to honor versions, but if it does, this is
               where the version is indicated.

 -start        Start of the segment relative to landmark.  Positions
               follow standard 1-based sequence rules.  If not specified,
               defaults to the beginning of the landmark.

 -end          End of the segment relative to the landmark.  If not specified,
               defaults to the end of the landmark.

The return value is a list of Bio::Das::SegmentI objects.  If the method
is called in a scalar context and there are no more than one segments
that satisfy the request, then it is allowed to return the segment.
Otherwise, the method must throw a "multiple segment exception".

=cut

#'

sub segment { shift->throw_not_implemented }

=head2 features

 Title   : features
 Usage   : $db->features(@args)
 Function: get all features, possibly filtered by type
 Returns : a list of Bio::SeqFeatureI objects
 Args    : see below
 Status  : public

This routine will retrieve features in the database regardless of
position.  It can be used to return all features, or a subset based on
their type

Arguments are -option=E<gt>value pairs as follows:

  -types     List of feature types to return.  Argument is an array
             of Bio::Das::FeatureTypeI objects or a set of strings
             that can be converted into FeatureTypeI objects.

  -callback   A callback to invoke on each feature.  The subroutine
              will be passed each Bio::SeqFeatureI object in turn.

  -attributes A hash reference containing attributes to match.

The -attributes argument is a hashref containing one or more attributes
to match against:

  -attributes => { Gene => 'abc-1',
                   Note => 'confirmed' }

Attribute matching is simple exact string matching, and multiple
attributes are ANDed together.  See L<Bio::DB::ConstraintsI> for a
more sophisticated take on this.

If one provides a callback, it will be invoked on each feature in
turn.  If the callback returns a false value, iteration will be
interrupted.  When a callback is provided, the method returns undef.

=cut

sub features { shift->throw_not_implemented }

=head2 get_feature_by_name

 Title   : get_feature_by_name
 Usage   : $db->get_feature_by_name(-class=>$class,-name=>$name)
 Function: fetch features by their name
 Returns : a list of Bio::SeqFeatureI objects
 Args    : the class and name of the desired feature
 Status  : public

This method can be used to fetch named feature(s) from the database.
The -class and -name arguments have the same meaning as in segment(),
and the method also accepts the following short-cut forms:

  1) one argument: the argument is treated as the feature name
  2) two arguments: the arguments are treated as the class and name
     (note: this uses _rearrange() so the first argument must not
     begin with a hyphen or it will be interpreted as a named
     argument).

This method may return zero, one, or several Bio::SeqFeatureI objects.
The implementor may allow the name to contain wildcards, in which case
standard C-shell glob semantics are expected.

=cut

sub get_feature_by_name {
  shift->throw_not_implemented();
}

=head2 get_feature_by_target

 Title   : get_feature_by_target
 Usage   : $db->get_feature_by_target($class => $name)
 Function: fetch features by their similarity target
 Returns : a list of Bio::SeqFeatureI objects
 Args    : the class and name of the desired feature
 Status  : public

This method can be used to fetch a named feature from the database
based on its similarity hit.  The arguments are the same as
get_feature_by_name().  If this is not implemented, the interface
defaults to using get_feature_by_name().

=cut

sub get_feature_by_target {
  shift->get_feature_by_name(@_);
}

=head2 get_feature_by_id

 Title   : get_feature_by_id
 Usage   : $db->get_feature_by_target($id)
 Function: fetch a feature by its ID
 Returns : a Bio::SeqFeatureI objects
 Args    : the ID of the feature
 Status  : public

If the database provides unique feature IDs, this can be used to
retrieve a single feature from the database.  If not overridden, this
interface calls get_feature_by_name() and returns the first element.

=cut

sub get_feature_by_id {
  (shift->get_feature_by_name(@_))[0];
}

=head2 get_feature_by_attribute

 Title   : get_feature_by_attribute
 Usage   : $db->get_feature_by_attribute(attribute1=>value1,attribute2=>value2)
 Function: fetch features by combinations of attribute values
 Returns : a list of Bio::SeqFeatureI objects
 Args    : the class and name of the desired feature
 Status  : public

This method can be used to fetch a set of features from the database.
Attributes are a list of name=E<gt>value pairs.  They will be
logically ANDed together.  If an attribute value is an array
reference, the list of values in the array is treated as an
alternative set of values to be ORed together.

=cut

sub get_feature_by_attribute {
  shift->throw_not_implemented();
}


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

=head2 get_seq_stream

 Title   : get_seq_stream
 Usage   : $seqio = $db->get_seq_stream(@args)
 Function: Performs a query and returns an iterator over it
 Returns : a Bio::SeqIO stream capable of returning Bio::SeqFeatureI objects
 Args    : As in features()
 Status  : public

This routine takes the same arguments as features(), but returns a
Bio::SeqIO::Stream-compliant object.  Use it like this:

  $stream = $db->get_seq_stream('exon');
  while (my $exon = $stream->next_seq) {
     print $exon,"\n";
  }

NOTE: In the interface this method is aliased to get_feature_stream(),
as the name is more descriptive.

=cut

sub get_seq_stream { shift->throw_not_implemented }
sub get_feature_stream {shift->get_seq_stream(@_) }

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
