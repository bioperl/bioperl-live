#
# BioPerl module for Bio::Das::SegmentI
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

Bio::Das::SegmentI - DAS-style access to a feature database

=head1 SYNOPSIS

  # Get a Bio::Das::SegmentI object from a Bio::DasI database...

  $segment = $das->segment(-name=>'Landmark',
                           -start=>$start,
                           -end => $end);

  @features = $segment->overlapping_features(-type=>['type1','type2']);
  # each feature is a Bio::SeqFeatureI-compliant object

  @features = $segment->contained_features(-type=>['type1','type2']);

  @features = $segment->contained_in(-type=>['type1','type2']);

  $stream = $segment->get_feature_stream(-type=>['type1','type2','type3'];
  while (my $feature = $stream->next_seq) {
     # do something with feature
  }

  $count = $segment->features_callback(-type=>['type1','type2','type3'],
                                       -callback => sub { ... { }
                                       );

=head1 DESCRIPTION

Bio::Das::SegmentI is a simplified alternative interface to sequence
annotation databases used by the distributed annotation system. In
this scheme, the genome is represented as a series of landmarks.  Each
Bio::Das::SegmentI object ("segment") corresponds to a genomic region
defined by a landmark and a start and end position relative to that
landmark.  A segment is created using the Bio::DasI segment() method.

Features can be filtered by the following attributes:

  1) their location relative to the segment (whether overlapping,
          contained within, or completely containing)

  2) their type

  3) other attributes using tag/value semantics

Access to the feature list uses three distinct APIs:

  1) fetching entire list of features at a time

  2) fetching an iterator across features

  3) a callback

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bio.perl.org

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Lincoln Stein

Email lstein@cshl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...

package Bio::Das::SegmentI;
use strict;


# Object preamble - inherits from Bio::Root::RootI;
use base qw(Bio::Root::RootI);

=head2 seq_id

 Title   : seq_id
 Usage   : $ref = $s->seq_id
 Function: return the ID of the landmark
 Returns : a string
 Args    : none
 Status  : Public

=cut

sub seq_id { shift->throw_not_implemented }

=head2 display_name

 Title   : seq_name
 Usage   : $ref = $s->seq_name
 Function: return the human-readable name for the landmark
 Returns : a string
 Args    : none
 Status  : Public

This defaults to the same as seq_id.

=cut

sub display_name { shift->seq_id }

=head2 start

 Title   : start
 Usage   : $s->start
 Function: start of segment
 Returns : integer
 Args    : none
 Status  : Public

This is a read-only accessor for the start of the segment.  Alias
to low() for Gadfly compatibility.

=cut

sub start  { shift->throw_not_implemented }
sub low    { shift->start }

=head2 end

 Title   : end
 Usage   : $s->end
 Function: end of segment
 Returns : integer
 Args    : none
 Status  : Public

This is a read-only accessor for the end of the segment. Alias to
high() for Gadfly compatibility.

=cut

sub end   { shift->throw_not_implemented  }
sub stop  { shift->end }
sub high  { shift->end }

=head2 length

 Title   : length
 Usage   : $s->length
 Function: length of segment
 Returns : integer
 Args    : none
 Status  : Public

Returns the length of the segment.  Always a positive number.

=cut

sub length { shift->throw_not_implemented; }

=head2 seq

 Title   : seq
 Usage   : $s->seq
 Function: get the sequence string for this segment
 Returns : a string
 Args    : none
 Status  : Public

Returns the sequence for this segment as a simple string.

=cut

sub seq {shift->throw_not_implemented}

=head2 ref

 Title   : ref
 Usage   : $ref = $s->ref([$newlandmark])
 Function: get/set the reference landmark for addressing
 Returns : a string
 Args    : none
 Status  : Public

This method is used to examine/change the reference landmark used to
establish the coordinate system.  By default, the landmark cannot be
changed and therefore this has the same effect as seq_id().  The new
landmark might be an ID, or another Das::SegmentI object.

=cut

sub ref    { shift->seq_id }
*refseq = \&ref;

=head2 absolute

 Title   : absolute
 Usage   : $s->absolute([$new_value])
 Function: get/set absolute addressing mode
 Returns : flag
 Args    : new flag (optional)
 Status  : Public

Turn on and off absolute-addressing mode.  In absolute addressing
mode, coordinates are relative to some underlying "top level"
coordinate system (such as a chromosome). ref() returns the identity
of the top level landmark, and start() and end() return locations
relative to that landmark.  In relative addressing mode, coordinates
are relative to the landmark sequence specified at the time of segment
creation or later modified by the ref() method.

The default is to return false and to do nothing in response to
attempts to set absolute addressing mode.

=cut

sub absolute { return }

=head2 features

 Title   : features
 Usage   : @features = $s->features(@args)
 Function: get features that overlap this segment
 Returns : a list of Bio::SeqFeatureI objects
 Args    : see below
 Status  : Public

This method will find all features that intersect the segment in a
variety of ways and return a list of Bio::SeqFeatureI objects.  The
feature locations will use coordinates relative to the reference
sequence in effect at the time that features() was called.

The returned list can be limited to certain types, attributes or
range intersection modes.  Types of range intersection are one of:

   "overlaps"      the default
   "contains"      return features completely contained within the segment
   "contained_in"  return features that completely contain the segment

Two types of argument lists are accepted.  In the positional argument
form, the arguments are treated as a list of feature types.  In the
named parameter form, the arguments are a series of -name=E<gt>value
pairs.

  Argument    Description
  --------   ------------

  -types      An array reference to type names in the format
	      "method:source"

  -attributes A hashref containing a set of attributes to match

  -rangetype  One of "overlaps", "contains", or "contained_in".

  -iterator   Return an iterator across the features.

  -callback   A callback to invoke on each feature

The -attributes argument is a hashref containing one or more
attributes to match against:

  -attributes => { Gene => 'abc-1',
                   Note => 'confirmed' }

Attribute matching is simple string matching, and multiple attributes
are ANDed together.  More complex filtering can be performed using the
-callback option (see below).

If -iterator is true, then the method returns an object reference that
implements the next_seq() method.  Each call to next_seq() returns a
new Bio::SeqFeatureI object.

If -callback is passed a code reference, the code reference will be
invoked on each feature returned.  The code will be passed two
arguments consisting of the current feature and the segment object
itself, and must return a true value. If the code returns a false
value, feature retrieval will be aborted.

-callback and -iterator are mutually exclusive options.  If -iterator
is defined, then -callback is ignored.

NOTE: the following methods all build on top of features(), and do not
need to be explicitly implemented.

    overlapping_features()
    contained_features()
    contained_in()
    get_feature_stream()

=cut

sub features {shift->throw_not_implemented}

=head2 overlapping_features

 Title   : overlapping_features
 Usage   : @features = $s->overlapping_features(@args)
 Function: get features that overlap this segment
 Returns : a list of Bio::SeqFeatureI objects
 Args    : see below
 Status  : Public

This method is identical to features() except that it defaults to
finding overlapping features.

=cut

sub overlapping_features {
  my $self = shift;
  my @args = $_[0] =~ /^-/ ? (@_,         -rangetype=>'overlaps')
                           : (-types=>\@_,-rangetype=>'overlaps');
  $self->features(@args);
}

=head2 contained_features

 Title   : contained_features
 Usage   : @features = $s->contained_features(@args)
 Function: get features that are contained in this segment
 Returns : a list of Bio::SeqFeatureI objects
 Args    : see below
 Status  : Public

This method is identical to features() except that it defaults to
a range type of 'contained'.

=cut

sub contained_features {
  my $self = shift;
  my @args = $_[0] =~ /^-/ ? (@_,         -rangetype=>'contained')
                           : (-types=>\@_,-rangetype=>'contained');
  $self->features(@args);
}

=head2 contained_in

 Title   : contained_in
 Usage   : @features = $s->contained_in(@args)
 Function: get features that contain this segment
 Returns : a list of Bio::SeqFeatureI objects
 Args    : see below
 Status  : Public

This method is identical to features() except that it defaults to
a range type of 'contained_in'.

=cut

sub contained_in {
  my $self = shift;
  my @args = $_[0] =~ /^-/ ? (@_,         -rangetype=>'contained_in')
                           : (-types=>\@_,-rangetype=>'contained_in');
  $self->features(@args);
}

=head2 get_feature_stream

 Title   : get_feature_stream
 Usage   : $iterator = $s->get_feature_stream(@args)
 Function: get an iterator across the segment
 Returns : an object that implements next_seq()
 Args    : see below
 Status  : Public

This method is identical to features() except that it always generates
an iterator.

NOTE: This is defined in the interface in terms of features().  You do not
have to implement it.

=cut

sub get_feature_stream {
  my $self = shift;
  my @args = defined $_[0] && $_[0] =~ /^-/ ? (@_,         -iterator=>1)
                                            : (-types=>\@_,-iterator=>1);
  $self->features(@args);
}

=head2 factory

 Title   : factory
 Usage   : $factory = $s->factory
 Function: return the segment factory
 Returns : a Bio::DasI object
 Args    : see below
 Status  : Public

This method returns a Bio::DasI object that can be used to fetch
more segments.  This is typically the Bio::DasI object from which
the segment was originally generated.

=cut

#'

sub factory {shift->throw_not_implemented}

=head2 primary_tag

 Title   : primary_tag
 Usage   : $tag = $s->primary_tag
 Function: identifies the segment as type "DasSegment"
 Returns : a string named "DasSegment"
 Args    : none
 Status  : Public, but see below

This method provides Bio::Das::Segment objects with a primary_tag()
field that identifies them as being of type "DasSegment".  This allows
the Bio::Graphics engine to render segments just like a feature in order
nis way useful.

This does not need to be implemented.  It is defined by the interface.

=cut

#'

sub primary_tag {"DasSegment"}

=head2 strand

 Title   : strand
 Usage   : $strand = $s->strand
 Function: identifies the segment strand as 0
 Returns : the number 0
 Args    : none
 Status  : Public, but see below

This method provides Bio::Das::Segment objects with a strand() field
that identifies it as being strandless.  This allows the Bio::Graphics
engine to render segments just like a feature in order nis way useful.

This does not need to be implemented.  It is defined by the interface.

=cut

sub strand      { 0 }

1;
