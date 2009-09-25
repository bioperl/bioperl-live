=head1 NAME

Bio::DB::GFF::RelSegment -- Sequence segment with relative coordinate support

=head1 SYNOPSIS

See L<Bio::DB::GFF>.

=head1 DESCRIPTION

Bio::DB::GFF::RelSegment is a stretch of sequence that can handle
relative coordinate addressing.  It inherits from
Bio::DB::GFF::Segment, and is the base class for
Bio::DB::GFF::Feature.

In addition to the source sequence, a relative segment has a
"reference sequence", which is used as the basis for its coordinate
system.  The reference sequence can be changed at will, allowing you
freedom to change the "frame of reference" for features contained
within the segment.  For example, by setting a segment's reference
sequence to the beginning of a gene, you can view all other features
in gene-relative coordinates.

The reference sequence and the source sequence must be on the same
physical stretch of DNA, naturally.  However, they do not have to be
on the same strand.  The strandedness of the reference sequence
determines whether coordinates increase to the right or the left.

Generally, you will not create or manipulate Bio::DB::GFF::RelSeg0ment
objects directly, but use those that are returned by the Bio::DB::GFF
module.

=head2 An Example

To understand how relative coordinates work, consider the following
example from the C. elegans database.  First we create the appropriate
GFF accessor object (the factory):

   my $db = Bio::DB::GFF->new(-dsn => 'dbi:mysql:elegans',
                              -adaptor=>'dbi:mysqlopt');

Now we fetch out a segment based on cosmid clone ZK909:

  my $seg = $db->segment('ZK909');

If we call the segment's refseq() method, we see that the base of the
coordinate system is the sequence "ZK154", and that its start and
stop positions are 1 and the length of the cosmid:

  print $seg->refseq;
  => ZK909

  print $seg->start,' - ',$seg->stop;
  => 1 - 33782

As a convenience, the "" operator is overloaded in this class, to give
the reference sequence, and start and stop positions:

  print $seg;
  => ZK909:1,33782

Internally, Bio::DB::GFF::RelSegment has looked up the absolute
coordinates of this segment and maintains the source sequence and the
absolute coordinates relative to the source sequence.  We can see this 
information using sourceseq() (inherited from Bio::DB::GFF::Segment)
and the abs_start() and abs_end() methods:

  print $seg->sourceseq;
  => CHROMOSOME_I

  print $seg->abs_start,' - ',$seg->abs_end;
  => 14839545 - 14873326

We can also put the segment into absolute mode, so that it behaves
like Bio::DB::Segment, and always represents coordinates on the source
sequence.  This is done by passing a true value to the absolute()
method:

  $seq->absolute(1);
  print $seg;
  => CHROMOSOME_I:14839545,14873326

We can change the reference sequence at any time.  One way is to call
the segment's ref() method, giving it the ID (and optionally the
class) of another landmark on the genome.  For example, if we know
that cosmid ZK337 is adjacent to ZK909, then we can view ZK909 in
ZK337-relative coordinates:

  $seg->refseq('ZK337');
  print $seg;
  => ZK337:-33670,111

We can call the segment's features() method in order to get the list
of contigs that overlap this segment (in the C. elegans database,
contigs have feature type "Sequence:Link"):

  @links = $seg->features('Sequence:Link');

We can now set the reference sequence to the first of these contigs like so:

  $seg->refseq($links[0]);
  print $seg;
  => Sequence:Link(LINK_Y95D11A):3997326,4031107

=cut

package Bio::DB::GFF::RelSegment;

use strict;

use Bio::DB::GFF::Feature;
use Bio::DB::GFF::Util::Rearrange;
use Bio::RangeI;

use base qw(Bio::DB::GFF::Segment);

use overload '""' => 'asString',
             'bool' => sub { overload::StrVal(shift) },
             fallback=>1;

=head1 API

The remainder of this document describes the API for
Bio::DB::GFF::Segment.

=cut

=head2 new

 Title   : new
 Usage   : $s = Bio::DB::GFF::RelSegment->new(@args)
 Function: create a new relative segment
 Returns : a new Bio::DB::GFF::RelSegment object
 Args    : see below
 Status  : Public

This method creates a new Bio::DB::GFF::RelSegment object.  Generally
this is called automatically by the Bio::DB::GFF module and
derivatives.

This function uses a named-argument style:

 -factory      a Bio::DB::GFF::Adaptor to use for database access
 -seq          ID of the source sequence
 -class        class of the source sequence
 -start        start of the desired segment relative to source sequence
 -stop         stop of the desired segment relative to source sequence
 -ref          ID of the reference sequence
 -refclass     class of the reference sequence
 -offset       0-based offset from source sequence to start of segment
 -length       length of desired segment
 -absolute, -force_absolute
               use absolute coordinates, rather than coordinates relative
               to the start of self or the reference sequence

The -seq argument accepts the ID of any landmark in the database.  The
stored source sequence becomes whatever the GFF file indicates is the
proper sequence for this landmark.  A class of "Sequence" is assumed
unless otherwise specified in the -class argument.

If the argument to -seq is a Bio::GFF::Featname object (such as
returned by the group() method), then the class is taken from that.

The optional -start and -stop arguments specify the end points for the
retrieved segment.  For those who do not like 1-based indexing,
-offset and -length are provided.  If both -start/-stop and
-offset/-length are provided, the latter overrides the former.
Generally it is not a good idea to mix metaphors.

-ref and -refclass together indicate a sequence to be used for
relative coordinates.  If not provided, the source sequence indicated
by -seq is used as the reference sequence.  If the argument to -ref is
a Bio::GFF::Featname object (such as returned by the group() method),
then the class is taken from that.

-force_absolute should be used if you wish to skip the lookup of the
absolute position of the source sequence that ordinarily occurs when
you create a relative segment.  In this case, the source sequence must
be a sequence that has been specified as the "source" in the GFF file.

=cut

# Create a new Bio::DB::GFF::RelSegment Object
# arguments are:
#      -factory    => factory and DBI interface
#      -seq        => $sequence_name
#      -start      => $start_relative_to_sequence
#      -stop       => $stop_relative_to_sequence
#      -ref        => $sequence which establishes coordinate system
#      -offset     => 0-based offset relative to sequence
#      -length     => length of segment
#      -nocheck    => turn off checking, force segment to be constructed
#      -absolute   => use absolute coordinate addressing

sub new {
  my $package = shift;
  my ($factory,$name,$start,$stop,$refseq,$class,$refclass,$offset,$length,$force_absolute,$nocheck) =
    rearrange([
	       'FACTORY',
	       [qw(NAME SEQ SEQUENCE SOURCESEQ)],
	       [qw(START BEGIN)],
	       [qw(STOP END)],
	       [qw(REFSEQ REF REFNAME)],
	       [qw(CLASS SEQCLASS)],
	       qw(REFCLASS),
	       [qw(OFFSET OFF)],
	       [qw(LENGTH LEN)],
	       [qw(ABSOLUTE)],
	       [qw(NOCHECK FORCE)],
	     ],@_);

  $package = ref $package if ref $package;
  $factory or $package->throw("new(): provide a -factory argument");

  # to allow people to use segments as sources
  if (ref($name) && $name->isa('Bio::DB::GFF::Segment')) {
    $start = 1              unless defined $start;
    $stop  = $name->length  unless defined $stop;
    return $name->subseq($start,$stop);
  }

  my @object_results;

  # support for Featname objects
  if (ref($name) && $name->can('class')) {
    $class = $name->class;
    $name  = $name->name;
  }

  # if the class of the landmark is not specified then default to 'Sequence'
  $class ||= eval{$factory->default_class} || 'Sequence';

  # confirm that indicated sequence is actually in the database!
  my @abscoords;

  # abscoords() will now return an array ref, each element of which is
  # ($absref,$absclass,$absstart,$absstop,$absstrand)

  if ($nocheck) {
    $force_absolute++;
    $start = 1;
  }

#  if ($force_absolute && defined($start)) { # absolute position is given to us
#    @abscoords = ([$name,$class,$start,$stop,'+']);
#  } else {
    my $result = $factory->abscoords($name,$class,$force_absolute ? $name : ()) or return;
    @abscoords = @$result;
#  }

  foreach (@abscoords) {
    my ($absref,$absclass,$absstart,$absstop,$absstrand,$sname) = @$_;
    $sname = $name unless defined $sname;
    my ($this_start,$this_stop,$this_length) = ($start,$stop,$length);

    # partially fill in object
    my $self = bless { factory => $factory },$package;

    $absstrand ||= '+';

    if ($absstart > $absstop) { # AAARGH!  DATA FORMAT ERROR!  FIX.
	($absstart,$absstop) = ($absstop,$absstart);
	$absstrand = $absstrand eq '+' ? '-' : '+';
    }

    # an explicit length overrides start and stop
    if (defined $offset) {
      warn "new(): bad idea to call new() with both a start and an offset"
	if defined $this_start;
      $this_start = $offset+1;
    }
    if (defined $this_length) {
      warn "new(): bad idea to call new() with both a stop and a length"
	if defined $this_stop;
      $this_stop = $this_start + $length - 1;
    }

    # this allows a SQL optimization way down deep
    $self->{whole}++ if $absref eq $sname and !defined($this_start) and !defined($this_stop);

    $this_start     = 1                    if !defined $this_start;
    $this_stop      = $absstop-$absstart+1 if !defined $this_stop;
    $this_length = $this_stop - $this_start + 1;

    # now offset to correct subsegment based on desired start and stop
    if ($force_absolute) {
#      ($this_start,$this_stop) = ($absstart,$absstop);
      $self->absolute(1);
    } elsif ($absstrand eq '+') {
      $this_start =  $absstart   + $this_start - 1;
      $this_stop  =  $this_start + $this_length - 1;
    } else {
      $this_start =  $absstop - ($this_start - 1);
      $this_stop  =  $absstop - ($this_stop - 1);
    }

    # handle truncation in either direction
    # This only happens if the segment runs off the end of
    # the reference sequence
    if ($factory->strict_bounds_checking &&
	(($this_start < $absstart) || ($this_stop > $absstop))) {
      # return empty if we are completely off the end of the ref se
      next unless $this_start<=$absstop && $this_stop>=$absstart;
      if (my $a = $factory->abscoords($absref,'Sequence')) {
	my $refstart = $a->[0][2];
	my $refstop  = $a->[0][3];
	if ($this_start < $refstart) {
	  $this_start = $refstart;
	  $self->{truncated}{start}++;
	}
	if ($this_stop > $refstop) {
	  $this_stop = $absstop;
	  $self->{truncated}{stop}++;
	}
      }
    }

    @{$self}{qw(sourceseq start stop strand class)}
      = ($absref,$this_start,$this_stop,$absstrand,$absclass);

    # handle reference sequence
    if (defined $refseq) {
      $refclass = $refseq->class if $refseq->can('class');
      $refclass ||= 'Sequence';
      my ($refref,$refstart,$refstop,$refstrand) = $factory->abscoords($refseq,$refclass);
      unless ($refref eq $absref) {
	$self->error("reference sequence is on $refref but source sequence is on $absref");
	return;
      }
      $refstart = $refstop if $refstrand eq '-';
      @{$self}{qw(ref refstart refstrand)} = ($refseq,$refstart,$refstrand);
    } else {
      $absstart = $absstop if $absstrand eq '-';
      @{$self}{qw(ref refstart refstrand)} = ($sname,$absstart,$absstrand);
    }
    push @object_results,$self;
  }

  return wantarray ? @object_results : $object_results[0];
}

# overridden methods
# start, stop, length
sub start {
  my $self = shift;
  return $self->strand < 0 ? $self->{stop} : $self->{start} if $self->absolute;
  $self->_abs2rel($self->{start});
}
sub end {
  my $self = shift;
  return $self->strand < 0 ? $self->{start} : $self->{stop} if $self->absolute;
  $self->_abs2rel($self->{stop});
}
*stop = \&end;

sub length {
  my $self = shift;
  return unless defined $self->abs_end;
  abs($self->abs_end - $self->abs_start) + 1;
}

sub abs_start {
  my $self = shift;
  if ($self->absolute) {
    my ($a,$b) = ($self->SUPER::abs_start,$self->SUPER::abs_end);
    return ($a<$b) ? $a : $b;
  }
  else {
    return $self->SUPER::abs_start(@_);
  }
}
sub abs_end {
  my $self = shift;
  if ($self->absolute) {
    my ($a,$b) = ($self->SUPER::abs_start,$self->SUPER::abs_end);
    return ($a>$b) ? $a : $b;
  }

  else {
    return $self->SUPER::abs_end(@_);
  }
}

*abs_stop = \&abs_end;

=head2 refseq

 Title   : refseq
 Usage   : $ref = $s->refseq([$newseq] [,$newseqclass])
 Function: get/set reference sequence
 Returns : current reference sequence
 Args    : new reference sequence and class (optional)
 Status  : Public

This method will get or set the reference sequence.  Called with no
arguments, it returns the current reference sequence.  Called with
either a sequence ID and class, a Bio::DB::GFF::Segment object (or
subclass) or a Bio::DB::GFF::Featname object, it will set the current
reference sequence and return the previous one.

The method will generate an exception if you attempt to set the
reference sequence to a sequence that isn't contained in the database,
or one that has a different source sequence from the segment.

=cut

#'
sub refseq {
  my $self = shift;
  my $g    = $self->{ref};
  if (@_) {
    my ($newref,$newclass);
    if (@_ == 2) {
      $newclass = shift;
      $newref   = shift;
    } else {
      $newref   = shift;
      $newclass = 'Sequence';
    }

    defined $newref or $self->throw('refseq() called with an undef reference sequence');

    # support for Featname objects
    $newclass = $newref->class if ref($newref) && $newref->can('class');

    # $self->throw("Cannot define a segment's reference sequence in terms of itself!")
    # if ref($newref) and overload::StrVal($newref) eq overload::StrVal($self);

    my ($refsource,undef,$refstart,$refstop,$refstrand);
    if ($newref->isa('Bio::DB::GFF::RelSegment')) {
      ($refsource,undef,$refstart,$refstop,$refstrand) =
	($newref->sourceseq,undef,$newref->abs_start,$newref->abs_end,$newref->abs_strand >= 0 ? '+' : '-');
    } else {
      my $coords = $self->factory->abscoords($newref,$newclass);
      foreach (@$coords) { # find the appropriate one
	($refsource,undef,$refstart,$refstop,$refstrand) = @$_;
	last if $refsource eq $self->{sourceseq};
      }
	
    }
    $self->throw("can't set reference sequence: $newref and $self are on different sequence segments")
      unless $refsource eq $self->{sourceseq};

    @{$self}{qw(ref refstart refstrand)} = ($newref,$refstart,$refstrand);
    $self->absolute(0);
  }
  return $self->absolute ? $self->sourceseq : $g;
}


=head2 abs_low

 Title   : abs_low
 Usage   : $s->abs_low
 Function: the absolute lowest coordinate of the segment
 Returns : an integer
 Args    : none
 Status  : Public

This is for GadFly compatibility, and returns the low coordinate in
absolute coordinates;

=cut

sub abs_low {
  my $self = shift;
  my ($a,$b) = ($self->abs_start,$self->abs_end);
  return ($a<$b) ? $a : $b;
}

=head2 abs_high

 Title   : abs_high
 Usage   : $s->abs_high
 Function: the absolute highest coordinate of the segment
 Returns : an integer
 Args    : none
 Status  : Public

This is for GadFly compatibility, and returns the high coordinate in
absolute coordinates;

=cut

sub abs_high {
  my $self = shift;
  my ($a,$b) = ($self->abs_start,$self->abs_end);
  return ($a>$b) ? $a : $b;
}


=head2 asString

 Title   : asString
 Usage   : $s->asString
 Function: human-readable representation of the segment
 Returns : a string
 Args    : none
 Status  : Public

This method will return a human-readable representation of the
segment.  It is the overloaded method call for the "" operator.

Currently the format is:

  refseq:start,stop

=cut

sub asString {
  my $self = shift;
  return $self->SUPER::asString if $self->absolute;
  my $label = $self->{ref};
  my $start = $self->start || '';
  my $stop  = $self->stop  || '';
  if (ref($label) && overload::StrVal($self) eq overload::StrVal($label->ref)) {
    $label = $self->abs_ref;
    $start = $self->abs_start;
    $stop  = $self->abs_end;
  }
  return "$label:$start,$stop";
}

=head2 name

 Title   : name
 Usage   : Alias for asString()

=cut

sub name { shift->asString }

=head2 absolute

 Title   : absolute
 Usage   : $abs = $s->absolute([$abs])
 Function: get/set absolute coordinates
 Returns : a boolean flag
 Args    : new setting for flag (optional)
 Status  : Public

Called with a boolean flag, this method controls whether to display
relative coordinates (relative to the reference sequence) or absolute
coordinates (relative to the source sequence).  It will return the
previous value of the setting.

=cut

sub absolute {
  my $self = shift;
  my $g = $self->{absolute};
  $self->{absolute} = shift if @_;
  $g;
}

=head2 features

 Title   : features
 Usage   : @features = $s->features(@args)
 Function: get features that overlap this segment
 Returns : a list of Bio::DB::GFF::Feature objects
 Args    : see below
 Status  : Public

This method will find all features that overlap the segment and return
a list of Bio::DB::GFF::Feature objects.  The features will use
coordinates relative to the reference sequence in effect at the time
that features() was called.

The returned list can be limited to certain types of feature by
filtering on their method and/or source.  In addition, it is possible
to obtain an iterator that will step through a large number of
features sequentially.

Arguments can be provided positionally or using the named arguments
format.  In the former case, the arguments are a list of feature types
in the format "method:source".  Either method or source can be
omitted, in which case the missing component is treated as a wildcard.
If no colon is present, then the type is treated as a method name.
Multiple arguments are ORed together.

Examples:

 @f = $s->features('exon:curated');           # all curated exons
 @f = $s->features('exon:curated','intron');  # curated exons and all introns
 @f = $s->features('similarity:.*EST.*');     # all similarities
                                              # having something to do
                                              # with ESTs

The named parameter form gives you control over a few options:

  -types      an array reference to type names in the format
	      "method:source"

  -merge     Whether to apply aggregators to the generated features (default yes)

  -rare      Turn on an optimization suitable for a relatively rare feature type,
             where it will be faster to filter by feature type first
             and then by position, rather than vice versa.

  -attributes a hashref containing a set of attributes to match

  -range_type One of 'overlapping', 'contains', or 'contained_in'

  -iterator  Whether to return an iterator across the features.

  -binsize   A true value will create a set of artificial features whose
             start and stop positions indicate bins of the given size, and
             whose scores are the number of features in the bin.  The
             class and method of the feature will be set to "bin",
             its source to "method:source", and its group to "bin:method:source".
             This is a handy way of generating histograms of feature density.

-merge is a boolean flag that controls whether the adaptor's
aggregators wll be applied to the features returned by this method.

If -iterator is true, then the method returns a single scalar value
consisting of a Bio::SeqIO object.  You can call next_seq() repeatedly
on this object to fetch each of the features in turn.  If iterator is
false or absent, then all the features are returned as a list.

The -attributes argument is a hashref containing one or more
attributes to match against:

  -attributes => { Gene => 'abc-1',
                   Note => 'confirmed' }

Attribute matching is simple string matching, and multiple attributes
are ANDed together.

=cut

#'

# return all features that overlap with this segment;
# optionally modified by a list of types to filter on
sub features {
  my $self = shift;
  my @args = $self->_process_feature_args(@_);
  return $self->factory->overlapping_features(@args);
}

=head2 get_SeqFeatures

 Title   : get_SeqFeatures
 Usage   :
 Function: returns the top level sequence features
 Returns : L<Bio::SeqFeatureI> objects
 Args    : none

Segments do not ordinarily return any subfeatures.

=cut

# A SEGMENT DOES NOT HAVE SUBFEATURES!
sub get_SeqFeatures { return }

=head2 feature_count

 Title   : feature_count
 Usage   : $seq->feature_count()
 Function: Return the number of SeqFeatures attached to a sequence
 Returns : integer representing the number of SeqFeatures
 Args    : none

This method comes through extension of Bio::FeatureHolderI. See
L<Bio::FeatureHolderI> for more information.

=cut

sub feature_count { 
    my $self = shift;
    my $ct = 0;
    my %type_counts = $self->types(-enumerate=>1);
    map { $ct += $_ } values %type_counts;
    $ct;
}

=head2 get_feature_stream

 Title   : features
 Usage   : $stream = $s->get_feature_stream(@args)
 Function: get a stream of features that overlap this segment
 Returns : a Bio::SeqIO::Stream-compliant stream
 Args    : see below
 Status  : Public

This is the same as features(), but returns a stream.  Use like this:

 $stream = $s->get_feature_stream('exon');
 while (my $exon = $stream->next_seq) {
    print $exon->start,"\n";
 }

=cut

sub get_feature_stream {
  my $self = shift;
  my @args = defined($_[0]) && $_[0] =~ /^-/ ? (@_,-iterator=>1) : (-types=>\@_,-iterator=>1);
  $self->features(@args);
}

=head2 get_seq_stream

 Title   : get_seq_stream
 Usage   : $stream = $s->get_seq_stream(@args)
 Function: get a stream of features that overlap this segment
 Returns : a Bio::SeqIO::Stream-compliant stream
 Args    : see below
 Status  : Public

This is the same as feature_stream(), and is provided for Bioperl
compatibility.  Use like this:

 $stream = $s->get_seq_stream('exon');
 while (my $exon = $stream->next_seq) {
    print $exon->start,"\n";
 }

=cut

*get_seq_stream = \&get_feature_stream;


=head2 overlapping_features

 Title   : overlapping_features
 Usage   : @features = $s->overlapping_features(@args)
 Function: get features that overlap this segment
 Returns : a list of Bio::DB::GFF::Feature objects
 Args    : see features()
 Status  : Public

This is an alias for the features() method, and takes the same
arguments.

=cut

*overlapping_features = \&features;

=head2 contained_features

 Title   : contained_features
 Usage   : @features = $s->contained_features(@args)
 Function: get features that are contained by this segment
 Returns : a list of Bio::DB::GFF::Feature objects
 Args    : see features()
 Status  : Public

This is identical in behavior to features() except that it returns
only those features that are completely contained within the segment,
rather than any that overlap.

=cut 

# return all features completely contained within this segment
sub contained_features {
  my $self = shift;
  local $self->{whole} = 0;
  my @args = $self->_process_feature_args(@_);
  return $self->factory->contained_features(@args);
}

# *contains = \&contained_features;

=head2 contained_in

 Title   : contained_in
 Usage   : @features = $s->contained_in(@args)
 Function: get features that contain this segment
 Returns : a list of Bio::DB::GFF::Feature objects
 Args    : see features()
 Status  : Public

This is identical in behavior to features() except that it returns
only those features that completely contain the segment.

=cut

# return all features completely contained within this segment
sub contained_in {
  my $self = shift;
  local $self->{whole} = 0;
  my @args = $self->_process_feature_args(@_);
  return $self->factory->contained_in(@args);
}

=head2 delete

 Title   : delete
 Usage   : $db->delete(@args)
 Function: delete features
 Returns : count of features deleted -- if available
 Args    : numerous, see below
 Status  : public

This method deletes all features that overlap the specified region or
are of a particular type.  If no arguments are provided and the -force
argument is true, then deletes ALL features.

Arguments:

 -type,-types  Either a single scalar type to be deleted, or an
               reference to an array of types.

 -range_type   Control the range type of the deletion.  One of "overlaps" (default)
               "contains" or "contained_in"

Examples:

  $segment->delete(-type=>['intron','repeat:repeatMasker']);  # remove all introns & repeats
  $segment->delete(-type=>['intron','repeat:repeatMasker']
		   -range_type => 'contains');                # remove all introns & repeats
                                                              # strictly contained in segment

IMPORTANT NOTE: This method only deletes features.  It does *NOT*
delete the names of groups that contain the deleted features.  Group
IDs will be reused if you later load a feature with the same group
name as one that was previously deleted.

NOTE ON FEATURE COUNTS: The DBI-based versions of this call return the
result code from the SQL DELETE operation.  Some dbd drivers return the
count of rows deleted, while others return 0E0.  Caveat emptor.

=cut

# return all features completely contained within this segment
sub delete {
  my $self = shift;
  my ($type,$range_type) =
    rearrange([[qw(TYPE TYPES)],'RANGE_TYPE'],@_);
  my $types = $self->factory->parse_types($type);  # parse out list of types
  $range_type ||= 'overlaps';
  return $self->factory->_delete({
                                  segments   => [$self],
                                  types      => $types,
                                  range_type => $range_type
                                  });
}

=head2 _process_feature_args

 Title   : _process_feature_args
 Usage   : @args = $s->_process_feature_args(@args)
 Function: preprocess arguments passed to features, 
           contained_features, and overlapping_features
 Returns : a list of parsed arguents
 Args    : see feature()
 Status  : Internal

This is an internal method that is used to check and format the
arguments to features() before passing them on to the adaptor.

=cut 

sub _process_feature_args {
  my $self       = shift;

  my ($ref,$class,$start,$stop,$strand,$whole)
    = @{$self}{qw(sourceseq class start stop strand whole)};

  ($start,$stop) = ($stop,$start) if defined $strand && $strand eq '-';

  my @args = (-ref=>$ref,-class=>$class);

  # indicating that we are fetching the whole segment allows certain
  # SQL optimizations.
  push @args,(-start=>$start,-stop=>$stop) unless $whole;

  if (@_) {
    if ($_[0] =~ /^-/) {
      push @args,@_;
    } else {
      my @types = @_;
      push @args,-types=>\@types;
    }
  }
  push @args,-parent=>$self;
  @args;
}

=head2 types

 Title   : types
 Usage   : @types = $s->types([-enumerate=>1])
 Function: list feature types that overlap this segment
 Returns : a list of Bio::DB::GFF::Typename objects or a hash
 Args    : see below
 Status  : Public

The types() method will return a list of Bio::DB::GFF::Typename
objects, each corresponding to a feature that overlaps the segment.
If the optional -enumerate parameter is set to a true value, then the
method will return a hash in which the keys are the type names and the 
values are the number of times a feature of that type is present on
the segment.  For example:

  %count = $s->types(-enumerate=>1);

=cut 

# wrapper for lower-level types() call.
sub types {
  my $self = shift;
  my ($ref,$class,$start,$stop,$strand) = @{$self}{qw(sourceseq class start stop strand)};
  ($start,$stop) = ($stop,$start) if $strand eq '-';

  my @args;
  if (@_ && $_[0] !~ /^-/) {
    @args = (-type => \@_)
  } else {
    @args = @_;
  }
  $self->factory->types(-ref  => $ref,
			-start=> $start,
			-stop => $stop,
			@args);
}

=head1 Internal Methods

The following are internal methods and should not be called directly.

=head2 new_from_segment

 Title   : new_from_segment
 Usage   : $s = $segment->new_from_segment(@args)
 Function: create a new relative segment
 Returns : a new Bio::DB::GFF::RelSegment object
 Args    : see below
 Status  : Internal

This constructor is used internally by the subseq() method.  It forces
the new segment into the Bio::DB::GFF::RelSegment package, regardless
of the package that it is called from.  This causes subclass-specfic
information, such as feature types, to be dropped when a subsequence
is created.

=cut

sub new_from_segment {
  my $package   = shift;
  $package      = ref $package if ref $package;
  my $segment   = shift;
  my $new = {};
  @{$new}{qw(factory sourceseq start stop strand class ref refstart refstrand)}
    = @{$segment}{qw(factory sourceseq start stop strand class ref refstart refstrand)};
  return bless $new,__PACKAGE__;
}

=head2 _abs2rel

 Title   : _abs2rel
 Usage   : @coords = $s->_abs2rel(@coords)
 Function: convert absolute coordinates into relative coordinates
 Returns : a list of relative coordinates
 Args    : a list of absolute coordinates
 Status  : Internal

This is used internally to map from absolute to relative
coordinates. It does not take the offset of the reference sequence
into account, so please use abs2rel() instead.

=cut

sub _abs2rel {
  my $self = shift;
  my @result;
  return unless defined $_[0];

  if ($self->absolute) {
    @result = @_;
  } else {
    my ($refstart,$refstrand) = @{$self}{qw(refstart refstrand)};
    @result = defined($refstrand) && $refstrand eq '-' ? map { $refstart - $_ + 1 } @_
                                                       : map { $_ - $refstart + 1 } @_;
  }
  # if called with a single argument, caller will expect a single scalar reply
  # not the size of the returned array!
  return $result[0] if @result == 1 and !wantarray;
  @result;
}

=head2 rel2abs

 Title   : rel2abs
 Usage   : @coords = $s->rel2abs(@coords)
 Function: convert relative coordinates into absolute coordinates
 Returns : a list of absolute coordinates
 Args    : a list of relative coordinates
 Status  : Public

This function takes a list of positions in relative coordinates to the
segment, and converts them into absolute coordinates.

=cut

sub rel2abs {
  my $self = shift;
  my @result;

  if ($self->absolute) {
    @result = @_;
  } else {
    my ($abs_start,$abs_strand) = ($self->abs_start,$self->abs_strand);
    @result = $abs_strand < 0 ? map { $abs_start - $_ + 1 } @_
                              : map { $_ + $abs_start - 1 } @_;
  }
  # if called with a single argument, caller will expect a single scalar reply
  # not the size of the returned array!
  return $result[0] if @result == 1 and !wantarray;
  @result;
}

=head2 abs2rel

 Title   : abs2rel
 Usage   : @rel_coords = $s->abs2rel(@abs_coords)
 Function: convert absolute coordinates into relative coordinates
 Returns : a list of relative coordinates
 Args    : a list of absolute coordinates
 Status  : Public

This function takes a list of positions in absolute coordinates
and returns a list expressed in relative coordinates.

=cut

sub abs2rel {
  my $self = shift;
  my @result;

  if ($self->absolute) {
    @result = @_;
  } else {
    my ($abs_start,$abs_strand) = ($self->abs_start,$self->abs_strand);
    @result = $abs_strand < 0 ? map { $abs_start - $_ + 1 } @_
                              : map { $_ - $abs_start + 1 } @_;
  }
  # if called with a single argument, caller will expect a single scalar reply
  # not the size of the returned array!
  return $result[0] if @result == 1 and !wantarray;
  @result;
}

sub subseq {
  my $self = shift;
  my $obj  = $self->SUPER::subseq(@_);
  bless $obj,__PACKAGE__;    # always bless into the generic RelSegment package
}

sub strand {
  my $self = shift;
  if ($self->absolute) {
    return _to_strand($self->{strand});
  }
  my $start = $self->start;
  my $stop  = $self->stop;
  return 0 unless defined $start and defined $stop;
  return $stop <=> $start;
}

sub _to_strand {
  my $s = shift;
  return -1 if $s eq '-';
  return +1 if $s eq '+';
  return 0;
}

=head2 Bio::RangeI Methods

The following Bio::RangeI methods are supported:

overlaps(), contains(), equals(),intersection(),union(),overlap_extent()

=cut

sub intersection {
  my $self     = shift;
  my (@ranges) = @_;
  unshift @ranges,$self if ref $self;
  $ranges[0]->isa('Bio::DB::GFF::RelSegment')
    or return $self->SUPER::intersection(@_);

  my $ref = $ranges[0]->abs_ref;
  my ($low,$high);
  foreach (@ranges) {
    return unless $_->can('abs_ref');
    $ref eq $_->abs_ref or return;
    $low  = $_->abs_low   if !defined($low)  or $low  < $_->abs_low;
    $high = $_->abs_high  if !defined($high) or $high > $_->abs_high;
  }

  return unless $low < $high;
  return Bio::DB::GFF::RelSegment->new(-factory => $self->factory,
				       -seq     => $ref,
				       -start   => $low,
				       -stop    => $high,
				      );
}

sub overlaps {
  my $self = shift;
  my($other,$so) = @_;
  return $self->SUPER::overlaps(@_) unless $other->isa('Bio::DB::GFF::RelSegment');
  return if $self->abs_ref ne $other->abs_ref;
  return if $self->abs_low  > $other->abs_high;
  return if $self->abs_high < $other->abs_low;
  1;
}

sub contains {
  my $self = shift;
  my($other,$so) = @_;
  return $self->SUPER::overlaps(@_) unless $other->isa('Bio::DB::GFF::RelSegment');
  return if $self->abs_ref ne $other->abs_ref;
  return unless $self->abs_low <= $other->abs_low;
  return unless $self->abs_high >= $other->abs_high;
  1;
}

sub union {
  my $self     = shift;
  my (@ranges) = @_;
  unshift @ranges,$self if ref $self;
  $ranges[0]->isa('Bio::DB::GFF::RelSegment')
    or return $self->SUPER::union(@_);

  my $ref = $ranges[0]->abs_ref;
  my ($low,$high);
  foreach (@ranges) {
    return unless $_->can('abs_ref');
    $ref eq $_->abs_ref or return;
    $low  = $_->abs_low  if !defined($low)  or $low  > $_->abs_low;
    $high = $_->abs_high if !defined($high) or $high < $_->abs_high;
  }
  $self->new(-factory=> $self->factory,
	     -seq    => $ref,
	     -start  => $low,
	     -stop   => $high);
}

sub version { 0 }


1;

__END__

=head1 BUGS

Schemas need some work.

=head1 SEE ALSO

L<bioperl>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.  

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

