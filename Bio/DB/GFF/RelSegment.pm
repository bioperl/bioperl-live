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

Generally, you will not create or manipulate Bio::DB::GFF::RelSegment
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
and the abs_start() and abs_stop() methods:

  print $seg->sourceseq;
  => CHROMOSOME_I

  print $seg->abs_start,' - ',$seg->abs_stop;
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
use Bio::DB::GFF::Segment;

use vars qw($VERSION @ISA);
@ISA = qw(Bio::DB::GFF::Segment);
$VERSION = '0.25';

use overload '""' => 'asString',
             'bool' => sub { overload::StrVal(shift) } ;

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

The -seq argument accepts the ID of any landmark in the database.  The
stored source sequence becomes whatever the GFF file indicates is the
proper sequence for this landmark.  A class of "Sequence" is assumed
unless otherwise specified in the -class argument.

The optional -start and -stop arguments specify the end points for the
retrieved segment.  For those who do not like 1-based indexing,
-offset and -length are provided.  If both -start/-stop and
-offset/-length are provided, the latter overrides the former.
Generally it is not a good idea to mix metaphors.

-ref and -refclass together indicate a sequence to be used for
relative coordinates.  If not provided, the source sequence indicated
by -seq is used as the reference sequence.

=cut

# Create a new Ace::Sequence::DBI::Segment object
# arguments are:
#      -factory    => factory and DBI interface
#      -seq        => $sequence_name
#      -start      => $start_relative_to_sequence
#      -stop       => $stop_relative_to_sequence
#      -ref        => $sequence which establishes coordinate system
#      -offset     => 0-based offset relative to sequence
#      -length     => length of segment
sub new {
  my $package = shift;
  my ($factory,$name,$start,$stop,$refseq,$class,$refclass,$offset,$length) =
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
	     ],@_);

  $package = ref $package if ref $package;
  $factory or $package->throw("new(): provide a -factory argument");

  # to allow people to use segments as sources
  if (ref($name) && $name->isa('Bio::DB::GFF::Segment')) {
    $start = 1              unless defined $start;
    $stop  = $name->length  unless defined $stop;
    return $name->subseq($start,$stop);
  }

  # partially fill in object
  my $self = bless { factory => $factory },$package;

  # if the class of the landmark is not specified then default to 'Sequence'
  $class ||= 'Sequence';

  # confirm that indicated sequence is actually in the database!
  my($absref,$absclass,$absstart,$absstop,$absstrand) = $factory->abscoords($name,$class)
    or return;

  # an explicit length overrides start and stop
  if (defined $offset) {
    warn "new(): bad idea to call new() with both a start and an offset"
      if defined $start;
    $start = $offset+1;
  }
  if (defined $length) {
    warn "new(): bad idea to call new() with both a stop and a length"
      if defined $stop;
    $stop = $start + $length - 1;
  }

  # this allows a SQL optimization way down deep
  $self->{whole}++ if $absref eq $name and !defined($start) and !defined($stop);

  $start = 1                    unless defined $start;
  $stop  = $absstop-$absstart+1 unless defined $stop;
  $length = $stop - $start + 1;

  # now offset to correct subsegment based on desired start and stop
  if ($absstrand eq '+') {
    $start =  $absstart + $start - 1;
    $stop  =  $start    + $length - 1;
  } else {
    $start =  $absstop - ($start - 1);
    $stop  =  $absstop - ($stop - 1);
  }
  @{$self}{qw(sourceseq start stop strand class)}
    = ($absref,$start,$stop,$absstrand,$absclass);

  # but what about the reference sequence?
  if (defined $refseq) {
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
    @{$self}{qw(ref refstart refstrand)} = ($name,$absstart,$absstrand);
  }

  return $self;
}

# overridden methods
# start, stop, length
sub start {
  my $self = shift;
  $self->abs2rel($self->{start});
}
sub stop {
  my $self = shift;
  $self->abs2rel($self->{stop});
}

sub length {
  my $self = shift;
  abs($self->abs_stop - $self->abs_start) + 1;
}

=head2 refseq

 Title   : refseq
 Usage   : $ref = $s->refseq([$newseq] [,$newseqclass])
 Function: get/set reference sequence
 Returns : current reference sequence
 Args    : new reference sequence and class (optional)
 Status  : Public

This method will get or set the reference sequence.  Called with no
arguments, it returns the current reference sequence.  Called with
either a sequence ID and class or a Bio::DB::GFF::Segment object (or
subclass), it will set the current reference sequence and return the
previous one.

The method will generate an exception if you attempt to set the
reference sequence to a sequence that isn't contained in the database,
or one that has a different source sequence from the segment.

=cut

sub refseq {
  my $self = shift;
  my $g    = $self->{ref};
  if (@_) {
    my $newref   = shift;
    my $newclass = shift || 'Sequence';

    $self->throw("Cannot define a segment's reference sequence in terms of itself!")
      if ref($newref) and overload::StrVal($newref) eq overload::StrVal($self);

    my ($refsource,undef,$refstart,$refstop,$refstrand);
    if ($newref->isa('Bio::DB::GFF::RelSegment')) {
      ($refsource,undef,$refstart,$refstop,$refstrand) = 
	($newref->sourceseq,undef,$newref->abs_start,$newref->abs_stop,$newref->abs_strand >= 0 ? '+' : '-');
    } else {
      ($refsource,undef,$refstart,$refstop,$refstrand) = 
	$self->factory->abscoords($newref,$newclass);
    }
    $self->throw("can't set reference sequence: $newref and $self are on different sequence segments")
      unless $refsource eq $self->{sourceseq};

    @{$self}{qw(ref refstart refstrand)} = ($newref,$refstart,$refstrand);
  }
  return $self->absolute ? $self->sourceseq : $g;
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
  my $label = $self->{absolute} ? $self->{sourceseq} : $self->{ref};
  my $start = $self->start || '';
  my $stop  = $self->stop  || '';
  return "$label:$start,$stop";
}

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
omitted, in which case the missing component is treated as a
wildcard.  Regular expressions are also allowed.  If no colon is
present, then the type is treated as a method name.  Multiple
arguments are ORed together.

Examples:

 @f = $s->features('exon:curated');           # all curated exons
 @f = $s->features('exon:curated','intron');  # curated exons and all introns
 @f = $s->features('similarity:.*EST.*');     # all similarities
                                              # having something to do
                                              # with ESTs

The named parameter form gives you control over a few options:

  -types      an array reference to type names in the format
	      "method:source"

  -merge     Whether to apply aggregators to the generated features.

  -iterator  Whether to return an iterator across the features.

-merge is a boolean flag that controls whether the adaptor's
aggregators wll be applied to the features returned by this method.

If -iterator is true, then the method returns a single scalar value
consisting of a Bio::SeqIO object.  You can call next_seq() repeatedly
on this object to fetch each of the features in turn.  If iterator is
false or absent, then all the features are returned as a list.

=cut

*features = \&overlapping_features;

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


# return all features that overlap with this segment;
# optionally modified by a list of types to filter on
sub overlapping_features {
  my $self = shift;
  my @args = $self->_process_feature_args(@_);
  return $self->factory->overlapping_features(@args);
}

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
  my @args = $self->_process_feature_args(@_);
  return $self->factory->contained_features(@args);
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
  my $self = shift;
  my ($ref,$class,$start,$stop,$strand,$whole)
    = @{$self}{qw(sourceseq class start stop strand whole)};

  ($start,$stop) = ($stop,$start) if $strand eq '-';

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
  $self->factory->types(-ref  => $ref,
			-class => $class,
			-start=> $start,
			-stop => $stop,
			@_);
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

=head2 abs2rel

 Title   : abs2rel
 Usage   : @coords = $s->abs2rel(@coords)
 Function: convert absolute coordinates into relative coordinates
 Returns : a list of relative coordinates
 Args    : a list of absolute coordinates
 Status  : Inernal

This is used internally to map from absolute to relative coordinates.

=cut 

sub abs2rel {
  my $self = shift;
  my @result;
  return unless defined $_[0];

  if ($self->absolute) {
    @result = @_;
  } else {
    my ($refstart,$refstrand) = @{$self}{qw(refstart refstrand)};
    @result = $refstrand eq '+' ? map { $_ - $refstart + 1 } @_
                                : map { $refstart - $_ + 1 } @_;
  }
  # if called with a single argument, caller will expect a single scalar reply
  # not the size of the returned array!
  return $result[0] if @result == 1 and !wantarray;
  @result;
}

=head2 rel2abs

 Title   : rel2abs
 Usage   : @coords = $s->re2absl(@coords)
 Function: convert relative coordinates into absolute coordinates
 Returns : a list of absolute coordinates
 Args    : a list of relative coordinates
 Status  : Inernal

This is used internally to map from relative to absolute coordinates.

=cut 

sub rel2abs {
  my $self = shift;
  my @result;

  if ($self->absolute) {
    @result = @_;
  } else {
    my ($refstart,$refstrand) = @{$self}{qw(refstart refstrand)};
    @result = $refstrand eq '+' ? map { $_ + $refstart + 1 } @_ 
                                : map { $refstart - $_ + 1 } @_;
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


1;
__END__

__END__

=head1 BUGS

Not completely Bio::SeqFeatureI compliant yet.

Schemas need some work.

=head1 SEE ALSO

L<bioperl>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.  

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

