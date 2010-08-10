package Bio::DB::SeqFeature::Segment;


=head1 NAME

Bio::DB::SeqFeature::Segment -- Location-based access to genome annotation data

=head1 SYNOPSIS

 use Bio::DB::SeqFeature::Store;
 # Open the sequence database
 my $db      = Bio::DB::SeqFeature::Store->new( -adaptor => 'DBI::mysql',
                                                -dsn     => 'dbi:mysql:test');
 my $segment  = $db->segment('Chr1',5000=>6000);
 my @features = $segment->features('mRNA','match');

=head1 DESCRIPTION

The segment object simplifies access to Bio::DB::SeqFeature store by
acting as a placeholder for a region of the genome. You can replace
this statement:

 @features = $db->features(-seq_id=>'Chr1',
                           -start=>5000,
                           -end=>6000,
                           -types=>['mRNA','match','repeat_region']);

with these statements:

 $segment = $db->segment('Chr1',5000=>6000);
 @features = $segment->features('mRNA','match','repeat_region');

You can also initialize a segment from an existing SeqFeature
object. The range will be picked up from the SeqFeature boundaries:

 $segment = Bio::DB::SeqFeature::Segment->new($feature);        # for Bio::DB::SeqFeature
 $segment = Bio::DB::SeqFeature::Segment->new($feature,$store); # for other Bio::SeqFeatureI objects

The segment object implements the full Bio::SeqFeature::CollectionI
interface, thereby allowing you to iterate over all features in the
range.

=cut

use strict;

use base 'Bio::SeqFeature::CollectionI','Bio::RangeI';
use Bio::DB::GFF::Util::Rearrange;
use overload '""' => \&as_string,
  fallback => 1;

=head1 PUBLIC METHODS

The following are public methods intended for external use.

=head2 new

 Title   : new
 Usage   : $segment = Bio::DB::SeqFeature::Segment->new(@options)
 Function: create a new Segment object
 Returns : A Bio::DB::SeqFeature::Segment object
 Args    : several - see below
 Status  : public

This class method creates a Bio::DB::SeqFeature::Segment object. You
must provide a Bio::DB::SeqFeature::Store as well as the coordinates
of the segment. These arguments can be provided explicitly or
indirectly.

First form:

 $segment = Bio::DB::SeqFeature::Segment->new($store,$seqid,$start,$end,$strand)

In this form a segment is defined by a Bio::DB::SeqFeature::Store, the
sequence ID, the start, end and strand. This is the form that is
invoked internally by Bio::DB::SeqFeature::Store when you call its
segment() method.

Second form:

 $segment = Bio::DB::SeqFeature::Segment->new($seqfeature [,$store]);

In this form, you pass new() a Bio::SeqFeatureI object. The segment is
constructed from the seq_id and coordinates are taken from the
object. If you pass a store-aware seqfeature object
(e.g. Bio::DB::SeqFeature) then the store database is also derived
from the feature. Otherwise you will have to pass the store as a
second argument.

=cut

###
# new()
#
# Call as Bio::DB::SeqFeature::Segment->new($seqfeature,$store)
#
# or
# Bio::DB::SeqFeature::Segment->new(-seqid=>$seqid,-start=>$start,-end=>$end,-strand=>$strand,-store=>$store)
#
sub new {
  my $class = shift;
  my ($store,$seqid,$start,$end,$strand,$id);
  if (ref $_[0] && UNIVERSAL::isa($_[0],'Bio::SeqFeatureI')) {
    my $seqfeature = shift;
    $store      = shift;
    $store       ||= eval {$seqfeature->object_store};
    $class->throw("I could not derive the Bio::DB::SeqFeature::Store object from the arguments passed to Bio::DB::SeqFeature::Segment->new(). Please pass the Store object as the second argument") unless $store;
    $seqid = $seqfeature->seq_id;
    $start = $seqfeature->start;
    $end   = $seqfeature->end;
    $strand= $seqfeature->strand;
    $id    = eval{$seqfeature->primary_id};
  }
  else {
    ($store,$seqid,$start,$end,$strand,$id) = @_;
  }
  return bless {
		store => $store,
		seqid => $seqid,
		start => $start,
		end   => $end,
		strand => $strand,
		primary_id => $id,
	       },ref($class) || $class;
}

=head2 features

 Title   : features
 Usage   : @features = $segment->features(@args)
 Function: fetch seqfeatures that overlap the segment
 Returns : list of features
 Args    : see below
 Status  : Public

This is the workhorse for feature query and retrieval. It takes a
series of -name=E<gt>$value arguments filter arguments. Features that
match all the filters are returned.

  Argument       Value
  --------       -----

 Location filters:
  -strand        Strand
  -range_type    Type of range match ('overlaps','contains','contained_in')

 Name filters:
  -name          Name of feature (may be a glob expression)
  -aliases       If true, match aliases as well as display names
  -class         Archaic argument for backward compatibility.
                  (-class=>'Clone',-name=>'ABC123') is equivalent
                  to (-name=>'Clone:ABC123')

 Type filters:
  -types         List of feature types (array reference) or one type (scalar)
  -type          Synonym for the above
  -primary_tag   Synonym for the above

  -attributes    Hashref of attribute=>value pairs as per
                    get_features_by_attribute(). Multiple alternative values
                    can be matched by providing an array reference.
  -attribute     synonym for -attributes

This is identical to the Bio::DB::SeqFeature::Store-E<gt>features()
method, except that the -seq_id, -start, and -end arguments are
provided by the segment object. If a simple list of arguments is
provided, then the list is taken to be the set of feature types
(primary tags) to filter on.

Examples:

All features that overlap the current segment:

 @features = $segment->features;

All features of type mRNA that overlap the current segment:

 @features = $segment->features('mRNA');

All features that are completely contained within the current segment:

 @features = $segment->features(-range_type=>'contains');

All "confirmed" mRNAs that overlap the current segment:

 @features = $segment->features(-attributes=>{confirmed=>1},-type=>'mRNA');

=cut

sub features {
  my $self = shift;
  my @args;
  if (@_ == 0) {
    @args = ();
  }
  elsif ($_[0] !~/^-/) {
    my @types = @_;
    @args = (-type=>\@types);
  } else {
    @args = @_;
  }
  $self->{store}->features(@args,-seqid=>$self->{seqid},-start=>$self->{start},-end=>$self->{end});
}

sub types {
    my $self = shift;
    my %types;
    my $iterator = $self->get_seq_stream(@_);
    while (my $f = $iterator->next_seq) {
	$types{$f->type}++;
    }
    return %types;
}

=head2 get_seq_stream

 Title   : get_seq_stream
 Usage   : $iterator = $segment->get_seq_stream(@args)
 Function: return an iterator across all features in the database
 Returns : a Bio::DB::SeqFeature::Store::Iterator object
 Args    : (optional) the feature() method
 Status  : public

This is identical to Bio::DB::SeqFeature::Store-E<gt>get_seq_stream()
except that the location filter is always automatically applied so
that the iterator you receive returns features that overlap the
segment's region.

When called without any arguments this method will return an iterator
object that will traverse all indexed features in the database that
overlap the segment's region. Call the iterator's next_seq() method to
step through them (in no particular order):

  my $iterator = $db->get_seq_stream;
  while (my $feature = $iterator->next_seq) {
    print $feature->primary_tag,' ',$feature->display_name,"\n";
  }

You can select a subset of features by passing a series of filter
arguments. The arguments are identical to those accepted by
$segment-E<gt>features().

get_feature_stream() ican be used as a synonym for this method.

=cut

#'

sub get_seq_stream {
  my $self = shift;
  $self->{store}->get_seq_stream(@_,-seqid=>$self->{seqid},-start=>$self->{start},-end=>$self->{end});
}

sub get_feature_stream { shift->get_seq_stream(@_) }

=head2 store

 Title   : store
 Usage   : $store = $segment->store
 Function: return the Bio::DB::SeqFeature::Store object associated with the segment
 Returns : a Bio::DB::SeqFeature::Store: object
 Args    : none
 Status  : public

=cut

sub factory { shift->{store} }
sub store   { shift->{store} }

=head2 primary_tag, type,

 Title   : primary_tag,type
 Usage   : $primary_tag = $segment->primary_tag
 Function: returns the string "region"
 Returns : "region"
 Args    : none
 Status  : public

The primary_tag method returns the constant tag "region". type() is a
synonym for this method.

=cut

sub type    { shift->primary_tag }

=head2 as_string

 Title   : as_string
 Usage   : $name = $segment->as_string
 Function: expands the object into a human-readable string
 Returns : "seq_id:start..end"
 Args    : none
 Status  : public

The as_string() method is overloaded into the "" operator so that the
object is represented as a human readable string in the form
"seq_id:start..end" when used in a string context.

=cut

sub as_string {
  my $self = shift;
  my $label = $self->seq_id;
  my $start = $self->start || '';
  my $end   = $self->end   || '';
  return "$label:$start..$end";
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

  my ($start,$strand) = ($self->start,$self->strand);
  @result = $strand < 0 ? map { $start - $_ + 1 } @_
                        : map { $_ + $start - 1 } @_;
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

  my ($start,$strand) = ($self->start,$self->abs_strand);
  @result = $strand < 0 ? map { $start - $_ + 1 } @_
                        : map { $_ - $start + 1 } @_;

  # if called with a single argument, caller will expect a single scalar reply
  # not the size of the returned array!
  return $result[0] if @result == 1 and !wantarray;
  @result;
}



=head2 Bio::SeqFeatureI compatibility methods

For convenience, segments are interchangeable with Bio::SeqFeature
objects in many cases. This means that segments can be passed to
BioPerl modules that expect Bio::SeqFeature objects and they should
work as expected. The primary tag of segment objects is "region"
(SO:0000001 "Continous sequence E<gt>=1 base pair").

All these methods are read-only except for the primary_id, which can
be get or set.

The following Bio::SeqFeatureI methods are supported:

=over 4

=item start

=item end

=item seq_id

=item strand

=item length

=item display_name

=item primary_id

=item primary_tag (always returns "region")

=item source_tag (always returns "Bio::DB::SeqFeature::Segment")

=item get_SeqFeatures (always returns an empty list)

=item seq

=item entire_seq

=item location

=item All Bio::RangeI methods

=back

=cut

sub start   { shift->{start}  }
sub end     { shift->{end}    }
sub seq_id  { shift->{seqid}  }
sub strand  { shift->{strand} }
sub ref     { shift->seq_id   }
*refseq = \&ref;

sub length  {
  my $self = shift;
  return abs($self->end - $self->start) +1;
}

sub primary_tag  { 'region' }
sub source_tag   { __PACKAGE__ }
sub display_name { shift->as_string }
sub name         { shift->display_name }
sub class        { 'region' }
sub abs_ref      { shift->ref}
sub abs_start    { shift->start}
sub abs_end      { shift->end}
sub abs_strand   { shift->strand}
sub get_SeqFeatures { }
sub get_all_tags { }
sub get_tag_values { }
sub add_tag_value { }
sub remove_tag { }
sub has_tag { }
sub seq {
  my $self = shift;
  require Bio::PrimarySeq unless Bio::PrimarySeq->can('new');
  my ($start,$end) = ($self->start,$self->end);
  if ($self->strand < 0) {
    ($start,$end) = ($end,$start);
  }
  return Bio::PrimarySeq->new(
			      -seq => $self->store->fetch_sequence($self->seq_id,$start,$end),
			      -id  => $self->display_name);
}
sub subseq {
  my $self = shift;
  my ($newstart,$newstop) = @_;
  my $store = $self->store or return;
  my $seq   = $store->fetch_sequence($self->seq_id,$self->start+$newstart-1,$self->end+$newstop-1);
  return Bio::PrimarySeq->new(-seq=>$seq);
}
sub dna {
  my $seq = shift->seq;
  $seq    = $seq->seq if CORE::ref($seq);
  return $seq;
}

sub entire_seq {
  my $self = shift;
  require Bio::PrimarySeq unless Bio::PrimarySeq->can('new');
  return Bio::PrimarySeq->new(
			      -seq => $self->store->fetch_sequence($self->seq_id),
			      -id  => $self->seq_id);
}

sub location {
  my $self = shift;
  require Bio::Location::Simple unless Bio::Location::Simple->can('new');  
  my $loc = Bio::Location::Simple->new(-start  => $self->start,
				       -end    => $self->end,
				       -strand => $self->strand);
  $loc->strand($self->strand);
  return $loc;
}
sub primary_id   {
  my $self = shift;
  my $d    = $self->{primary_id};
  $self->{primary_id} = shift if @_;
  $d;
}

sub target { return }
sub score  { return }
sub stop   { shift->end }
sub absolute { return 1 }
sub desc   { shift->as_string }
sub display_id { shift->display_name }
sub primary_seq { shift->seq }
sub accession_number { return undef }  # intended return undef
sub alphabet { return undef }          # intended return undef

1;

__END__

=head1 BUGS

This is an early version, so there are certainly some bugs. Please
use the BioPerl bug tracking system to report bugs.

=head1 SEE ALSO

L<bioperl>,
L<Bio::DB::SeqFeature::Store>,
L<Bio::DB::SeqFeature::GFF3Loader>,
L<Bio::DB::SeqFeature::Store::DBI::mysql>,
L<Bio::DB::SeqFeature::Store::bdb>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2006 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut


