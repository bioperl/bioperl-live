=head1 NAME

Bio::DB::GFF::Segment -- Simple DNA segment object

=head1 SYNOPSIS

See L<Bio::DB::GFF>.

=head1 DESCRIPTION

Bio::DB::GFF::Segment provides the basic representation of a range of
DNA contained in a GFF database.  It is the base class from which the
Bio::DB::GFF::RelSegment and Bio::DB::GFF::Feature classes are
derived.

Generally, you will not create or manipulate Bio::DB::GFF::Segment
objects directly, but use those that are returned by the Bio::DB::GFF
module.

=cut

package Bio::DB::GFF::Segment;

use strict;
use Bio::Root::RootI;

use vars qw($VERSION @ISA);
@ISA = qw(Bio::Root::RootI);
$VERSION = '0.25';

use overload 
  '""'     => 'asString',
  eq       => 'equals',
  fallback => 1;

=head1 API

The remainder of this document describes the API for
Bio::DB::GFF::Segment.

=cut

=head2 new

 Title   : new
 Usage   : $s = Bio::DB::GFF::Segment->new(@args)
 Function: create a new segment
 Returns : a new Bio::DB::GFF::Segment object
 Args    : see below
 Status  : Public

This method creates a new Bio::DB::GFF::Segment object.  Generally
this is called automatically by the Bio::DB::GFF module and
derivatives.

There are five positional arguments:

 $factory      a Bio::DB::GFF::Adaptor to use for database access
 $sourceseq    ID of the source sequence
 $sourceclass  class of the source sequence
 $start        start of the desired segment relative to source sequence
 $stop         stop of the desired segment relative to source sequence

=cut

sub new {
  my $class = shift;
  my ($factory,$segclass,$segname,$start,$stop) = @_;

  $factory or $class->throw("->new(): provide a factory argument");
  $class = ref $class if ref $class;
  return bless { factory   => $factory,
		 sourceseq => $segname,
		 class     => $segclass,
		 start     => $start,
		 stop      => $stop,
		 strand    => 0,
	       },$class;
}

# read-only accessors

=head2 factory

 Title   : factory
 Usage   : $s->factory
 Function: get the factory object
 Returns : a Bio::DB::GFF::Adaptor
 Args    : none
 Status  : Public

This is a read-only accessor for the Bio::DB::GFF::Adaptor object used 
to create the segment.

=cut

sub factory { shift->{factory} }

# start, stop, length
=head2 start

 Title   : start
 Usage   : $s->start
 Function: start of segment
 Returns : integer
 Args    : none
 Status  : Public

This is a read-only accessor for the start of the segment.

=cut

sub start  { shift->{start} }

=head2 stop

 Title   : stop
 Usage   : $s->stop
 Function: stop of segment
 Returns : integer
 Args    : none
 Status  : Public

This is a read-only accessor for the stop of the segment.

=cut

sub stop   { shift->{stop}  }

=head2 end

 Title   : end
 Usage   : $s->end
 Function: stop of segment
 Returns : integer
 Args    : none
 Status  : Public

This is an alias for stop(), provided for BioSeq compatibility.

=cut

*end = \&stop;

=head2 length

 Title   : length
 Usage   : $s->length
 Function: length of segment
 Returns : integer
 Args    : none
 Status  : Public

Returns the length of the segment.  Always a positive number.

=cut

sub length { abs($_[0]->{start} - $_[0]->{stop})+1 }


=head2 strand

 Title   : strand
 Usage   : $s->strand
 Function: strand of segment
 Returns : +1,-1
 Args    : none
 Status  : Public

Returns the strand on which the segment resides, either +1 or -1.

=cut

sub strand {
  my $self = shift; 
  return $self->{stop} >= $self->{start} ? +1 : -1;
}

=head2 sourceseq

 Title   : sourceseq
 Usage   : $s->sourceseq
 Function: get the segment source
 Returns : a string
 Args    : none
 Status  : Public

Returns the name of the source sequence for this segment.

=cut

sub sourceseq { shift->{sourceseq} }

=head2 class

 Title   : class
 Usage   : $s->class
 Function: get the source sequence class
 Returns : a string
 Args    : none
 Status  : Public

Returns the class for the source sequence for this segment.

=cut

sub class     { shift->{class}     }

=head2 subseq

 Title   : subseq
 Usage   : $s->subseq($start,$stop)
 Function: generate a subsequence
 Returns : a Bio::DB::GFF::Segment object
 Args    : start and end of subsequence
 Status  : Public

This method generates a new segment from the start and end positions
given in the arguments.  If stop < start, then the strand is reversed.

=cut

sub subseq {
  my $self = shift;
  my ($newstart,$newstop) = @_;
  my ($refseq,$start,$stop,$class) = ($self->{sourceseq},
				      $self->{start},$self->{stop},
				      $self->class);

  # We deliberately force subseq to return objects of type RelSegment
  # Otherwise, when we get a subsequence from a Feature object,
  # its method and source go along for the ride, which is incorrect.
  my $new = $self->new_from_segment($self);
  if ($start <= $stop) {
    @{$new}{qw(start stop)} = ($start + $newstart - 1, $start + $newstop  - 1);
  } else {
    @{$new}{qw(start stop)} = ($start - ($newstart - 1), $start - ($newstop  - 1)),

  }

  $new;
}

=head2 dna

 Title   : dna
 Usage   : $s->dna
 Function: get the DNA string for this segment
 Returns : a string
 Args    : none
 Status  : Public

Returns the DNA for this segment as a simple string.  (-) strand
segments are automatically reverse complemented

=cut

sub dna {
  my $self = shift;
  my ($ref,$class,$start,$stop,$strand) 
    = @{$self}{qw(sourceseq class start stop strand)};
  ($start,$stop) = ($stop,$start) if $strand eq '-';
  $self->factory->dna($ref,$class,$start,$stop);
}

=head2 equals

 Title   : equals
 Usage   : $s->equals($d)
 Function: segment equality
 Returns : true, if two segments are equal
 Args    : another segment
 Status  : Public

Returns true if the two segments have the same source sequence, start and stop.

=cut

sub equals {
  my $self = shift;
  my $peer = shift;
  return $self->{start} eq $peer->{start}
         && $self->{stop}  eq $peer->{stop}
         && $self->{sourceseq} eq $peer->{sourceseq};
}

=head2 asString

 Title   : asString
 Usage   : $s->asString
 Function: human-readable string for segment
 Returns : a string
 Args    : none
 Status  : Public

Returns a human-readable string representing this sequence.  Format
is:
    
   sourceseq/start,stop

=cut

sub asString {
  my $self = shift;
  my $label = $self->refseq;
  my $start = $self->start;
  my $stop  = $self->stop;
  return "$label:$start,$stop";
}

=head2 clone

 Title   : clone
 Usage   : $copy = $s->clone
 Function: make a copy of this segment
 Returns : a Bio::DB::GFF::Segment object
 Args    : none
 Status  : Public

This method creates a copy of the segment and returns it.

=cut

# deep copy of the thing
sub clone {
  my $self = shift;
  my %h = %$self;
  return bless \%h,ref($self);
}

=head2 error

 Title   : error
 Usage   : $error = $s->error([$new_error])
 Function: get or set the last error
 Returns : a string
 Args    : an error message (optional)
 Status  : Public

In case of a fault, this method can be used to obtain the last error
message.  Internally it is called to set the error message.

=cut

sub error {
  my $self = shift;
  my $g = $self->{error};
  $self->{error} = shift if @_;
  $g;
}

=head1 Relative Addressing Methods

The following methods are provided for compatibility with
Bio::DB::GFF::RelSegment, which provides relative addressing
functions.

=head2 abs_start

 Title   : abs_start
 Usage   : $s->abs_start
 Function: the absolute start of the segment
 Returns : an integer
 Args    : none
 Status  : Public

This is an alias to start(), and provided for API compatibility with
Bio::DB::GFF::RelSegment.

=cut

*abs_start  = \&start;


=head2 abs_stop

 Title   : abs_stop
 Usage   : $s->abs_stop
 Function: the absolute stop of the segment
 Returns : an integer
 Args    : none
 Status  : Public

This is an alias to stop(), and provided for API compatibility with
Bio::DB::GFF::RelSegment.

=cut

*abs_stop   = \&stop;


=head2 abs_strand

 Title   : abs_strand
 Usage   : $s->abs_strand
 Function: the absolute strand of the segment
 Returns : +1,-1
 Args    : none
 Status  : Public

This is an alias to strand(), and provided for API compatibility with
Bio::DB::GFF::RelSegment.

=cut

*abs_strand = \&strand;


=head2 abs_ref

 Title   : abs_ref
 Usage   : $s->abs_ref
 Function: the reference sequence for this segment
 Returns : a string
 Args    : none
 Status  : Public

This is an alias to sourceseq(), and is here to provide API
compatibility with Bio::DB::GFF::RelSegment.

=cut

*abs_ref    = \&sourceseq;

=head2 refseq

 Title   : refseq
 Usage   : $s->refseq
 Function: get or set the reference sequence
 Returns : a string
 Args    : none
 Status  : Public

Examine or change the reference sequence. This is an alias to
sourceseq(), provided here for API compatibility with
Bio::DB::GFF::RelSegment.

=cut

*refseq     = \&sourceseq;

=head2 ref

 Title   : ref
 Usage   : $s->refseq
 Function: get or set the reference sequence
 Returns : a string
 Args    : none
 Status  : Public

An alias for refseq()

=cut

sub ref { shift->refseq(@_) }

1;
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

