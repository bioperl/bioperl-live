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
use Bio::Annotation::Collection;

use base qw(Bio::Root::Root Bio::RangeI Bio::SeqI Bio::Das::SegmentI);

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
  $segclass = $segname->class if ref($segname) && $segname->can('class');
  $segclass ||= 'Sequence';

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

=head2 end

 Title   : end
 Usage   : $s->end
 Function: end of segment
 Returns : integer
 Args    : none
 Status  : Public

This is a read-only accessor for the end of the segment.

=cut

sub end   { shift->{stop}  }

=head2 stop

 Title   : stop
 Usage   : $s->stop
 Function: stop of segment
 Returns : integer
 Args    : none
 Status  : Public

This is an alias for end(), provided for AcePerl compatibility.

=cut

*stop = \&end;

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
 Returns : +1,0,-1
 Args    : none
 Status  : Public

Returns the strand on which the segment resides, either +1, 0 or -1.

=cut

sub strand {
  my $self = shift;
  0;
}

=head2 low

 Title   : low
 Usage   : $s->low
 Function: return lower coordinate
 Returns : lower coordinate
 Args    : none
 Status  : Public

Returns the lower coordinate, either start or end.

=cut

sub low {
  my $self = shift;
  my ($start,$stop) = ($self->start,$self->stop);
  return $start < $stop ? $start : $stop;
}
*abs_low = \&low;

=head2 high

 Title   : high
 Usage   : $s->high
 Function: return higher coordinate
 Returns : higher coordinate
 Args    : none
 Status  : Public

Returns the higher coordinate, either start or end.

=cut

sub high {
  my $self = shift;
  my ($start,$stop) = ($self->start,$self->stop);
  return $start > $stop ? $start : $stop;
}
*abs_high = \&high;

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
 Usage   : $s->class([$newclass])
 Function: get the source sequence class
 Returns : a string
 Args    : new class (optional)
 Status  : Public

Gets or sets the class for the source sequence for this segment.

=cut

sub class     { 
  my $self = shift;
  my $d = $self->{class};
  $self->{class} = shift if @_;
  $d;
}

=head2 subseq

 Title   : subseq
 Usage   : $s->subseq($start,$stop)
 Function: generate a subsequence
 Returns : a Bio::DB::GFF::Segment object
 Args    : start and end of subsequence
 Status  : Public

This method generates a new segment from the start and end positions
given in the arguments.  If stop E<lt> start, then the strand is reversed.

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

=head2 seq

 Title   : seq
 Usage   : $s->seq
 Function: get the sequence string for this segment
 Returns : a Bio::PrimarySeq
 Args    : none
 Status  : Public

Returns the sequence for this segment as a Bio::PrimarySeq.  (-)
strand segments are automatically reverse complemented

The method is called dna() return the data as a simple sequence
string.

=cut

sub seq {
  my $self = shift;
  my $dna = $self->dna;
  require Bio::PrimarySeq unless Bio::PrimarySeq->can('new');
  return Bio::PrimarySeq->new(-id => $self->display_name) unless $dna;
  return Bio::PrimarySeq->new(-seq => $dna,
			      -id  => $self->display_name);
}

=head2 dna

 Title   : dna
 Usage   : $s->dna
 Function: get the DNA string for this segment
 Returns : a string
 Args    : none
 Status  : Public

Returns the sequence for this segment as a simple string. (-) strand
segments are automatically reverse complemented

The method is also called protein().

=cut

sub dna {
  my $self = shift;
  my ($ref,$class,$start,$stop,$strand) 
    = @{$self}{qw(sourceseq class start stop strand)};
  return $self->factory->dna($ref,$start,$stop,$class);
}

*protein = \&dna;


=head2 primary_seq

 Title   : primary_seq
 Usage   : $s->primary_seq
 Function: returns a Bio::PrimarySeqI compatible object
 Returns : a Bio::PrimarySeqI object
 Args    : none
 Status  : Public

This is for compatibility with BioPerl's separation of SeqI
from PrimarySeqI.  It just returns itself.

=cut

#'

sub primary_seq { shift }

=head2 type

 Title   : type
 Usage   : $s->type
 Function: return the string "feature"
 Returns : the string "feature"
 Args    : none
 Status  : Public

This is for future sequence ontology-compatibility and
represents the default type of a feature on the genome

=cut

sub type { "feature" }

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
  return unless defined $peer;
  return $self->asString eq $peer unless ref($peer) && $peer->isa('Bio::DB::GFF::Segment');
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

=head2 abs_end

 Title   : abs_end
 Usage   : $s->abs_end
 Function: the absolute stop of the segment
 Returns : an integer
 Args    : none
 Status  : Public

This is an alias to stop(), and provided for API compatibility with
Bio::DB::GFF::RelSegment.

=cut

*abs_stop   = \&stop;
*abs_end    = \&stop;

=head2 abs_strand

 Title   : abs_strand
 Usage   : $s->abs_strand
 Function: the absolute strand of the segment
 Returns : +1,0,-1
 Args    : none
 Status  : Public

This is an alias to strand(), and provided for API compatibility with
Bio::DB::GFF::RelSegment.

=cut

sub abs_strand {
  my $self = shift;
  return $self->abs_end <=> $self->abs_start;
}

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

=head2 seq_id

 Title   : seq_id
 Usage   : $ref = $s->seq_id
 Function: get the reference sequence in a LocationI-compatible way
 Returns : a string
 Args    : none
 Status  : Public

An alias for refseq() but only allows reading.

=cut

sub seq_id { shift->refseq }
*seqname = \&seq_id;

=head2 truncated

 Title   : truncated
 Usage   : $truncated = $s->truncated
 Function: Flag indicating that the segment was truncated during creation
 Returns : A boolean flag
 Args    : none
 Status  : Public

This indicates that the sequence was truncated during creation.  The
returned flag is undef if no truncation occured.  If truncation did
occur, the flag is actually an array ref in which the first element is
true if truncation occurred on the left, and the second element
occurred if truncation occurred on the right.

=cut

sub truncated {
  my $self = shift;
  my $hash = $self->{truncated} or return;
  CORE::ref($hash) eq 'HASH' or return [1,1];  # paranoia -- not that this would ever happen ;-)
  return [$hash->{start},$hash->{stop}];
}

=head2 Bio::RangeI Methods

The following Bio::RangeI methods are supported:

overlaps(), contains(), equals(),intersection(),union(),overlap_extent()

=cut

sub overlaps {
  my $self  = shift;
  my($other,$so) = @_;
  if ($other->isa('Bio::DB::GFF::RelSegment')) {
    return if $self->abs_ref ne $other->abs_ref;
  }
  $self->SUPER::overlaps(@_);
}

sub contains {
  my $self  = shift;
  my($other,$so) = @_;
  if ($other->isa('Bio::DB::GFF::RelSegment')) {
    return if $self->abs_ref ne $other->abs_ref;
  }
  $self->SUPER::contains(@_);
}
#sub equals {
#  my $self  = shift;
#  my($other,$so) = @_;
#  if ($other->isa('Bio::DB::GFF::RelSegment')) {
#    return if $self->abs_ref ne $other->abs_ref;
#  }
#  $self->SUPER::equals(@_);
#}
sub intersection {
  my $self  = shift;
  my($other,$so) = @_;
  if ($other->isa('Bio::DB::GFF::RelSegment')) {
    return if $self->abs_ref ne $other->abs_ref;
  }
  $self->SUPER::intersection(@_);
}
sub union {
  my $self  = shift;
  my($other) = @_;
  if ($other->isa('Bio::DB::GFF::RelSegment')) {
    return if $self->abs_ref ne $other->abs_ref;
  }
  $self->SUPER::union(@_);
}

sub overlap_extent {
  my $self  = shift;
  my($other) = @_;
  if ($other->isa('Bio::DB::GFF::RelSegment')) {
    return if $self->abs_ref ne $other->abs_ref;
  }
  $self->SUPER::overlap_extent(@_);
}


=head2 Bio::SeqI implementation

=cut

=head2 primary_id

 Title   : primary_id
 Usage   : $unique_implementation_key = $obj->primary_id;
 Function: Returns the unique id for this object in this
           implementation. This allows implementations to manage their
           own object ids in a way the implementaiton can control
           clients can expect one id to map to one object.

           For sequences with no accession number, this method should
           return a stringified memory location.

 Returns : A string
 Args    : None
 Status  : Virtual


=cut

sub primary_id {
   my ($obj,$value) = @_;

   if( defined $value) {
      $obj->{'primary_id'} = $value;
    }
   if( ! exists $obj->{'primary_id'} ) {
       return "$obj";
   }
   return $obj->{'primary_id'};
}


=head2 display_name

 Title   : display_name
 Usage   : $id = $obj->display_name or $obj->display_name($newid);
 Function: Gets or sets the display id, also known as the common name of
           the Seq object.

           The semantics of this is that it is the most likely string
           to be used as an identifier of the sequence, and likely to
           have "human" readability.  The id is equivalent to the LOCUS
           field of the GenBank/EMBL databanks and the ID field of the
           Swissprot/sptrembl database. In fasta format, the >(\S+) is
           presumed to be the id, though some people overload the id
           to embed other information. Bioperl does not use any
           embedded information in the ID field, and people are
           encouraged to use other mechanisms (accession field for
           example, or extending the sequence object) to solve this.

           Notice that $seq->id() maps to this function, mainly for
           legacy/convenience issues.
 Returns : A string
 Args    : None or a new id

Note, this used to be called display_id(), and this name is preserved for
backward compatibility.  The default is to return the seq_id().

=cut

sub display_name { shift->seq_id }
*display_id = \&display_name;

=head2 accession_number

 Title   : accession_number
 Usage   : $unique_biological_key = $obj->accession_number;
 Function: Returns the unique biological id for a sequence, commonly
           called the accession_number. For sequences from established
           databases, the implementors should try to use the correct
           accession number. Notice that primary_id() provides the
           unique id for the implemetation, allowing multiple objects
           to have the same accession number in a particular implementation.

           For sequences with no accession number, this method should return
           "unknown".
 Returns : A string
 Args    : None


=cut

sub accession_number {
    return 'unknown';
}

=head2 alphabet

 Title   : alphabet
 Usage   : if( $obj->alphabet eq 'dna' ) { /Do Something/ }
 Function: Returns the type of sequence being one of
           'dna', 'rna' or 'protein'. This is case sensitive.

           This is not called <type> because this would cause
           upgrade problems from the 0.5 and earlier Seq objects.

 Returns : a string either 'dna','rna','protein'. NB - the object must
           make a call of the type - if there is no type specified it
           has to guess.
 Args    : none
 Status  : Virtual


=cut

sub alphabet{
    return 'dna'; # no way this will be anything other than dna!
}

=head2 desc

 Title   : desc
 Usage   : $seqobj->desc($string) or $seqobj->desc()
 Function: Sets or gets the description of the sequence
 Example :
 Returns : The description
 Args    : The description or none


=cut

sub desc { shift->asString }

*description = \&desc;

=head2 species

 Title   : species
 Usage   : $species = $seq->species() or $seq->species($species)
 Function: Gets or sets the species
 Example :
 Returns : Bio::Species object
 Args    : None or Bio::Species object

See L<Bio::Species> for more information

=cut

sub species {
    my ($self, $species) = @_;
    if ($species) {
        $self->{'species'} = $species;
    } else {
        return $self->{'species'};
    }
}

=head2 annotation

 Title   : annotation
 Usage   : $ann = $seq->annotation or $seq->annotation($annotation)
 Function: Gets or sets the annotation
 Example :
 Returns : Bio::Annotation object
 Args    : None or Bio::Annotation object

See L<Bio::Annotation> for more information

=cut

sub annotation {
   my ($obj,$value) = @_;
   if( defined $value || ! defined $obj->{'annotation'} ) {
       $value = Bio::Annotation::Collection->new() unless defined $value;
      $obj->{'annotation'} = $value;
    }
    return $obj->{'annotation'};

}

=head2 is_circular

 Title   : is_circular
 Usage   : if( $obj->is_circular) { /Do Something/ }
 Function: Returns true if the molecule is circular
 Returns : Boolean value
 Args    : none

=cut

sub is_circular{
    return 0;
}


1;
__END__

=head1 BUGS

Report them please.

=head1 SEE ALSO

L<bioperl>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.  

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 CONTRIBUTORS

Jason Stajich E<lt>jason@bioperl.orgE<gt>.

=cut

