# $Id$
#
# BioPerl module for Bio::Seq::SeqWithQuality
#
# Cared for by Chad Matsalla <bioinformatics@dieselwurks.com
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::SequenceTrace - Bioperl object packaging a sequence with its trace

=head1 SYNOPSIS

=head1 DESCRIPTION

This object stores a sequence with its trace.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org             - General discussion
  http://bio.perl.org/MailList.html - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Chad Matsalla

Email bioinformatics@dieselwurks.com

=head1 CONTRIBUTORS 

Jason Stajich, jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


package Bio::Seq::SequenceTrace;

use vars qw(@ISA);

use strict;
use Bio::Root::Root;
use Bio::Seq::QualI;
use Bio::PrimarySeqI;
use Bio::PrimarySeq;
use Bio::Seq::PrimaryQual;
use Bio::Seq::TraceI;

@ISA = qw(Bio::Root::Root Bio::Seq::SeqWithQuality Bio::Seq::TraceI);

=head2 new()

 Title   : new()
 Usage   : $st = Bio::Seq::SequenceTrace->new
     (    -sequencewithquality =>   Bio::Seq::SequenceWithQuality,
          -trace_a  =>   \@trace_values_for_a_channel,
          -trace_t  =>   \@trace_values_for_t_channel,
          -trace_g  =>   \@trace_values_for_g_channel,
          -trace_c  =>   \@trace_values_for_c_channel,
          -trace_indices    => '0 5 10 15 20 25 30 35'
     );
 Function: Returns a new Bio::Seq::SequenceTrace object from basic
        constructors.
 Returns : a new Bio::Seq::SequenceTrace object
Arguments: I think that these are all describes in the usage above.

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
	# default: turn OFF the warnings
	$self->{supress_warnings} = 1;
    my($sequence_with_quality,$trace_indices,$trace_a,$trace_t,
          $trace_g,$trace_c) =
          $self->_rearrange([qw(
               SEQUENCEWITHQUALITY
               TRACE_INDICES
               TRACE_A
               TRACE_T
               TRACE_G)], @args);
          # first, deal with the sequence and quality information
     if ($sequence_with_quality && ref($sequence_with_quality) eq "Bio::Seq::SeqWithQuality") {
          $self->{swq} = $sequence_with_quality;
     }
     else {
          $self->throw("A Bio::Seq::SequenceTrace object must be created with a
               Bio::Seq::SeqWithQuality object.");
     }
     $self->{trace_a} = $trace_a ? $trace_a : undef;
     $self->{trace_t} = $trace_t ? $trace_t : undef;
     $self->{trace_g} = $trace_g ? $trace_g : undef;
     $self->{trace_c} = $trace_c ? $trace_c : undef;
     $self->{trace_indices} = $trace_indices ? $trace_indices : undef;
    return $self;
}

=head2 trace($base,\@new_values)

 Title   : trace($base,\@new_values)
 Usage   : @trace_Values  = @{$obj->trace($base,\@new_values)};
 Function: Returns the trace values as a reference to an array containing the
     trace values. The individual elements of the trace array are not validated
     and can be any numeric value.
 Returns : A reference to an array.
 Status  : 
Arguments: $base : which color channel would you like the trace values for?
               - $base must be one of "A","T","G","C"
          \@new_values : a reference to an array of values containing trace
               data for this base

=cut

sub trace {
   my ($self,$base_channel,$values) = @_;
     $base_channel =~ tr/A-Z/a-z/;
     if (length($base_channel) > 1 && $base_channel !~ /a|t|g|c/) {
          $self->throw("The base channel must be a, t, g, or c");
     }
     if ( $values && ref($values) eq "ARRAY") {
          $self->{trace_$base_channel} = $values;
     }
     elsif ($values) {
          $self->warn("You tried to change the traces for the $base_channel but
               the values you wave were not a reference to an array.");
     }
     return $self->{trace_$base_channel};
}


=head2 trace_indices($new_indices)

 Title   : trace_indices($new_indices)
 Usage   : $indices = $obj->trace_indices($new_indices);
 Function: Return the trace iindex points for this object.
 Returns : A scalar
 Args    : If used, the trace indices will be set to the provided value.

=cut

sub trace_indices {
   my ($self,$trace_indices)= @_;
     if ($trace_indices) { $self->{trace_indices} = $trace_indices; }
     return $self->{trace_indices};
}













=head2 _common_id()

 Title   : _common_id()
 Usage   : $common_id = $self->_common_id();
 Function: Compare the display_id of {qual_ref} and {seq_ref}.
 Returns : Nothing if they don't match. If they do return
	   {seq_ref}->display_id()
 Args    : None.

=cut

#'
sub _common_id {
	my $self = shift;
	return if (!$self->{seq_ref} || !$self->{qual_ref});
	my $sid = $self->{seq_ref}->display_id();
	return if (!$sid);
	return if (!$self->{qual_ref}->display_id());
	return $sid if ($sid eq $self->{qual_ref}->display_id());
		# should this become a warning?
		# print("ids $sid and $self->{qual_ref}->display_id() do not match. Bummer.\n");
}

=head2 _common_display_id()

 Title   : _common_id()
 Usage   : $common_id = $self->_common_display_id();
 Function: Compare the display_id of {qual_ref} and {seq_ref}.
 Returns : Nothing if they don't match. If they do return
	   {seq_ref}->display_id()
 Args    : None.

=cut

#'
sub _common_display_id {
	my $self = shift;
	$self->common_id();
}

=head2 _common_accession_number()

 Title   : _common_accession_number()
 Usage   : $common_id = $self->_common_accession_number();
 Function: Compare the accession_number() of {qual_ref} and {seq_ref}.
 Returns : Nothing if they don't match. If they do return
	   {seq_ref}->accession_number()
 Args    : None.

=cut

#'
sub _common_accession_number {
	my $self = shift;
	return if ($self->{seq_ref} || $self->{qual_ref});
	my $acc = $self->{seq_ref}->accession_number();
		# if (!$acc) { print("the seqref has no acc.\n"); }
	return if (!$acc);
		# if ($acc eq $self->{qual_ref}->accession_number()) { print("$acc matches ".$self->{qual_ref}->accession_number()."\n"); }
	return $acc if ($acc eq $self->{qual_ref}->accession_number());
		# should this become a warning?
		# print("accession numbers $acc and $self->{qual_ref}->accession_number() do not match. Bummer.\n");
}

=head2 _common_primary_id()

 Title   : _common_primary_id()
 Usage   : $common_primard_id = $self->_common_primary_id();
 Function: Compare the primary_id of {qual_ref} and {seq_ref}.
 Returns : Nothing if they don't match. If they do return
	   {seq_ref}->primary_id()
 Args    : None.

=cut

#'
sub _common_primary_id {
	my $self = shift;
	return if ($self->{seq_ref} || $self->{qual_ref});
	my $pid = $self->{seq_ref}->primary_id();
	return if (!$pid);
	return $pid if ($pid eq $self->{qual_ref}->primary_id());
		# should this become a warning?
		# print("primary_ids $pid and $self->{qual_ref}->primary_id() do not match. Bummer.\n");

}

=head2 _common_desc()

 Title   : _common_desc()
 Usage   : $common_desc = $self->_common_desc();
 Function: Compare the desc of {qual_ref} and {seq_ref}.
 Returns : Nothing if they don't match. If they do return
	   {seq_ref}->desc()
 Args    : None.

=cut

#'
sub _common_desc {
	my $self = shift;
	return if ($self->{seq_ref} || $self->{qual_ref});
	my $des = $self->{seq_ref}->desc();
	return if (!$des);
	return $des if ($des eq $self->{qual_ref}->desc());
		# should this become a warning?
		# print("descriptions $des and $self->{qual_ref}->desc() do not match. Bummer.\n");

}

=head2 set_common_descriptors()

 Title   : set_common_descriptors()
 Usage   : $self->set_common_descriptors();
 Function: Compare the descriptors (id,accession_number,display_id,
	primary_id, desc) for the PrimarySeq and PrimaryQual objects
	within the SeqWithQuality object. If they match, make that
	descriptor the descriptor for the SeqWithQuality object.
 Returns : Nothing.
 Args    : None.

=cut

sub set_common_descriptors {
	my $self = shift;
	return if ($self->{seq_ref} || $self->{qual_ref});
	&_common_id();
	&_common_display_id();
	&_common_accession_number();
	&_common_primary_id();
	&_common_desc();
}

=head2 alphabet()

 Title   : alphabet();
 Usage   : $molecule_type = $obj->alphabet();
 Function: Get the molecule type from the PrimarySeq object.
 Returns : What what PrimarySeq says the type of the sequence is.
 Args    : None.

=cut

sub alphabet {
	my $self = shift;
	return $self->{seq_ref}->alphabet();	
}

=head2 display_id()

 Title   : display_id()
 Usage   : $id_string = $obj->display_id();
 Function: Returns the display id, aka the common name of the Quality
        object.
        The semantics of this is that it is the most likely string to be
        used as an identifier of the quality sequence, and likely to have
        "human" readability.  The id is equivalent to the ID field of the
        GenBank/EMBL databanks and the id field of the Swissprot/sptrembl
        database. In fasta format, the >(\S+) is presumed to be the id,
        though some people overload the id to embed other information.
        Bioperl does not use any embedded information in the ID field,
        and people are encouraged to use other mechanisms (accession
	field for example, or extending the sequence object) to solve
	this. Notice that $seq->id() maps to this function, mainly for
        legacy/convience issues.
	This method sets the display_id for the SeqWithQuality object.
 Returns : A string
 Args    : If a scalar is provided, it is set as the new display_id for
	the SeqWithQuality object.
 Status  : Virtual

=cut

sub display_id {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'display_id'} = $value;
    }
    return $obj->{'display_id'};

}

=head2 accession_number()

 Title   : accession_number()
 Usage   : $unique_biological_key = $obj->accession_number();
 Function: Returns the unique biological id for a sequence, commonly
        called the accession_number. For sequences from established
        databases, the implementors should try to use the correct
        accession number. Notice that primary_id() provides the unique id
        for the implemetation, allowing multiple objects to have the same
        accession number in a particular implementation. For sequences
        with no accession number, this method should return "unknown".
	This method sets the accession_number for the SeqWithQuality
	object. 
 Returns : A string (the value of accession_number)
 Args    : If a scalar is provided, it is set as the new accession_number
	for the SeqWithQuality object.
 Status  : Virtual


=cut

sub accession_number {
    my( $obj, $acc ) = @_;

    if (defined $acc) {
        $obj->{'accession_number'} = $acc;
    } else {
        $acc = $obj->{'accession_number'};
        $acc = 'unknown' unless defined $acc;
    }
    return $acc;
}

=head2 primary_id()

 Title   : primary_id()
 Usage   : $unique_implementation_key = $obj->primary_id();
 Function: Returns the unique id for this object in this implementation.
        This allows implementations to manage their own object ids in a
        way the implementaiton can control clients can expect one id to
        map to one object. For sequences with no accession number, this
        method should return a stringified memory location.
	This method sets the primary_id for the SeqWithQuality
	object.
 Returns : A string. (the value of primary_id)
 Args    : If a scalar is provided, it is set as the new primary_id for
	the SeqWithQuality object.

=cut

sub primary_id {
   my ($obj,$value) = @_;
   if ($value) {
      $obj->{'primary_id'} = $value;
    }
   return $obj->{'primary_id'};

}

=head2 desc()

 Title   : desc()
 Usage   : $qual->desc($newval); _or_ 
           $description = $qual->desc();
 Function: Get/set description text for this SeqWithQuality object.
 Returns : A string. (the value of desc)
 Args    : If a scalar is provided, it is set as the new desc for the
	SeqWithQuality object.

=cut

sub desc {
	# a mechanism to set the disc for the SeqWithQuality object.
	# probably will be used most often by set_common_features()
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'desc'} = $value;
    }
    return $obj->{'desc'};
}

=head2 id()

 Title   : id()
 Usage   : $id = $qual->id();
 Function: Return the ID of the quality. This should normally be (and
        actually is in the implementation provided here) just a synonym
        for display_id().
 Returns : A string. (the value of id)
 Args    : If a scalar is provided, it is set as the new id for the
	SeqWithQuality object.

=cut

sub id {
   my ($self,$value) = @_;
   if (!$self) { $self->throw("no value for self in $value"); }
   if( defined $value ) {
       return $self->display_id($value);
   }
   return $self->display_id();
}

=head2 seq

 Title   : seq()
 Usage   : $string    = $obj->seq(); _or_
	$obj->seq("atctatcatca");
 Function: Returns the sequence that is contained in the imbedded in the
	PrimarySeq object within the SeqWithQuality object
 Returns : A scalar (the seq() value for the imbedded PrimarySeq object.)
 Args    : If a scalar is provided, the SeqWithQuality object will
	attempt to set that as the sequence for the imbedded PrimarySeq
	object. Otherwise, the value of seq() for the PrimarySeq object
	is returned.
 Notes   : This is probably not a good idea because you then should call
	length() to make sure that the sequence and quality are of the
	same length. Even then, how can you make sure that this sequence
	belongs with that quality? I provided this to give you rope to
	hang yourself with. Tie it to a strong device and use a good
	knot.

=cut

sub seq {
	my ($self,$value) = @_;
	if( defined $value) {
		$self->{seq_ref}->seq($value);
		$self->length();
	}
	return $self->{seq_ref}->seq();
}

=head2 qual()

 Title   : qual()
 Usage   : @quality_values  = @{$obj->qual()}; _or_
	$obj->qual("10 10 20 40 50");
 Function: Returns the quality as imbedded in the PrimaryQual object
	within the SeqWithQuality object.
 Returns : A reference to an array containing the quality values in the 
	PrimaryQual object.
 Args    : If a scalar is provided, the SeqWithQuality object will
	attempt to set that as the quality for the imbedded PrimaryQual
	object. Otherwise, the value of qual() for the PrimaryQual
	object is returned.
 Notes   : This is probably not a good idea because you then should call
	length() to make sure that the sequence and quality are of the
	same length. Even then, how can you make sure that this sequence
	belongs with that quality? I provided this to give you a strong
	board with which to flagellate yourself. 

=cut

sub qual {
	my ($self,$value) = @_;

	if( defined $value) {
		$self->{qual_ref}->qual($value);
			# update the lengths
		$self->length();
	}
	return $self->{qual_ref}->qual();
}




=head2 length()

 Title   : length()
 Usage   : $length = $seqWqual->length();
 Function: Get the length of the SeqWithQuality sequence/quality.
 Returns : Returns the length of the sequence and quality if they are
	both the same. Returns "DIFFERENT" if they differ.
 Args    : None.

=cut

sub length {
    my $self = shift;
     # what do I return here? Whew. Ambiguity...
     ########

}


=head2 qual_obj

 Title   : qual_obj($different_obj)
 Usage   : $qualobj = $seqWqual->qual_obj(); _or_
	$qualobj = $seqWqual->qual_obj($ref_to_primaryqual_obj);
 Function: Get the PrimaryQual object that is imbedded in the
	SeqWithQuality object or if a reference to a PrimaryQual object
	is provided, set this as the PrimaryQual object imbedded in the
	SeqWithQuality object.
 Returns : A reference to a Bio::Seq::SeqWithQuality object.

=cut

sub qual_obj {
    my ($self,$value) = @_;
    return $self->{swq}->qual_obj($value);
}


=head2 seq_obj

 Title   : seq_obj()
 Usage   : $seqobj = $seqWqual->qual_obj(); _or_
	$seqobj = $seqWqual->seq_obj($ref_to_primary_seq_obj);
 Function: Get the PrimarySeq object that is imbedded in the
	SeqWithQuality object or if a reference to a PrimarySeq object is
	provided, set this as the PrimarySeq object imbedded in the
	SeqWithQuality object.
 Returns : A reference to a Bio::PrimarySeq object.

=cut

sub seq_obj {
    my ($self,$value) = @_;
    return $self->{swq}->seq_obj($value);
}

=head2 _set_descriptors

 Title   : _set_descriptors()
 Usage   : $seqWqual->_qual_obj($qual,$seq,$id,$acc,$pid,$desc,$given_id,
	$alphabet);
 Function: Set the descriptors for the SeqWithQuality object. Try to
	match the descriptors in the PrimarySeq object and in the
	PrimaryQual object if descriptors were not provided with
	construction.
 Returns : Nothing.
 Args    : $qual,$seq,$id,$acc,$pid,$desc,$given_id,$alphabet as found
	in the new() method.
 Notes   : Really only intended to be called by the new() method. If
	you want to invoke a similar function try
	set_common_descriptors().

=cut


sub _set_descriptors {
    my ($self,$qual,$seq,$id,$acc,$pid,$desc,$given_id,$alphabet) = @_;
     $self->{swq}->_seq_descriptors($qual,$seq,$id,$acc,$pid,$desc,$given_id,$alphabet);
}

=head2 subseq($start,$end)

 Title   : subseq($start,$end)
 Usage   : $subsequence = $obj->subseq($start,$end);
 Function: Returns the subseq from start to end, where the first base
           is 1 and the number is inclusive, ie 1-2 are the first two
           bases of the sequence.
 Returns : A string.
 Args    : Two positions.

=cut

sub subseq {
    my ($self,@args) = @_;
    # does a single value work?
    return $self->{swq}->subseq(@args);	
}

=head2 baseat($position)

 Title   : baseat($position)
 Usage   : $base_at_position_6 = $obj->baseat("6");
 Function: Returns a single base at the given position, where the first
	base is 1 and the number is inclusive, ie 1-2 are the first two
	bases of the sequence.
 Returns : A scalar.
 Args    : A position.

=cut

sub baseat {
    my ($self,$val) = @_;
    return $self->{swq}->subseq($val,$val);
}

=head2 subqual($start,$end)

 Title   : subqual($start,$end)
 Usage   : @qualities = @{$obj->subqual(10,20);
 Function: returns the quality values from $start to $end, where the
        first value is 1 and the number is inclusive, ie 1-2 are the
	first two bases of the sequence. Start cannot be larger than
	end but can be equal.
 Returns : A reference to an array.
 Args    : a start position and an end position

=cut

sub subqual {
    my ($self,@args) = @_;
    return $self->{swq}->subqual(@args);
}

=head2 qualat($position)

 Title   : qualat($position)
 Usage   : $quality = $obj->qualat(10);
 Function: Return the quality value at the given location, where the
        first value is 1 and the number is inclusive, ie 1-2 are the
	first two bases of the sequence. Start cannot be larger than
	end but can be equal.
 Returns : A scalar.
 Args    : A position.

=cut

sub qualat {
    my ($self,$val) = @_;
    return $self->{swq}->qualat($val);
}

=head2 sub_trace_index($start,$end)

 Title   : sub_trace_index($start,$end)
 Usage   : @trace_indices = @{$obj->sub_trace_index(10,20);
 Function: returns the trace index values from $start to $end, where the
        first value is 1 and the number is inclusive, ie 1-2 are the
	first two bases of the sequence. Start cannot be larger than
	end but can be e_trace_index.
 Returns : A reference to an array.
 Args    : a start position and an end position

=cut

sub sub_trace_index {
   my ($self,$start,$end) = @_;

   if( $start > $end ){
       $self->throw("in sub_trace_index, start [$start] has to be greater than end [$end]");
   }

   if( $start <= 0 || $end > $self->length ) {
       $self->throw("You have to have start positive and length less than the total length of sequence [$start:$end] Total ".$self->length."");
   }

   # remove one from start, and then length is end-start

   $start--;
     $end--;
     my @sub_trace_index_array = @{$self->{trace_indices}}[$start..$end];

     #   return substr $self->seq(), $start, ($end-$start);
     return \@sub_trace_index_array;

}




=head2 trace_index_at($position)

 Title   : trace_index_at($position)
 Usage   : $trace_index = $obj->trace_index_at(10);
 Function: Return the trace_index value at the given location, where the
        first value is 1 and the number is inclusive, ie 1-2 are the
	first two bases of the sequence. Start cannot be larger than
	end but can be etrace_index_.
 Returns : A scalar.
 Args    : A position.

=cut

sub trace_index_at {
    my ($self,$val) = @_;
    my @trace_index_at = @{$self->sub_trace_index($val,$val)};
    if (scalar(@trace_index_at) == 1) {
	return $trace_index_at[0];
    }
    else {
	$self->throw("AAAH! trace_index_at provided more then one quality.");
    }
}


1;
