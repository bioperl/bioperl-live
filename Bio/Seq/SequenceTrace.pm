#
# BioPerl module for Bio::Seq::SequenceTrace
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

  # example code here

=head1 DESCRIPTION

This object stores a sequence with its trace.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

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

=head1 AUTHOR - Chad Matsalla

Email bioinformatics@dieselwurks.com


The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


package Bio::Seq::SequenceTrace;


use strict;
use Bio::Seq::QualI;
use Bio::PrimarySeqI;
use Bio::PrimarySeq;
use Bio::Seq::PrimaryQual;

use base qw(Bio::Root::Root Bio::Seq::Quality Bio::Seq::TraceI);

=head2 new()

 Title   : new()
 Usage   : $st = Bio::Seq::SequenceTrace->new
     (    -swq =>   Bio::Seq::SequenceWithQuality,
          -trace_a  =>   \@trace_values_for_a_channel,
          -trace_t  =>   \@trace_values_for_t_channel,
          -trace_g  =>   \@trace_values_for_g_channel,
          -trace_c  =>   \@trace_values_for_c_channel,
          -accuracy_a    =>   \@a_accuracies,
          -accuracy_t    =>   \@t_accuracies,
          -accuracy_g    =>   \@g_accuracies,
          -accuracy_c    =>   \@c_accuracies,
          -peak_indices    => '0 5 10 15 20 25 30 35'
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
    my($swq,$peak_indices,$trace_a,$trace_t,
          $trace_g,$trace_c,$acc_a,$acc_t,$acc_g,$acc_c) =
          $self->_rearrange([qw(
               SWQ
               PEAK_INDICES
               TRACE_A
               TRACE_T
               TRACE_G
               TRACE_C
               ACCURACY_A
               ACCURACY_T
               ACCURACY_G
               ACCURACY_C )], @args);
          # first, deal with the sequence and quality information
     if ($swq && ref($swq) eq "Bio::Seq::Quality") {
          $self->{swq} = $swq;
     }
     else {
          $self->throw("A Bio::Seq::SequenceTrace object must be created with a
               Bio::Seq::Quality object. You provided this type of object: "
               .ref($swq));
     }
     if (!$acc_a) {
          # this means that you probably did not provide traces and accuracies
          # and that they need to be synthesized
          $self->set_accuracies();
     }
     else {
          $self->accuracies('a',$acc_a);
          $self->accuracies('t',$acc_t);
          $self->accuracies('g',$acc_g);
          $self->accuracies('c',$acc_c);
     }
     if (!$trace_a) {
          $self->_synthesize_traces();
     }
     else {
          $self->trace('a',$trace_a);
          $self->trace('t',$trace_t);
          $self->trace('g',$trace_g);
          $self->trace('c',$trace_c);
          $self->peak_indices($peak_indices);
     }
     $self->id($self->seq_obj->id);
    return $self;
}

sub swq_obj {
     my $self = shift;
     $self->warn('swq_obj() is deprecated: use seq_obj()');
     return $self->{swq};
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
     if (!$base_channel) {
          $self->throw('You must provide a valid base channel (atgc) to use trace()');
     }
     $base_channel =~ tr/A-Z/a-z/;
     if ($base_channel !~ /[acgt]/) {
          $self->throw('You must provide a valid base channel (atgc) to use trace()');
     }
     if ($values) {
             if (ref($values) eq "ARRAY") {
                  $self->{trace}->{$base_channel} = $values;
             }
             else {
                    my @trace = split(' ',$values);
                  $self->{trace}->{$base_channel} = \@trace;
             }
     }
     if ($self->{trace}->{$base_channel}) {
          return $self->{trace}->{$base_channel};
     }
     else {
          return;
     }
}


=head2 peak_indices($new_indices)

 Title   : peak_indices($new_indices)
 Usage   : $indices = $obj->peak_indices($new_indices);
 Function: Return the trace index points for this object.
 Returns : A scalar
 Args    : If used, the trace indices will be set to the provided value.

=cut

sub peak_indices {
   my ($self,$peak_indices)= @_;
     if ($peak_indices) {
          if (ref($peak_indices) eq "ARRAY") {
               $self->{peak_indices} = $peak_indices;
          }
          else {
               my @indices = split(' ',$peak_indices);
               $self->{peak_indices} = \@indices;
         } 
     }
     if (!$self->{peak_indices}) {
          my @temp = ();
          $self->{peak_indices} = \@temp;
     }
     return $self->{peak_indices};
}


=head2 _reset_peak_indices()

 Title   : _rest_peak_indices()
 Usage   : $obj->_reset_peak_indices();
 Function: Reset the peak indices.
 Returns : Nothing.
 Args    : None.
 Notes   : When you create a sub_trace_object, the peak indices
     will still be pointing to the apporpriate location _in the
     original trace_. In order to fix this, the initial value must
     be subtracted from each value here. ie. The first peak index
     must be "1".

=cut

sub _reset_peak_indices {
   my $self = shift;
     my $length = $self->length();
     my $subtractive = $self->peak_index_at(1);
     my ($original,$new);
     $self->peak_index_at(1,"null");
     for (my $counter=2; $counter<= $length; $counter++) {
          my $original = $self->peak_index_at($counter);
          $new = $original - $subtractive;
          $self->peak_index_at($counter,$new);
     }
     return;
}





=head2 peak_index_at($position)

 Title   : peak_index_at($position)
 Usage   : $peak_index = $obj->peak_index_at($postition);
 Function: Return the trace iindex point at this position
 Returns : A scalar
 Args    : If used, the trace index at this position will be 
     set to the provided value.

=cut

sub peak_index_at {
   my ($self,$position,$value)= @_;
   if ($value) {
          if ($value eq "null") {
               $self->peak_indices->[$position-1] = "0";
          }
          else {
               $self->peak_indices->[$position-1] = $value;
          }
   }
    return $self->peak_indices()->[$position-1];
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
	return $self->{swq}->alphabet;
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
	This method sets the display_id for the Quality object.
 Returns : A string
 Args    : If a scalar is provided, it is set as the new display_id for
	the Quality object.
 Status  : Virtual

=cut

sub display_id {
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{swq}->display_id($value);
    }
    return $self->{swq}->display_id();

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
	This method sets the accession_number for the Quality
	object. 
 Returns : A string (the value of accession_number)
 Args    : If a scalar is provided, it is set as the new accession_number
	for the Quality object.
 Status  : Virtual


=cut

sub accession_number {
    my( $self, $acc ) = @_;
    if (defined $acc) {
        $self->{swq}->accession_number($acc);
    } else {
        $acc = $self->{swq}->accession_number();
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
	This method sets the primary_id for the Quality
	object.
 Returns : A string. (the value of primary_id)
 Args    : If a scalar is provided, it is set as the new primary_id for
	the Quality object.

=cut

sub primary_id {
   my ($self,$value) = @_;
   if ($value) {
      $self->{swq}->primary_id($value);
    }
   return $self->{swq}->primary_id();

}

=head2 desc()

 Title   : desc()
 Usage   : $qual->desc($newval); _or_ 
           $description = $qual->desc();
 Function: Get/set description text for this Quality object.
 Returns : A string. (the value of desc)
 Args    : If a scalar is provided, it is set as the new desc for the
	   Quality object.

=cut

sub desc {
	# a mechanism to set the desc for the Quality object.
	# probably will be used most often by set_common_features()
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{swq}->desc($value);
    }
      return $self->{swq}->desc();
}

=head2 id()

 Title   : id()
 Usage   : $id = $qual->id();
 Function: Return the ID of the quality. This should normally be (and
        actually is in the implementation provided here) just a synonym
        for display_id().
 Returns : A string. (the value of id)
 Args    : If a scalar is provided, it is set as the new id for the
	   Quality object.

=cut

sub id {
   my ($self,$value) = @_;
   if (!$self) { $self->throw("no value for self in $value"); }
   if( defined $value ) {
       $self->{swq}->display_id($value);
   }
   return $self->{swq}->display_id();
}

=head2 seq

 Title   : seq()
 Usage   : $string    = $obj->seq(); _or_
	$obj->seq("atctatcatca");
 Function: Returns the sequence that is contained in the imbedded in the
	PrimarySeq object within the Quality object
 Returns : A scalar (the seq() value for the imbedded PrimarySeq object.)
 Args    : If a scalar is provided, the Quality object will
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
		$self->{swq}->seq($value);
	}
	return $self->{swq}->seq();
}

=head2 qual()

 Title   : qual()
 Usage   : @quality_values  = @{$obj->qual()}; _or_
	$obj->qual("10 10 20 40 50");
 Function: Returns the quality as imbedded in the PrimaryQual object
	within the Quality object.
 Returns : A reference to an array containing the quality values in the 
	PrimaryQual object.
 Args    : If a scalar is provided, the Quality object will
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
		$self->{swq}->qual($value);
	}
	return $self->{swq}->qual();
}




=head2 length()

 Title   : length()
 Usage   : $length = $seqWqual->length();
 Function: Get the length of the Quality sequence/quality.
 Returns : Returns the length of the sequence and quality
 Args    : None.

=cut

sub length {
    my $self = shift;
    return $self->seq_obj->length;
}


=head2 qual_obj

 Title   : qual_obj($different_obj)
 Usage   : $qualobj = $seqWqual->qual_obj(); _or_
	$qualobj = $seqWqual->qual_obj($ref_to_primaryqual_obj);
 Function: Get the Qualilty object that is imbedded in the
	Quality object or if a reference to a PrimaryQual object
	is provided, set this as the PrimaryQual object imbedded in the
	Quality object.
 Returns : A reference to a Bio::Seq::Quality object.

Identical to L<seq_obj>.

=cut

sub qual_obj {
    my ($self,$value) = @_;
#    return $self->{swq}->qual_obj($value);
    return $self->{swq};
}


=head2 seq_obj

 Title   : seq_obj()
 Usage   : $seqobj = $seqWqual->seq_obj(); _or_
	$seqobj = $seqWqual->seq_obj($ref_to_primary_seq_obj);
 Function: Get the PrimarySeq object that is imbedded in the
	Quality object or if a reference to a PrimarySeq object is
	provided, set this as the PrimarySeq object imbedded in the
	Quality object.
 Returns : A reference to a Bio::PrimarySeq object.

=cut

sub seq_obj {
    my ($self,$value) = @_;
    return $self->{swq};
}

=head2 _set_descriptors

 Title   : _set_descriptors()
 Usage   : $seqWqual->_qual_obj($qual,$seq,$id,$acc,$pid,$desc,$given_id,
	$alphabet);
 Function: Set the descriptors for the Quality object. Try to
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
    $self->{swq}->_seq_descriptors($qual,$seq,$id,$acc,$pid,
				   $desc,$given_id,$alphabet);
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

=head2 sub_peak_index($start,$end)

 Title   : sub_peak_index($start,$end)
 Usage   : @peak_indices = @{$obj->sub_peak_index(10,20);
 Function: returns the trace index values from $start to $end, where the
        first value is 1 and the number is inclusive, ie 1-2 are the
	first two trace indices for this channel.
 Returns : A reference to an array.
 Args    : a start position and an end position

=cut

sub sub_peak_index {
   my ($self,$start,$end) = @_;
   if( $start > $end ){
       $self->throw("in sub_peak_index, start [$start] has to be greater than end [$end]");
   }

   if( $start <= 0 || $end > $self->length ) {
       $self->throw("You have to have start positive and length less than the total length of sequence [$start:$end] Total ".$self->length."");
   }

   # remove one from start, and then length is end-start

   $start--;
     $end--;
     my @sub_peak_index_array = @{$self->{peak_indices}}[$start..$end];

     #   return substr $self->seq(), $start, ($end-$start);
     return \@sub_peak_index_array;

}

=head2 sub_trace($start,$end)

 Title   : sub_trace($base_channel,$start,$end)
 Usage   : @trace_values = @{$obj->sub_trace('a',10,20)};
 Function: returns the trace values from $start to $end, where the
        first value is 1 and the number is inclusive, ie 1-2 are the
	first two bases of the sequence. Start cannot be larger than
	end but can be e_peak_index.
 Returns : A reference to an array.
 Args    : a start position and an end position

=cut

sub sub_trace {
   my ($self,$base_channel,$start,$end) = @_;
   if( $start > $end ){
       $self->throw("in sub_trace, start [$start] has to be greater than end [$end]");
   }

   if( $start <= 0 || $end > $self->trace_length() ) {
       $self->throw("You have to have start positive and length less than the total length of traces [$start:$end] Total ".$self->trace_length."");
   }

   # remove one from start, and then length is end-start

   $start--;
     $end--;
     my @sub_peak_index_array = @{$self->trace($base_channel)}[$start..$end];

     #   return substr $self->seq(), $start, ($end-$start);
     return \@sub_peak_index_array;

}

=head2 trace_length()

 Title   : trace_length()
 Usage   : $trace_length = $obj->trace_length();
 Function: Return the length of the trace if all four traces (atgc)
     are the same. Otherwise, throw an error.
 Returns : A scalar.
 Args    : none

=cut

sub trace_length {
    my $self = shift;
     if ( !$self->trace('a') || !$self->trace('t') || !$self->trace('g') || !$self->trace('c') ) {
           $self->warn("One or more of the trace channels are missing. Cannot give you a length.");

     } 
     my $lengtha = scalar(@{$self->trace('a')});
     my $lengtht = scalar(@{$self->trace('t')});
     my $lengthg = scalar(@{$self->trace('g')});
     my $lengthc = scalar(@{$self->trace('c')});
     if (($lengtha == $lengtht) && ($lengtha == $lengthg) && ($lengtha == $lengthc) ) {
          return $lengtha;
     }
     $self->warn("Not all of the trace indices are the same length".
          " Here are their lengths: a: $lengtha t:$lengtht ".
          " g: $lengthg c: $lengthc");
}



=head2 sub_trace_object($start,$end)

 Title   : sub_trace_object($start,$end)
 Usage   : $smaller_object = $object->sub_trace_object('1','100');
 Function: Get a subset of the sequence, its quality, and its trace.
 Returns : A reference to a Bio::Seq::SequenceTrace object
 Args    : a start position and an end position
 Notes   : 
     - the start and end position refer to the positions of _bases_.
     - for example, to get a sub SequenceTrace for bases 5-10,
          use this routine.
          - you will get the bases, qualities, and the trace values
          - you can then use this object to synthesize a new scf
               using seqIO::scf.

=cut

sub sub_trace_object {
     my ($self,$start,$end) = @_;
          my ($start2,$end2);
        my @subs = @{$self->sub_peak_index($start,$end)};
        $start2 = shift(@subs);
        $end2 =  pop(@subs);
     my $new_object =  Bio::Seq::SequenceTrace->new(
               -swq =>   Bio::Seq::Quality->new(
                             -seq => $self->subseq($start,$end),
                             -qual     =>   $self->subqual($start,$end),
                             -id    =>   $self->id()
                         ),
             -trace_a  => $self->sub_trace('a',$start2,$end2),
             -trace_t  => $self->sub_trace('t',$start2,$end2),
             -trace_g  => $self->sub_trace('g',$start2,$end2),
             -trace_c  => $self->sub_trace('c',$start2,$end2),
             -peak_indices =>   $self->sub_peak_index($start,$end)

        );
     $new_object->set_accuracies();
     $new_object->_reset_peak_indices();
     return $new_object;
}

=head2 _synthesize_traces()

 Title   : _synthesize_traces()
 Usage   : $obj->_synthesize_traces();
 Function: Synthesize false traces for this object.
 Returns : Nothing.
 Args    : None.
 Notes   : This method is intended to be invoked when this
     object is created with a SWQ object- that is to say that
     there is a sequence and a set of qualities but there was
     no actual trace data.

=cut

sub _synthesize_traces {
     my ($self) = shift;
     $self->peak_indices(qw());
#ml     my $version = 2;
          # the user should be warned if traces already exist
          #
          #
#ml     ( my $sequence = $self->seq() ) =~ tr/a-z/A-Z/;
#ml     my @quals = @{$self->qual()};
#ml     my $info;
         # build the ramp for the first base.
         # a ramp looks like this "1 4 13 29 51 71 80 71 51 29 13 4 1" times the quality score.
         # REMEMBER: A C G T
         # note to self-> smooth this thing out a bit later
     my $ramp_data;
    @{$ramp_data->{'ramp'}} = qw( 1 4 13 29 51 75 80 75 51 29 13 4 1 );
         # the width of the ramp
    $ramp_data->{'ramp_width'} = scalar(@{$ramp_data->{'ramp'}});
         # how far should the peaks overlap?
    $ramp_data->{'ramp_overlap'} = 1;
          # where should the peaks be located?
    $ramp_data->{'peak_at'} = 7;
    $ramp_data->{'ramp_total_length'} =
          $self->seq_obj()->length() * $ramp_data->{'ramp_width'}
          - $self->seq_obj()->length() * $ramp_data->{'ramp_overlap'};
    my $pos;
    my $total_length = $ramp_data->{ramp_total_length};
     $self->initialize_traces("0",$total_length+2);
         # now populate them
    my ($current_base,$place_base_at,$peak_quality,$ramp_counter,$current_ramp,$ramp_position);
#ml    my $sequence_length = $self->length();
    my $half_ramp = int($ramp_data->{'ramp_width'}/2);
    for ($pos = 0; $pos<$self->length();$pos++) {
          $current_base = uc $self->seq_obj()->subseq($pos+1,$pos+1);
               # print("Synthesizing the ramp for $current_base\n");
          my $all_bases = "ATGC";
          $peak_quality = $self->qual_obj()->qualat($pos+1);
                    # where should the peak for this base be placed? Modeled after a mktrace scf
          $place_base_at = ($pos * $ramp_data->{'ramp_width'}) -
                      ($pos * $ramp_data->{'ramp_overlap'}) -
                   $half_ramp + $ramp_data->{'ramp_width'} - 1;
               # print("Placing this base at this position: $place_base_at\n");
          push @{$self->peak_indices()},$place_base_at;
          $ramp_position = $place_base_at - $half_ramp;
          if ($current_base =~ "N" ) {
               $current_base = "A";
          }
          for ($current_ramp = 0; $current_ramp < $ramp_data->{'ramp_width'};  $current_ramp++) {
                    # print("Placing a trace value here: $current_base ".($ramp_position+$current_ramp+1)." ".$peak_quality*$ramp_data->{'ramp'}->[$current_ramp]."\n");
             $self->trace_value_at($current_base,$ramp_position+$current_ramp+1,$peak_quality*$ramp_data->{'ramp'}->[$current_ramp]);
          }
          $self->peak_index_at($pos+1,
              $place_base_at+1
          );
#ml          my $other_bases = $self->_get_other_bases($current_base);
          # foreach ( split('',$other_bases) ) {
          #          push @{$self->{'text'}->{"v3_base_accuracy"}->{$_}},0;
          #}
    }
}





=head2 _dump_traces($transformed)

 Title   : _dump_traces("transformed")
 Usage   : &_dump_traces($ra,$rc,$rg,$rt);
 Function: Used in debugging. Prints all traces one beside each other.
 Returns : Nothing.
 Args    : References to the arrays containing the traces for A,C,G,T.
 Notes   : Beats using dumpValue, I'll tell ya. Much better then using
           join' ' too.
     - if a scalar is included as an argument (any scalar), this
     procedure will dump the _delta'd trace. If you don't know what
     that means you should not be using this.

=cut

#'
sub _dump_traces {
    my ($self) = @_;
    my (@sA,@sT,@sG,@sC);
    print ("Count\ta\tc\tg\tt\n");
     my $length = $self->trace_length();
    for (my $curr=1; $curr <= $length; $curr++) {
     print(($curr-1)."\t".$self->trace_value_at('a',$curr).
                "\t".$self->trace_value_at('c',$curr).
                "\t".$self->trace_value_at('g',$curr).
                "\t".$self->trace_value_at('t',$curr)."\n");
    }
    return;
}

=head2 _initialize_traces()

 Title   : _initialize_traces()
 Usage   : $trace_object->_initialize_traces();
 Function: Creates empty arrays to hold synthetic trace values.
 Returns : Nothing.
 Args    : None.

=cut

sub initialize_traces {
     my ($self,$value,$length) = @_;
     foreach (qw(a t g c)) {
          my @temp;
          for (my $count=0; $count<$length; $count++) {
               $temp[$count] = $value;
          }
          $self->trace($_,\@temp);
     }
}

=head2 trace_value_at($channel,$position)

 Title   : trace_value_at($channel,$position)
 Usage   : $value = $trace_object->trace_value_at($channel,$position);
 Function: What is the value of the trace for this base at this position?
 Returns : A scalar represnting the trace value here.
 Args    : a base channel (a,t,g,c)
           a position ( < $trace_object->trace_length() )

=cut

sub trace_value_at {
     my ($self,$channel,$position,$value) = @_;
     if ($value) {
          $self->trace($channel)->[$position] = $value;
     }
     return $self->sub_trace($channel,($position),($position))->[0];
}

sub _deprecated_get_scf_version_2_base_structure {
          # this sub is deprecated- check inside SeqIO::scf
     my $self = shift;
     my (@structure,$current);
     my $length = $self->length();
     for ($current=1; $current <= $self->length() ; $current++) {
           my $base_here = $self->seq_obj()->subseq($current,$current);
          $base_here = lc($base_here);
          my $probabilities;
          $probabilities->{$base_here} = $self->qual_obj()->qualat($current);
          my $other_bases = "atgc";
          my $empty = "";
          $other_bases =~ s/$base_here/$empty/e;
          foreach ( split('',$other_bases) ) {
               $probabilities->{$_} = "0";
          }
          @structure = (
               @structure,
              $self->peak_index_at($current),
              $probabilities->{'a'},
              $probabilities->{'t'},
              $probabilities->{'g'},
              $probabilities->{'c'}
         ); 
          
     }
     return \@structure;
}

sub _deprecated_get_scf_version_3_base_structure {
     my $self = shift;
     my $structure;
     $structure = join('',$self->peak_indices());
     return $structure;
}


=head2 accuracies($channel,$position)

 Title   : trace_value_at($channel,$position)
 Usage   : $value = $trace_object->trace_value_at($channel,$position);
 Function: What is the value of the trace for this base at this position?
 Returns : A scalar representing the trace value here.
 Args    : a base channel (a,t,g,c)
           a position ( < $trace_object->trace_length() )

=cut


sub accuracies {
     my ($self,$channel,$value) = @_;
     if ($value) {
          if (ref($value) eq "ARRAY") {
               $self->{accuracies}->{$channel} = $value;
          }
          else {
               my @acc = split(' ',$value);
               $self->{accuracies}->{$channel} = \@acc;
          }
     }
     return $self->{accuracies}->{$channel};
}


=head2 set_accuracies()

 Title   : set_sccuracies()
 Usage   : $trace_object->set_accuracies();
 Function: Take a sequence's quality and synthesize proper scf-style
     base accuracies that can then be accessed with
     accuracies("a") or something like it.
 Returns : Nothing.
 Args    : None.

=cut

sub set_accuracies {
     my $self = shift;
     my $count = 0;
     my $length = $self->length();
     for ($count=1; $count <= $length; $count++) {
          my $base_here = $self->seq_obj()->subseq($count,$count);
          my $qual_here = $self->qual_obj()->qualat($count);
          $self->accuracy_at($base_here,$count,$qual_here);
          my $other_bases = $self->_get_other_bases($base_here);
          foreach (split('',$other_bases)) {
               $self->accuracy_at($_,$count,"null");
          }
     }
}


=head2 scf_dump()

 Title   : scf_dump()
 Usage   : $trace_object->scf_dump();
 Function: Prints out the contents of the structures representing
     the SequenceTrace in a manner similar to io_lib's scf_dump.
 Returns : Nothing. Prints out the contents of the structures
     used to represent the sequence and its trace.
 Args    : None.
 Notes   : Used in debugging, obviously.

=cut

sub scf_dump {
     my $self = shift;
     my $count;
     for ($count=1;$count<=$self->length();$count++) {
          my $base_here = lc($self->seq_obj()->subseq($count,$count));
          print($base_here." ".sprintf("%05d",$self->peak_index_at($count))."\t");
          foreach (sort qw(a c g t)) {
               print(sprintf("%03d",$self->accuracy_at($_,$count))."\t");
          }
          print("\n");
     }
     $self->_dump_traces();
}

=head2 _get_other_bases($this_base)

 Title   : _get_other_bases($this_base)
 Usage   : $other_bases = $trace_object->_get_other_bases($this_base);
 Function: A utility routine to return bases other then the one provided.
     I was doing this over and over so I put it here.
 Returns : Three of a,t,g and c.
 Args    : A base (atgc)
 Notes   : $obj->_get_other_bases("a") returns "tgc"

=cut

sub _get_other_bases {
     my ($self,$this_base) = @_;
     $this_base = lc($this_base);
     my $all_bases = "atgc";
     my $empty = "";
     $all_bases =~ s/$this_base/$empty/e;
     return $all_bases;
}


=head2 accuracy_at($base,$position)

 Title   : accuracy_at($base,$position)
 Usage   : $accuracy = $trace_object->accuracy_at($base,$position);
 Function: 
 Returns : Returns the accuracy of finding $base at $position.
 Args    : 1. a base channel (atgc) 2. a value to _set_ the accuracy
 Notes   : $obj->_get_other_bases("a") returns "tgc"

=cut


sub accuracy_at {
     my ($self,$base,$position,$value) = @_;
     $base = lc($base);
     if ($value) {
          if ($value eq "null") {
               $self->{accuracies}->{$base}->[$position-1] = "0";
          }
          else {
               $self->{accuracies}->{$base}->[$position-1] = $value;
          }
     }
     return $self->{accuracies}->{$base}->[$position-1];
}

1;

