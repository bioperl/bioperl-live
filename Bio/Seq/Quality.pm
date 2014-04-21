#
# BioPerl module for Bio::Seq::Quality
#
# Please direct questions and support issues to
# <bioperl-l@bioperl.org>
#
# Cared for by Heikki Lehvaslaiho
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::Quality - Implementation of sequence with residue quality
                    and trace values

=head1 SYNOPSIS

  use Bio::Seq::Quality;

  # input can be space delimited string or array ref
  my $qual = '0 1 2 3 4 5 6 7 8 9 11 12';
  my $trace = '0 5 10 15 20 25 30 35 40 45 50 55';

  my $seq = Bio::Seq::Quality->new
      ( -qual => $qual,
        -trace_indices => $trace,
        -seq =>  'atcgatcgatcg',
        -id  => 'human_id',
        -accession_number => 'S000012',
        -verbose => -1   # to silence deprecated methods
  );

  my $quals = $seq->qual; # array ref
  my $traces = $seq->trace;  # array ref

  my $quals = $seq->qual_text; # string
  my $traces = $seq->trace_text; # string


  # get sub values
  $quals = $seq->subqual(2, 3);  # array ref
  $traces = $seq->subtrace(2, 3); # array ref
  $quals = $seq->subqual_text(2, 3); # string
  $quals = $seq->subtrace_text(2, 3); # string

  # set sub values 
  $seq->subqual(2, 3, "9 9");
  $seq->subtrace(2, 3, "9 9");



=head1 DESCRIPTION

This object stores base quality values together with the sequence
string.

It is a reimplementation of Chad Matsalla's Bio::Seq::SeqWithQuality
module using Bio::Seq::MetaI.

The implementation is based on Bio::Seq::Meta::Array. qual() and
trace() are base methods to store and retrieve information that have
extensions to retrieve values as a scalar (e.g. qual_text() ), or get
or set subvalues (e.g. subqual() ). See L<Bio::Seq::MetaI> for more
details.

All the functional code is in Bio::Seq::Meta::Array.

There deprecated methods that are included for compatibility with
Bio::Seq::SeqWithQuality. These will print a warning unless verbosity
of the object is set to be less than zero.

=head2 Differences from Bio::Seq::SeqWithQuality

It is not possible to fully follow the interface of
Bio::Seq::SeqWithQuality since internally a Bio::Seq::SeqWithQuality
object is a composite of two independent objects: a Bio::PrimarySeq
object and Bio::Seq::PrimaryQual object. Both of these objects can be
created separately and merged into Bio::Seq::SeqWithQuality.

This implementation is based on Bio::Seq::Meta::Array that is a
subclass of Bio::PrimarySeq that stores any number of meta information
in unnamed arrays.

Here we assume that two meta sets, called 'qual' and 'trace_indices'
are attached to a sequence. (But there is nothing that prevents you to
add as many named meta sets as you need using normal meta() methods).

qual() is an alias to meta(), qualat($loc) is an alias to
submeta($loc,$loc).

trace_indices() in Bio::Seq::SeqWithQuality has been abbreviated to
trace() and is an alias to named_meta('trace').

You can create an object without passing any arguments to the
constructor (Bio::Seq::SeqWithQuality fails without alphabet). It will
warn about not being able to set alphabet unless you set verbosity of
the object to a negative value.

After the latest rewrite, the meta information sets (quality and
trace) no longer cover all the residues automatically. Methods to
check the length of meta information (L<quality_length>,
L<trace_length>)and to see if the ends are flushed to the sequence
have been added (L<quality_is_flush>, L<trace_is_flush>). To force the
old functinality, set L<force_flush> to true.

qual_obj() and seq_obj() methods do not exist!

Finally, there is only one set of descriptors (primary_id, display_id,
accession_number) for the object.


=head1 SEE ALSO

L<Bio::Seq::MetaI>,
L<Bio::Seq::Meta::Array>

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

=head1 AUTHOR - Heikki Lehvaslaiho

Email heikki-at-bioperl-dot-org

=head1 CONTRIBUTORS

Chad Matsalla, bioinformatics at dieselwurks dot com

Dan Bolser, dan dot bolser at gmail dot com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::Quality;
use strict;

use base qw(Bio::Seq::Meta::Array);

## Is this the right place (and way) to define this?
our $MASK_CHAR = 'X';
our $DEFAULT_NAME = 'DEFAULT';
our $GAP = '-';
our $META_GAP = ' ';

=head2 new

 Title   : new
 Usage   : $metaseq = Bio::Seq::Quality->new
	        ( -qual => '0 1 2 3 4 5 6 7 8 9 11 12',
                  -trace => '0 5 10 15 20 25 30 35 40 45 50 55',
                  -seq =>  'atcgatcgatcg',
	          -id  => 'human_id',
	          -accession_number => 'S000012',
	        );
 Function: Constructor for Bio::Seq::Quality class.  Note that you can
           provide an empty quality and trace strings.

 Returns : a new Bio::Seq::Quality object

=cut


sub new {
    my ($class, @args) = @_;

    my $self = $class->SUPER::new(@args);

    my($meta, $qual, $trace, $trace_indices, $trace_data) =
        $self->_rearrange([qw(META
                              QUAL
                              TRACE
                              TRACE_INDICES
                              TRACE_DATA)],
                          @args);

    $self->{'_meta'}->{$DEFAULT_NAME} = [];
    $self->{'_meta'}->{'trace'} = [];
    $self->{trace_data} = $trace_data;

    $meta && $self->meta($meta);
    $qual && $self->qual($qual);
    $trace && $self->named_meta('trace', $trace);
    $trace_indices && $self->named_meta('trace', $trace_indices);

    return $self;
}



## QUAL

=head2 qual

 Title   : qual
 Usage   : $qual_values  = $obj->qual($values_string);
 Function:

           Get and set method for the meta data starting from residue
           position one. Since it is dependent on the length of the
           sequence, it needs to be manipulated after the sequence.

           The length of the returned value always matches the length
           of the sequence.

 Returns : reference to an array of meta data
 Args    : new value, string or array ref or Bio::Seq::PrimaryQual, optional

Setting quality values resets the clear range.

=cut

sub qual {
    my $self = shift;
    my $value = shift;
    $value = $value->qual
        if ref($value) and ref($value) ne 'ARRAY' and
           $value->isa('Bio::Seq::PrimaryQual');
    $self->_empty_cache if $value;
    return $self->named_meta($DEFAULT_NAME, $value);
}

=head2 qual_text

 Title   : qual_text
 Usage   : $qual_values  = $obj->qual_text($values_arrayref);
 Function: Variant of meta() and qual()  guarantied to return a string
           representation  of meta data. For details, see L<meta>.
 Returns : a string
 Args    : new value, optional

=cut

sub qual_text {
    return join ' ', @{shift->named_submeta($DEFAULT_NAME, @_)};
}

=head2 subqual

 Title   : subqual
 Usage   : $subset_of_qual_values = $obj->subqual(10, 20, $value_string);
           $subset_of_qual_values = $obj->subqual(10, undef, $value_string);
 Function:

           Get and set method for meta data for subsequences.

           Numbering starts from 1 and the number is inclusive, ie 1-2
           are the first two residue of the sequence. Start cannot be
           larger than end but can be equal.

           If the second argument is missing the returned values
           should extend to the end of the sequence.

 Returns : A reference to an array
 Args    : integer, start position
           integer, end position, optional when a third argument present
           new value, optional

=cut

sub subqual {
    shift->named_submeta($DEFAULT_NAME, @_);
}

=head2 subqual_text

 Title   : subqual_text
 Usage   : $meta_values  = $obj->subqual_text(20, $value_string);
 Function: Variant of subqual() returning a stringified
           representation  of meta data. For details, see L<Bio::Seq::MetaI>.
 Returns : a string
 Args    : new value, optional

=cut

sub subqual_text {
    return join ' ', @{shift->named_submeta($DEFAULT_NAME, @_)};
}

=head2 quality_length

 Title   : quality_length()
 Usage   : $qual_len  = $obj->quality_length();
 Function: return the number of elements in the quality array
 Returns : integer
 Args    : -

=cut

sub quality_length {
   my ($self) = @_;
   return $self->named_meta_length($DEFAULT_NAME);
}

=head2 quality_is_flush

 Title   : quality_is_flush
 Usage   : $quality_is_flush  = $obj->quality_is_flush()
 Function: Boolean to tell if the trace length equals the sequence length.
           Returns true if force_flush() is set.
 Returns : boolean 1 or 0
 Args    : none

=cut

sub quality_is_flush {
    return shift->is_flush('quality');
}



## TRACE

=head2 trace

 Title   : trace
 Usage   : $trace_values  = $obj->trace($values_string);
 Function:

           Get and set method for the meta data starting from residue
           position one. Since it is dependent on the length of the
           sequence, it needs to be manipulated after the sequence.

           The length of the returned value always matches the length
           of the sequence.

 Returns : reference to an array of meta data
 Args    : new value, string or array ref, optional

=cut

sub trace {
    my $self = shift;
    my $value = shift;
    return $self->named_meta('trace', $value);
}

=head2 trace_text

 Title   : trace_text
 Usage   : $trace_values  = $obj->trace_text($values_arrayref);
 Function: Variant of meta() and trace()  guarantied to return a string
           representation  of meta data. For details, see L<meta>.
 Returns : a string
 Args    : new value, optional

=cut

sub trace_text {
    return join ' ', @{shift->named_submeta('trace', @_)};
}

=head2 subtrace

 Title   : subtrace
 Usage   : $subset_of_trace_values = $obj->subtrace(10, 20, $value_string);
           $subset_of_trace_values = $obj->subtrace(10, undef, $value_string);
 Function:

           Get and set method for meta data for subsequences.

           Numbering starts from 1 and the number is inclusive, ie 1-2
           are the first two residue of the sequence. Start cannot be
           larger than end but can be equal.

           If the second argument is missing the returned values
           should extend to the end of the sequence.

 Returns : A reference to an array
 Args    : integer, start position
           integer, end position, optional when a third argument present
           new value, optional


=cut

sub subtrace {
    return shift->named_submeta('trace', @_);
}

=head2 subtrace_text

 Title   : subtrace_text
 Usage   : $meta_values  = $obj->subtrace_text(20, $value_string);
 Function: Variant of subtrace() returning a stringified
           representation  of meta data. For details, see L<Bio::Seq::MetaI>.
 Returns : a string
 Args    : new value, optional

=cut

sub subtrace_text {
    return join ' ', @{shift->named_submeta('trace', @_)};
}

=head2 trace_length

 Title   : trace_length()
 Usage   : $trace_len  = $obj->trace_length();
 Function: return the number of elements in the trace set
 Returns : integer
 Args    : -

=cut

sub trace_length {
   my ($self) = @_;
   return $self->named_meta_length('trace');
}

=head2 trace_is_flush

 Title   : trace_is_flush
 Usage   : $trace_is_flush  = $obj->trace_is_flush()
 Function: Boolean to tell if the trace length equals the sequence length.
           Returns true if force_flush() is set.
 Returns : boolean 1 or 0
 Args    : none

=cut

sub trace_is_flush {
    return shift->is_flush('trace');
}



=head2 get_trace_graph

 Title    : get_trace_graph
 Usage    : @trace_values = $obj->get_trace_graph( -trace => 'a',
                                                   -scale => 100)
 Function : Returns array of raw trace values for a trace file, or
            false if no trace data exists.  Requires a value for trace
            to get either the a, g, c or t trace information, and an
            optional value for scale, which rescales the data between
            0 and the provided value, a scale value of '0' performs no
            scaling
 Returns  : Array or 0
 Args     : string, trace to retrieve, one of a, g, c or t integer,
            scale, for scaling of trace between 0 and scale, or 0 for
            no scaling, optional

=cut

sub get_trace_graph
{
	my $self = shift;
        my($trace, $scale) =
	$self->_rearrange([qw(TRACE
			      SCALE
			      )],
			  @_);
	unless (defined($self->{trace_data})) { return 0 }
	unless (grep { lc($trace) eq $_ } ('a', 'g', 'c', 't')) { return 0 }
	$trace = lc($trace) . "_trace";
	my @trace_data = exists $self->{trace_data}->{$trace} &&
	    ref $self->{trace_data}->{$trace} eq 'ARRAY' ?
	    @{$self->{trace_data}->{$trace}} : ();
	my $max = $self->{trace_data}->{max_height};
	if (defined($scale) and $scale != 0)
	{
	    @trace_data = map { $_ / $max * $scale } @trace_data;
	}
	return @trace_data;
}


=head2 threshold

  Title   : threshold
  Usage   : $qual->threshold($value);
  Function: Sets the quality threshold.
  Returns : an integer
  Args    : new value, optional

Value used by *clear_range* method below.

=cut

sub threshold {
    my $self = shift;
    my $value = shift;
    if (defined $value) {
	$self->throw("Threshold needs to be an integer [$value]")
	    unless $value =~ /^[-+]?\d+$/;
	$self->_empty_cache 
	    if defined $self->{_threshold} and $self->{_threshold} ne $value;
	$self->{_threshold} = $value;
    }
    return $self->{_threshold};
}


=head2 mask_below_threshold

  Title   : mask_below_threshold
  Usage   : $count = $obj->count_clear_ranges($threshold);
  Function: Counts number of ranges in the sequence where quality
            values are above the threshold
  Returns : count integer
  Args    : threshold integer, optional

Set threshold first using method L<threshold>.

=cut

sub mask_below_threshold {
    my $self = shift;
    my $threshold = shift;

    $self->threshold($threshold) if defined $threshold;

    # populate the cache if needed
    $self->_find_clear_ranges unless defined $self->{_ranges};

    my $maskSeq = $self->seq;
    my $maskQual = $self->qual;

    ## There must be a more efficient way than this!
    for(my $i=0; $i<length($maskSeq); $i++){
	#print join ("\t", $i, $maskQual->[$i]), "\n";
	substr($maskSeq, $i, 1, $MASK_CHAR)
	    if $maskQual->[$i] < $self->{_threshold};
    }

    ## This is the *wrong* way to do it!
    #for my $r (@{$self->{_ranges}} ){
    #	substr($maskSeq, $r->{start}, $r->{length}, $MASK_CHAR x $r->{length});
    #}

    return $maskSeq;
}

=head2 count_clear_ranges

  Title   : count_clear_ranges
  Usage   : $count = $obj->count_clear_ranges($threshold);
  Function: Counts number of ranges in the sequence where quality
            values are above the threshold
  Returns : count integer
  Args    : threshold integer, optional

Set threshold first using method L<threshold>.

=cut

sub count_clear_ranges {
    my $self = shift;
    my $threshold = shift;
    $self->threshold($threshold) if defined $threshold;

    # populate the cache if needed
    $self->_find_clear_ranges unless defined $self->{_ranges};

    return scalar @{$self->{_ranges}};
}

=head2 clear_ranges_length

  Title   : clear_ranges_length
  Usage   : $total_lenght = $obj->clear_ranges_length($threshold);
  Function: Return number of residues with quality values above
            the threshold in all clear ranges
  Returns : an integer
  Args    : threshold, optional

Set threshold first using method L<threshold>.

I think this method needs a better name! count_high_quality_bases? or
sum_clear_ranges?

=cut

sub clear_ranges_length {
    my $self = shift;
    my $threshold = shift;
    $self->threshold($threshold) if defined $threshold;

    # populate the cache if needed
    $self->_find_clear_ranges unless defined $self->{_ranges};

    my $sum;
    map {$sum += $_->{length}}  @{$self->{_ranges}};
    return $sum;
}

=head2 get_clear_range

  Title   : get_clear_range
  Usage   : $newqualobj = $obj->get_clear_range($threshold);
  Function: Return longest subsequence that has quality values above
            the given threshold, or a default value of 13
  Returns : a new Bio::Seq::Quality object
  Args    : threshold, optional

Set threshold first using method L<threshold>.

Note, this method could be implemented using some gaussian smoothing
of the quality scores. Currently one base below the threshold is
enough to end the clear range.

=cut

sub get_clear_range {
    my $self = shift;
    my $threshold = shift;
    $self->threshold($threshold) if defined $threshold;

    # populate the cache if needed
    $self->_find_clear_ranges unless defined $self->{_ranges};

    # fix for bug 2847
    return unless defined $self->{_ranges};

    # pick the longest
    for (sort {$b->{length} <=> $a->{length} } @{$self->{_ranges}} ){
        my $newqualobj = Bio::Seq::Quality->new(
                        -seq  => $self->subseq(  $_->{start}, $_->{end}),
                        -qual => $self->subqual($_->{start}, $_->{end}),
                        -id   => $self->id);
        
        $newqualobj->threshold($threshold);
    
        return $newqualobj;
    }
}



=head2 get_all_clean_ranges

  Title   : get_all_clean_ranges
  Usage   : @ranges = $obj->get_all_clean_ranges($minlength);
  Function: Return all ranges where quality values are above
            the threshold. Original ordering.
  Returns : an ordered array of new Bio::Seq::Quality objects
  Args    : minimum length , optional

Set threshold first using method L<threshold>.

=cut

sub get_all_clean_ranges {
    my $self = shift;
    my $minl = shift || 0;

    $self->throw("Mimimum length needs to be zero or a positive integer [$minl]")
        unless $minl =~ /^\+?\d+$/;

    # populate the cache if needed
    $self->_find_clear_ranges unless defined $self->{_ranges};

    # return in the order of occurrence
    my @ranges;
    for my $r (sort {$b->{start} <=> $a->{start} } @{$self->{_ranges}} ){
	next if $r->{length} < $minl;
	
	## Constructor should allow "-threshold => ..."!
	push @ranges, Bio::Seq::Quality->new
	    ( -seq => $self->subseq(  $r->{start}, $r->{end}),
	      -qual => $self->subqual($r->{start}, $r->{end}),
          -id   => $self->id
	    );
    }
    return @ranges;
}


#
# _find_clear_ranges: where range/threshold calculations happen
#

sub _find_clear_ranges {
    my $self = shift;
    my $qual = $self->qual;
    
    $self->throw("You need to set the threshold value first")
        unless defined $self->threshold;
    
    my $threshold = $self->threshold;
    
    my $rangeFlag = 0;
    
    for(my $i=0; $i<@$qual; $i++){
	## Are we currently within a clear range or not?
	if($rangeFlag){
	    ## Did we just leave the clear range?
	    if($qual->[$i]<$threshold){
		## Log the range
		my $range;
		$range->{end} = $i-1;
		$range->{start}  = $rangeFlag;
		$range->{length} = $i - $rangeFlag;
		push @{$self->{_ranges}}, $range;
		## and reset the range flag.
		$rangeFlag = 0;
	    }
	    ## else nothing changes
	}
	else{
	    ## Did we just enter a clear range?
	    if($qual->[$i]>=$threshold){
		## then set the range flag!
		$rangeFlag = $i;
	    }
	    ## else nothing changes
	}
    }
    ## Did we exit the last clear range?
    if($rangeFlag){
	my $i = scalar(@$qual);
	## Log the range
	my $range;
	$range->{end} = $i-1;
	$range->{start}  = $rangeFlag;
	$range->{length} = $i - $rangeFlag;
	push @{$self->{_ranges}}, $range;
    }
    
    1;
}


sub _empty_cache {
    my $self = shift;
    undef $self->{_ranges};
}




################## deprecated methods ##################


sub trace_indices {
    my $self = shift;
    return $self->named_meta('trace');
}

sub trace_index_at {
    my ($self, $val) =@_;
    return shift @{$self->named_submeta('trace', $val, $val)};
}


sub sub_trace_index {
    my $self = shift; 
    return $self->named_submeta('trace', @_);
}


sub qualat {
    my ($self, $val) =@_;
    return shift @{$self->submeta($val, $val)};
}


sub baseat {
    my ($self,$val) = @_;
    return $self->subseq($val,$val);
}



1;

