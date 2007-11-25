# $Id$
#
# BioPerl module for Bio::Seq::Quality
#
# Cared for by Heikki Lehvaslaiho
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::Quality - Implementation of sequence with residue quality and trace values

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
or set subvalues (e.g. subqual() ). See L<Bio::Seq::MetaI> for more details.

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

This implementation is based on Bio::Seq::Meta::Array that is a subclass
of Bio::PrimarySeq that stores any number of meta information in
unnamed arrays.

Here we assume that two meta sets, called 'qual' and 'trace_indices' are
attached to a sequence. (But there is nothing that prevents you to add
as many named meta sets as you need using normal meta() methods).

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
have been added (L<quality_is_flush>, L<trace_is_flush>). To force
the old functinality, set L<force_flush> to true.

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

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Heikki Lehvaslaiho

Email heikki-at-bioperl-dot-org

=head1 CONTRIBUTORS

Chad Matsalla, bioinformatics@dieselwurks.com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::Quality;
use vars qw($DEFAULT_NAME $GAP $META_GAP);
use strict;

#use overload '""' => \&to_string;

use base qw(Bio::Seq::Meta::Array);


BEGIN {

    $DEFAULT_NAME = 'DEFAULT';
    $GAP = '-';
    $META_GAP = ' ';
}

=head2 new

 Title   : new
 Usage   : $metaseq = Bio::Seq::Quality->new
	        ( -qual => '0 1 2 3 4 5 6 7 8 9 11 12',
                  -trace => '0 5 10 15 20 25 30 35 40 45 50 55',
                  -seq =>  'atcgatcgatcg',
	          -id  => 'human_id',
	          -accession_number => 'S000012',
	        );
 Function: Constructor for Bio::Seq::Quality class.
           Note that you can provide an empty quality and trace strings.
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
			      TRACE_DATA
                              )],
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

=cut

sub qual {
    my $self = shift;
    my $value = shift;
    $value = $value->qual if ref($value) and ref($value) ne 'ARRAY' and $value->isa('Bio::Seq::PrimaryQual');
    $self->named_meta($DEFAULT_NAME, $value);
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
    return join ' ', @{shift->submeta(@_)};
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
    return join ' ', @{shift->submeta(@_)};
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
   return $self->named_meta_length('DEFAULT');
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
    return shift->named_meta('trace', shift);
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
 Usage    : @trace_values = $obj->get_trace_graph( -trace => 'a', -scale => 100)
 Function : Returns array of raw trace values for a trace file, or false if no trace data exists.
            Requires a value for trace to get either the a, g, c or t trace information, and an
	    optional value for scale, which rescales the data between 0 and the provided value, 
	    a scale value of '0' performs no scaling
 Returns  : Array or 0
 Args     : string, trace to retrieve, one of a, g, c or t
            integer, scale, for scaling of trace between 0 and scale, or 0 for no scaling, optional

=cut

sub get_trace_graph
{
	my $self = shift;
        my($trace, $scale) =
	$self->_rearrange([qw(TRACE
			      SCALE
			      )],
			  @_);
	unless (grep { lc($trace) eq $_ } ('a', 'g', 'c', 't')) { return 0 }
	unless (defined($self->{trace_data})) { return 0 }
	$trace = lc($trace) . "_trace";
	my @trace_data = @{$self->{trace_data}->{$trace}};
	my $max = $self->{trace_data}->{max_height};
	if (defined($scale) and $scale != 0)
	{
		@trace_data = map { $_ / $max * $scale } @trace_data;
	}
	return @trace_data;
}	


################## deprecated methdods ##################


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

