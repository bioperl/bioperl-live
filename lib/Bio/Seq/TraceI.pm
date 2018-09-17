# BioPerl module for Bio::Seq::TraceI
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

Bio::Seq::TraceI - Interface definition for a Bio::Seq::Trace

=head1 SYNOPSIS

    # get a Bio::Seq::Qual compliant object somehow
          $st = &get_object_somehow();

    # to test this is a seq object
          $st->isa("Bio::Seq::TraceI") 
               || $obj->throw("$obj does not implement the Bio::Seq::TraceI interface");

          # set the trace for T to be @trace_points
     my $arrayref = $st->trace("T",\@trace_points);
          # get the trace points for "C"
     my $arrayref = $st->trace("C");
          # get a subtrace for "G" from 10 to 100 
     $arrayref = $st->subtrace("G",10,100);
          # what is the trace value for "A" at position 355?
     my $trace_calue = $st->traceat("A",355);
          # create a false trace for "A" with $accuracy
     $arrayref = $st->false_trace("A",Bio::Seq::Quality, $accuracy);
          # does this trace have entries for each base?
     $bool = $st->is_complete();
          # how many entries are there in this trace?
     $length = $st->length();



=head1 DESCRIPTION

This object defines an abstract interface to basic trace information. This
information may have come from an ABI- or scf- formatted file or may have been
made up.

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

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::TraceI;
use strict;
use Carp;
use Dumpvalue;
use Bio::Root::RootI;

=head1 Implementation Specific Functions

These functions are the ones that a specific implementation must
define.

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
   my ($self) = @_;
   if( $self->can('throw') ) {
       $self->throw("Bio::Seq::TraceI definition of trace - implementing class did not provide this method");
   } else {
       confess("Bio::Seq::TraceI definition of trace - implementing class did not provide this method");
   }
}

=head2 subtrace($base,$start,$end)

 Title   : subtrace($base,$start,$end)
 Usage   : @subset_of_traces = @{$obj->subtrace("A",10,40)};
 Function: returns the trace values from $start to $end, where the
	first value is 1 and the number is inclusive, ie 1-2 are the first
	two trace values of this base. Start cannot be larger than end but can
	be equal.
 Returns : A reference to an array.
 Args    : $base: "A","T","G" or "C"
          $start: a start position
          $end  : an end position

=cut

sub subtrace {
   my ($self) = @_;

   if( $self->can('throw') ) {
       $self->throw("Bio::Seq::TraceI definition of subtrace - implementing class did not provide this method");
   } else {
       confess("Bio::Seq::TraceI definition of subtrace - implementing class did not provide this method");
   }

}

=head2 can_call_new()

 Title   : can_call_new()
 Usage   : if( $obj->can_call_new ) {
             $newobj = $obj->new( %param );
	 }
 Function: can_call_new returns 1 or 0 depending on whether an
	implementation allows new constructor to be called. If a new
	constructor is allowed, then it should take the followed hashed
	constructor list.
           $myobject->new( -qual => $quality_as_string,
			   -display_id  => $id,
			   -accession_number => $accession,
			   );
 Example :
 Returns : 1 or 0
 Args    :


=cut

sub can_call_new{
   my ($self,@args) = @_;
   # we default to 0 here
   return 0;
}

=head2 traceat($channel,$position)

 Title   : qualat($channel,$position)
 Usage   : $trace = $obj->traceat(500);
 Function: Return the trace value at the given location, where the
        first value is 1 and the number is inclusive, ie 1-2 are the first
        two bases of the sequence. Start cannot be larger than end but can
        be equal.
 Returns : A scalar.
 Args    : A base and a position.

=cut

sub traceat {
   my ($self,$value) = @_;
   if( $self->can('warn') ) {
       $self->warn("Bio::Seq::TraceI definition of traceat - implementing class did not provide this method");
   } else {
       warn("Bio::Seq::TraceI definition of traceat - implementing class did not provide this method");
   }
   return '';
} 

=head2 length()

 Title   : length()
 Usage   : $length = $obj->length("A");
 Function: Return the length of the array holding the trace values for the "A"
     channel. A check should be done to make sure that this Trace object
     is_complete() before doing this to prevent hazardous results.
 Returns : A scalar (the number of elements in the quality array).
 Args    : If used, get the traces from that channel. Default to "A"

=cut

sub length {
   my ($self)= @_;
   if( $self->can('throw') ) {
       $self->throw("Bio::Seq::TraceI definition of length - implementing class did not provide this method");
   } else {
       confess("Bio::Seq::TraceI definition of length - implementing class did not provide this method");
   }
}

=head2 trace_indices($new_indices)

 Title   : trace_indices($new_indices)
 Usage   : $indices = $obj->trace_indices($new_indices);
 Function: Return the trace iindex points for this object.
 Returns : A scalar
 Args    : If used, the trace indices will be set to the provided value.

=cut

sub trace_indices {
   my ($self)= @_;
   if( $self->can('throw') ) {
       $self->throw("Bio::Seq::TraceI definition of trace_indices - implementing class did not provide this method");
   } else {
       confess("Bio::Seq::TraceI definition of trace_indices - implementing class did not provide this method");
   }
}




1;
