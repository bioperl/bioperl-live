# BioPerl module for parsing MUMmer output (mummer tool only).
#
# You may distribute this module under the same terms as perl itself.
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::RefLocatableSeq - subtype of L<Bio::LocatableSeq|Bio::LocatableSeq>
to store partial locatable sequences within the reference sequence.

=head1 SYNOPSIS

 $obj = Bio::Seq::RefLocatableSeq->new();

 #obtain the reference sequence information
 ($refstart, $refend) = ($obj->ref_start, $obj->ref_end);

 #obtain the reference sequence match
 ($refseq, $queryseq) = ($obj->ref_seq, $obj->seq);

 #all other methods from Bio::LocatableSeq can be used to obtain location
 #information.

=head1 DESCRIPTION

Bio::Seq::RefLocatableSeq is a L<Bio::LocatableSeq|Bio::LocatableSeq> object
that stores partial locatable sequences within both the reference sequence
and the query sequence. This is typically used for alignment algorithms
that find partial alignments. For example, the MUMmer tool will find exact
matches within the reference and query sequences that may start anywhere in
both sequences.

Beyond that, this is only a minor modification to
L<Bio::LocatableSeq|Bio::LocatableSeq> to enable storing information about
the reference sequence.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an intergral part of the evolution of this and other Bioperl
modules. Send your comments and suggestions preferable to one of the Bioperl
mailing lists. Your participation is much appreciated.

 bioperl-l@bioperl.org                  - General discussion
 http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track the
bugs and their resolution. Bug reports can be submitted via the web:

 http://bugzilla.open-bio.org

=head1 AUTHORS

 Jason Switzer  - jswitzer@gmail.com
 Joshua Wu      - nike284@gmail.com
 Aimee Seuffer  - aimeeseuf@msn.com

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a '_'.

=cut

package Bio::Seq::RefLocatableSeq;

use base qw(Bio::LocatableSeq);
use strict;

=head2 new

 Title      : new
 Usage      : $obj = Bio::Seq::RefLocatableSeq->new(...)
 Function   : creates a new Bio::Seq::RefLocatableSeq object from a supplied
              sequence, query, and reference location.
 Args       : ref_id - ID of the reference sequence (for multiple reference
                       sequence alignments).
              ref_seq - Subsequence of the reference sequence.
              ref_start - Starting position of the reference sequence alignment.
              ref_end - Stopping position of the reference sequence alignment.

=cut

sub new {
   my ($self, @args) = @_;
   $self = $self->SUPER::new(@args, '-alphabet' => 'dna');

   my ($start, $end, $ref_id, $ref_seq, $ref_start, $ref_end, $strand) =
       $self->_rearrange([qw(START END REF_ID REF_SEQ REF_START REF_END STRAND)],
                         @args);
   $self->start($start);
   $self->end($end);
   $self->strand($strand);
   $self->ref_id($ref_id);
   $self->ref_seq($ref_seq);
   $self->ref_start($ref_start);
   #not overly sure what this did in LocatableSeq, but need to call it last
   $self->ref_end($ref_end);
   return $self;
}

=head2

 Title      : ref_start
 Usage      : $obj->ref_start(1234);
              $start = $obj->ref_start;
 Function   : get/set the reference sequence starting position
 Returns    : the starting position of the current alignment in the reference
              sequence.
 Args       : starting position to set, if provided

=cut

sub ref_start {
   my $self = shift;
   if(@_) {
      my $value = shift;
      $self->{'ref_start'} = $value;
   }
   return $self->{'ref_start'};
}

=head2

 Title      : ref_id
 Usage      : $obj->ref_id("SequnceId1");
              $id = $obj->ref_id;
 Function   : get/set the reference sequence ID
 Returns    : the ID of the current alignment in the reference sequence.
 Args       : ID to set, if provided

=cut

sub ref_id {
   my $self = shift;
   if(@_) {
      my $value = shift;
      $self->{'ref_id'} = $value;
   }
   return $self->{'ref_id'};
}

=head2

 Title      : ref_seq
 Usage      : $obj->ref_seq("ACCGTTC");
              $id = $obj->ref_seq;
 Function   : get/set the reference sequence alignment match.
 Returns    : the subsequence alignment in the reference sequence.
 Args       : subsequence alignment match to set, if provided

=cut

sub ref_seq {
   my $self = shift;
   if(@_) {
      my $value = shift;
      $self->{'ref_seq'} = $value;
   }
   return $self->{'ref_seq'};
}

=head2

 Title      : ref_end
 Usage      : $obj->ref_end(5234);
              $start = $obj->ref_end;
 Function   : get/set the reference sequence stopping position
 Returns    : the stopping position of the current alignment in the reference
              sequence.
 Args       : stopping position to set, if provided

=cut

sub ref_end {
   my $self = shift;
   if(@_) {
      my $value = shift;
      my $string = $self->ref_seq;
      if($string && $self->ref_start) {
         my $s2 = $string;
         $string =~ s/[.-]+//g;
         my $len = CORE::length $string;
         my $new_end = $self->ref_start + $len - 1;
         my $ref_id = $self->ref_id;
         if($new_end != $value && $self->verbose > 0) {
            my $msg = "In sequence $ref_id residue count gives value $len.\n";
            $msg .= "Overriding value [$value] with value $new_end";
            $msg .= " for Bio::RefLocatableSeq::ref_end().";
            $self->warn($msg);
            $value = $new_end;
         }
      }
      $self->{'ref_end'} = $value;
      return $self->{'ref_end'};
   }
}

1;

