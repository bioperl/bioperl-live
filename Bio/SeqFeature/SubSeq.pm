#
# BioPerl module for Bio::SeqFeature::SubSeq
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Copyright Florent Angly
#
# You may distribute this module under the same terms as perl itself


=head1 NAME

Bio::SeqFeature::SubSeq - Feature representing a subsequence

=head1 SYNOPSIS

  # SubSeq with implicit sequence
  use Bio::Seq;
  my $template = Bio::Seq->new( -seq => 'AAAAACCCCCGGGGGTTTTT' );
  $subseq = Bio::SeqFeature::Amplicon->new(
      -start    => 6,
      -end      => 15,
      -template => $template,
  );
  print "Subsequence is: ".$amplicon->seq->seq."\n"; # Should be 'CCCCCGGGGG'

  # SubSeq with explicit sequence
  use Bio::SeqFeature::Subseq;
  my $subseq = Bio::SeqFeature::Amplicon->new( 
      -seq => $seq_object,
  );

=head1 DESCRIPTION

Bio::SeqFeature::SubSeq extends L<Bio::SeqFeature::Generic> features to
represent a subsequence. When this feature is attached to a template sequence,
the sequence of feature is the subsequence of the template at this location. The
purpose of this class is to represent a sequence as a feature without having to
explictly store its sequence string.

Of course, you might have reasons to explicitly set a sequence. In that case,
note that the length of the sequence is allowed to not match the position of the
feature. For example, you can set sequence of length 10 in a SubSeq feature that
spans positions 30 to 50 of the template if you so desire.

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
the bugs and their resolution.  Bug reports can be submitted via 
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Florent Angly <florent.angly@gmail.com>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


package Bio::SeqFeature::SubSeq;

use strict;

use base qw(Bio::SeqFeature::Generic);

=head2 new

 Title   : new()
 Usage   : my $subseq = Bio::SeqFeature::SubSeq( -start => 1, -end => 10, -strand => -1);
 Function: Instantiate a new Bio::SeqFeature::SubSeq feature object
 Args    : -seq      , the sequence object or sequence string of the feature (optional)
           -template , attach the feature to the provided parent template sequence or feature (optional).
                       Note that you must specify the feature location to do this.
           -start, -end, -location, -strand and all other L<Bio::SeqFeature::Generic> argument can be used.
 Returns : A Bio::SeqFeature::SubSeq object

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($seq, $template) = $self->_rearrange([qw(SEQ TEMPLATE)], @args);
    if (defined $seq) {
        # Set the subsequence explicitly
        if (not ref $seq) {
            # Convert string to sequence object
            $seq = Bio::PrimarySeq->new( -seq => $seq );
        } else {
            # Sanity check
            if (not $seq->isa('Bio::PrimarySeqI')) {
                $self->throw("Expected a sequence object but got a '".ref($seq)."'\n");
            }
        }
        $self->seq($seq);
    }
    if ($template) {
        if ( not($self->start) || not($self->end) ) {
            $self->throw('Could not attach feature to template $template because'.
                         ' the feature location was not specified.');
        }

        # Need to attach to parent sequence and then add sequence feature
        my $template_seq;
        if ($template->isa('Bio::SeqFeature::Generic')) {
           $template_seq = $template->entire_seq;
        } elsif ($template->isa('Bio::SeqI')) {
           $template_seq = $template;
        } else {
           $self->throw("Expected a Bio::SeqFeature::Generic or Bio::SeqI object".
                        " as template, but got '$template'.");
        }
        $self->attach_seq($template_seq);
        $template->add_SeqFeature($self);

    }
    return $self;
}  


=head2 seq

 Title   : seq()
 Usage   : my $seq = $subseq->seq();
 Function: Get or set the sequence object of this SubSeq feature. If no sequence
           was provided, but the subseq is attached to a sequence, get the
           corresponding subsequence.
 Returns : A sequence object or undef
 Args    : None.

=cut

sub seq {
    my ($self, $value) = @_;
    if (defined $value) {
        # The sequence is explicit
        if ( not(ref $value) || not $value->isa('Bio::PrimarySeqI') ) {
            $self->throw("Expected a sequence object but got a '".ref($value)."'\n");
        }
        $self->{seq} = $value;
    }
    my $seq = $self->{seq};
    if (not defined $seq) {
        # The sequence is implied
        $seq = $self->SUPER::seq;
    }
    return $seq;
}


=head2 length

 Title   : seq()
 Usage   : my $length = $subseq->seq();
 Function: Get the length of the SubSeq feature. It is similar to the length()
           method of L<Bio::Generic::SeqFeature>, which computes length based
           on the location of the feature. However, if the feature was not
           given a location, return the length of the subsequence if possible.
 Returns : integer or undef
 Args    : None.

=cut

sub length {
    my ($self) = @_;
    # Try length from location first
    if ($self->start && $self->end) {
        return $self->SUPER::length();
    }
    # Then try length from subsequence
    my $seq = $self->seq;
    if (defined $seq) {
        return length $seq->seq;
    }
    # We failed
    return undef;
}



1;
