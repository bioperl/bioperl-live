#
# BioPerl module for Bio::AlignIO::largemultifasta

#   based on the Bio::SeqIO::largefasta module
#       by Ewan Birney <birney@ebi.ac.uk>
#       and Lincoln Stein  <lstein@cshl.org>
#
#       and the SimpleAlign.pm module of Ewan Birney
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself
# _history
# January 20, 2004
# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::largemultifasta - Largemultifasta MSA Sequence
input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::AlignIO> class.

=head1 DESCRIPTION

This object can transform L<Bio::SimpleAlign> objects to and from
largemultifasta flat file databases.  This is for the fasta sequence
format NOT FastA analysis program.  To process the pairwise alignments
from a FastA (FastX, FastN, FastP, tFastA, etc) use the Bio::SearchIO
module.

Reimplementation of Bio::AlignIO::fasta modules so that creates
temporary files instead of keeping the whole sequences in memory.

=head1 FEEDBACK

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

=head1 AUTHORS - Albert Vilella, Heikki Lehvaslaiho

Email: avilella-at-gmail-dot-com, heikki-at-bioperl-dot-org


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::largemultifasta;
use strict;

use Bio::Seq::LargeLocatableSeq;
use Bio::Seq::SeqFactory;

use base qw(Bio::AlignIO Bio::SeqIO Bio::SimpleAlign);


sub _initialize {
    my($self,@args) = @_;
    $self->SUPER::_initialize(@args);
    if( ! defined $self->sequence_factory ) {
        $self->sequence_factory(Bio::Seq::SeqFactory->new(
            -verbose => $self->verbose(),
            -type => 'Bio::Seq::LargeLocatableSeq'
        ));
    }
}

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream while taking care
           of the length
 Returns : Bio::Seq object
 Args    : NONE

=cut

sub next_seq {
    my ($self) = @_;
    my $largeseq = $self->sequence_factory->create(-alphabet=>$self->alphabet);
    my ($id,$fulldesc,$entry);
    my $count = 0;
    my $seen = 0;
    while( defined ($entry = $self->_readline) ) {
        if( $seen == 1 && $entry =~ /^\s*>/ ) {
            $self->_pushback($entry);
            return $largeseq;
        }
        if ( $entry eq '>' ) {
            $seen = 1; next;
        } elsif( $entry =~ /\s*>(.+?)$/ ) {
            $seen = 1;
            ($id,$fulldesc) = ($1 =~ /^\s*(\S+)\s*(.*)$/)
            or $self->warn("Can't parse fasta header");
            $largeseq->display_id($id);
            $largeseq->primary_id($id);
            $largeseq->desc($fulldesc);
        } else {
            $entry =~ s/\s+//g;
            $largeseq->add_sequence_as_string($entry);
        }
        (++$count % 1000 == 0 && $self->verbose() > 0) && print "line $count\n";
    }
    if( ! $seen ) { return; }
    return $largeseq;
}


=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream.
 Returns : L<Bio::Align::AlignI> object - returns 0 on end of file
        or on error
 Args    : NONE

=cut

sub next_aln {
    my $self = shift;
    my $largeseq;
    my $aln = Bio::SimpleAlign->new();
    while (defined ($largeseq = $self->next_seq) ) {
        $aln->add_seq($largeseq);
        $self->debug("sequence readed\n");
    }

    my $alnlen = $aln->length;
    foreach my $largeseq ( $aln->each_seq ) {
        if( $largeseq->length < $alnlen ) {
            my ($diff) = ($alnlen - $largeseq->length);
            $largeseq->seq("-" x $diff);
        }
    }

    return $aln if $aln->num_sequences;
    return;

}

=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in largemultifasta format
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object


=cut

sub write_aln {
    my ($self,@aln) = @_;
    my ($seq,$desc,$rseq,$name,$count,$length,$seqsub);

    foreach my $aln (@aln) {
        if( ! $aln || ! $aln->isa('Bio::Align::AlignI')  ) {
            $self->warn("Must provide a Bio::Align::AlignI object when calling write_aln");
            next;
        }
        foreach $rseq ( $aln->each_seq() ) {
            $name = $aln->displayname($rseq->get_nse());
            $seq  = $rseq->seq();
            $desc = $rseq->description || '';
            $self->_print (">$name $desc\n") or return ;
            $count =0;
            $length = length($seq);
            while( ($count * 60 ) < $length ) {
                $seqsub = substr($seq,$count*60,60);
                $self->_print ("$seqsub\n") or return ;
                $count++;
            }
        }
    }
    $self->flush if $self->_flush_on_write && defined $self->_fh;
    return 1;
}

1;
