#
# BioPerl module for Bio::Align::PairwiseStatistics
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Align::PairwiseStatistics - Base statistic object for Pairwise Alignments

=head1 SYNOPSIS

  use strict;
  my $stats = Bio::Align::PairwiseStatistics->new();

  # get alignment object of two sequences somehow
  my $pwaln;
  print $stats->number_of_comparable_bases($pwaln);
  my $score = $stats->score_nuc($pwaln);


=head1 DESCRIPTION

Calculate pairwise statistics.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Align::PairwiseStatistics;
use vars qw($GapChars);
use strict;


BEGIN { $GapChars = '(\.|\-)'; }

use base qw(Bio::Root::Root Bio::Align::StatisticsI);

=head2 number_of_comparable_bases

 Title   : number_of_comparable_bases
 Usage   : my $bases = $stat->number_of_comparable_bases($aln);
 Function: Returns the count of the number of bases that can be
           compared (L) in this alignment ( length - gaps)
 Returns : integer
 Args    : L<Bio::Align::AlignI>


=cut

sub number_of_comparable_bases{
   my ($self,$aln) = @_;
  if ( ! defined $aln || ! $aln->isa('Bio::Align::AlignI') ) {
    $self->throw("Must provide a Bio::Align::AlignI compliant object to ".
      "Bio::Align::PairwiseStatistics");
       return 0;
  } elsif ( $aln->num_sequences != 2 ) { 
    $self->throw("Only pairwise calculations supported. Found ".
      $aln->num_sequences." sequences in alignment\n");
   }
   my $L = $aln->length - $self->number_of_gaps($aln);
   return $L;
}

=head2 number_of_differences

 Title   : number_of_differences
 Usage   : my $nd = $stat->number_of_distances($aln);
 Function: Returns the number of differences between two sequences
 Returns : integer
 Args    : L<Bio::Align::AlignI>


=cut

sub number_of_differences{
   my ($self,$aln) = @_;
    if( ! defined $aln || ! $aln->isa('Bio::Align::AlignI') ) {
    $self->throw("Must provide a Bio::Align::AlignI compliant object to ".
      "Bio::Align::PairwiseStatistics");
  } elsif ( $aln->num_sequences != 2 ) { 
    $self->throw("Only pairwise calculations supported. Found ".
      $aln->num_sequences." sequences in alignment\n");
    }
   my (@seqs);
  foreach my $seq ( $aln->each_seq ) {
       push @seqs, [ split(//,$seq->seq())];
   }
   my $firstseq = shift @seqs;
  #my $secondseq = shift @seqs;
   my $diffcount = 0;
   for (my $i = 0;$i<$aln->length; $i++ ) {
    next if ( $firstseq->[$i]  =~ /^$GapChars$/ );
       foreach my $seq ( @seqs ) {
      next if ( $seq->[$i]  =~ /^$GapChars$/ );
	   if( $firstseq->[$i] ne $seq->[$i] ) {
	       $diffcount++;
	   }
       }
   }
   return $diffcount;
}

=head2 number_of_gaps

 Title   : number_of_gaps
 Usage   : my $nd = $stat->number_of_gaps($aln);
 Function: Returns the number of gapped positions among sequences in alignment
 Returns : integer
 Args    : L<Bio::Align::AlignI>


=cut

sub number_of_gaps{
   my ($self,$aln) = @_;
  if ( ! defined $aln || ! $aln->isa('Bio::Align::AlignI') ) {
    $self->throw("Must provide a Bio::Align::AlignI compliant object to ".
      "Bio::Align::PairwiseStatistics");
  } elsif ( $aln->num_sequences != 2 ) { 
    $self->throw("Only pairwise calculations supported. Found ".
      $aln->num_sequences." sequences in alignment\n");
    }
   my $gapline = $aln->gap_line;
   # this will count the number of '-' characters
   return $gapline =~ tr/-/-/;
}


=head2 score_nuc

 Title   : score_nuc
 Usage   : my $score = $stat->score_nuc($aln);
             or
           my $score = $stat->score_nuc(
             -aln =>$aln,
             -match    => 1,
             -mismatch => -1,
             -gap_open => -1,
             -gap_ext  => -1
           );
 Function: Calculate the score of an alignment of 2 nucleic acid sequences. The
           scoring parameters can be specified. Otherwise the blastn default
           parameters are used: match = 2, mismatch = -3, gap opening = -5, gap
           extension = -2
 Returns : alignment score (number)
 Args    : L<Bio::Align::AlignI>
           match score [optional]
           mismatch score [optional]
           gap opening score [optional]
           gap extension score [optional]

=cut

sub score_nuc {
  my ($self, @args) = @_;
  my ( $aln, $match, $mismatch, $gap_open, $gap_ext) = $self->_rearrange( [qw(
    ALN MATCH MISMATCH GAP_OPEN GAP_EXT)], @args );
  if ( ! defined $aln || ! $aln->isa('Bio::Align::AlignI') ) {
    $self->throw("Must provide a Bio::Align::AlignI compliant object to ".
      "Bio::Align::PairwiseStatistics");
  } elsif ( $aln->num_sequences != 2 ) { 
    $self->throw("Only pairwise calculations supported. Found ".
      $aln->num_sequences." sequences in alignment\n");
  }
  my $seq1 = $aln->get_seq_by_pos(1);
  my $seq2 = $aln->get_seq_by_pos(2);
  if (! ( ($seq1->alphabet eq 'dna' || $seq1->alphabet eq 'rna') &&
    ($seq2->alphabet eq 'dna' || $seq2->alphabet eq 'rna') )) {
    $self->throw("Can only score nucleic acid alignments");
  }
  $match    ||=  2; # Blastn scoring defaults
  $mismatch ||= -3;
  $gap_open ||= -5;
  $gap_ext  ||= -2;
  my $score = 0;
  my $prevres1 = '-';
  my $prevres2 = '-';
  for (my $pos = 1 ; $pos <= $aln->length ; $pos++) {
    my $res1 = $seq1->subseq($pos, $pos);
    my $res2 = $seq2->subseq($pos, $pos);
    if (!($res1 eq '-' || $res2 eq '-')) { # no gap
      if ($res1 eq $res2) { # same residue
        $score += $match;
      } else { # other residue
        $score += $mismatch;
      }
    } else { # open or ext gap?
      my $open = 0;
      if (!($res1 eq '-' && $res2 eq '-')) { # exactly one gap
        my $prevres = $prevres1;
        $prevres = $prevres2 if $res2 eq '-';
        $open = 1 unless $prevres eq '-';
      } else { # 2 gaps
        $open = 1 unless $prevres1 eq '-' && $prevres2 eq '-';
      }
      if ($open) {
        $score += $gap_open; # gap opening
      } else {
        $score += $gap_ext; # gap extension
      }
    }
    $prevres1 = $res1;
    $prevres2 = $res2;
  }
  return $score;
}

1;
