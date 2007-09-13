# $Id$
#
# BioPerl module for Bio::Align::PairwiseStatistics
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

=head1 DESCRIPTION

Calculate pairwise statistics.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

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
   if( ! defined $aln || ! $aln->isa('Bio::Align::AlignI') ) {
       $self->warn("Must provide a Bio::Align::AlignI compliant object to Bio::Align::PairwiseStatistics");
       return 0;
   } elsif( $aln->no_sequences != 2 ) { 
       $self->warn("only pairwise calculations currently supported ". $aln->no_sequences."\n");
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
	$self->warn("Must provide a Bio::Align::AlignI compliant object to Bio::Align::PairwiseStatistics");
	return 0;
    } elsif( $aln->no_sequences != 2 ) { 
	$self->warn("only pairwise calculations currently supported");
    }
   my (@seqs);
   foreach my $seq ( $aln->each_seq) {
       push @seqs, [ split(//,$seq->seq())];
   }
   my $firstseq = shift @seqs;
#    my $secondseq = shift @seqs;
   my $diffcount = 0;
   for (my $i = 0;$i<$aln->length; $i++ ) {
       next if( $firstseq->[$i]  =~ /^$GapChars$/);
       foreach my $seq ( @seqs ) {
	   next if( $seq->[$i]  =~ /^$GapChars$/);
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
    if( ! defined $aln || ! $aln->isa('Bio::Align::AlignI') ) {
	$self->warn("Must provide a Bio::Align::AlignI compliant object to Bio::Align::PairwiseStatistics");
	return 0;
    }
   my $gapline = $aln->gap_line;
   # this will count the number of '-' characters
   return $gapline =~ tr/-/-/;
}

1;
