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

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Align::PairwiseStatistics;
use vars qw(@ISA $GapChars);
use strict;

use Bio::Align::StatisticsI;
use Bio::Root::Root;

BEGIN { $GapChars = '(\.|\-)'; }

@ISA = qw(Bio::Root::Root Bio::Align::StatisticsI );

=head2 number_of_comparable_bases

 Title   : number_of_comparable_bases
 Usage   : my $bases = $stat->number_of_comparable_bases($aln);
 Function: Returns the count of the number of bases that can be
           compared (L) in this alignment ( length - gaps)
 Returns : integer
 Args    : Bio::Align::AlignI


=cut

sub number_of_comparable_bases{
   my ($self,$aln) = @_;
   if( ! defined $aln || ! $aln->isa('Bio::Align::AlignI') ) {
       $self->warn("Must provide a Bio::Align::AlignI compliant object to Bio::Align::PairwiseStatistics");
       return 0;
   } elsif( $aln->no_sequences != 2 ) { 
       $self->warn("only pairwise calculations currently supported");
   }
   my $L = $aln->length - $self->number_of_gaps($aln);
   return $L;
}

=head2 number_of_differences

 Title   : number_of_differences
 Usage   : my $nd = $stat->number_of_distances($aln);
 Function: Returns the number of differences between two  
 Returns : integer
 Args    : Bio::Align::AlignI


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
 Function: Returns the number of differences between two  
 Example :
 Returns : 
 Args    :


=cut

sub number_of_gaps{
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
   my $gapcount = 0;
    for (my $i = 0;$i<$aln->length; $i++ ) { 
	($gapcount++) && next if( $firstseq->[$i]  =~ /^$GapChars$/);
	foreach my $seq ( @seqs ) {
	    ($gapcount++) && next if( $seq->[$i]  =~ /^$GapChars$/);
	}
    }
    return $gapcount;
}

1;
