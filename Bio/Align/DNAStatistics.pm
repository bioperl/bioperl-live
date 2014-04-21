#
# BioPerl module for Bio::Align::DNAStatistics
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-AT-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Align::DNAStatistics - Calculate some statistics for a DNA alignment

=head1 SYNOPSIS

  use Bio::AlignIO;
  use Bio::Align::DNAStatistics;

  my $stats = Bio::Align::DNAStatistics->new();
  my $alignin = Bio::AlignIO->new(-format => 'emboss',
                                 -file   => 't/data/insulin.water');
  my $aln = $alignin->next_aln;
  my $jcmatrix = $stats->distance(-align => $aln, 
                                  -method => 'Jukes-Cantor');

  print $jcmatrix->print_matrix;
  ## and for measurements of synonymous /nonsynonymous substitutions ##

  my $in = Bio::AlignIO->new(-format => 'fasta',
                            -file   => 't/data/nei_gojobori_test.aln');
  my $alnobj = $in->next_aln;
  my ($seq1id,$seq2id) = map { $_->display_id } $alnobj->each_seq;
  my $results = $stats->calc_KaKs_pair($alnobj, $seq1id, $seq2id);
  print "comparing ".$results->[0]{'Seq1'}." and ".$results->[0]{'Seq2'}."\n";
  for (sort keys %{$results->[0]} ){
      next if /Seq/;
      printf("%-9s %.4f \n",$_ , $results->[0]{$_});
  }

  my $results2 = $stats->calc_all_KaKs_pairs($alnobj);
  for my $an (@$results2){
      print "comparing ". $an->{'Seq1'}." and ". $an->{'Seq2'}. " \n";
      for (sort keys %$an ){
	  next if /Seq/;
	  printf("%-9s %.4f \n",$_ , $an->{$_});
      }
      print "\n\n";
  }

  my $result3 = $stats->calc_average_KaKs($alnobj, 1000);
  for (sort keys %$result3 ){
      next if /Seq/;
      printf("%-9s %.4f \n",$_ , $result3->{$_});
  }

=head1 DESCRIPTION

This object contains routines for calculating various statistics and
distances for DNA alignments.  The routines are not well tested and do
contain errors at this point.  Work is underway to correct them, but
do not expect this code to give you the right answer currently!  Use
dnadist/distmat in the PHLYIP or EMBOSS packages to calculate the
distances.


Several different distance method calculations are supported.  Listed
in brackets are the pattern which will match

=over 3

=item *

JukesCantor [jc|jukes|jukescantor|jukes-cantor]

=item *

Uncorrected [jcuncor|uncorrected]

=item *

F81 [f81|felsenstein]

=item *

Kimura [k2|k2p|k80|kimura]

=item *

Tamura [t92|tamura|tamura92]

=item *

F84 [f84|felsenstein84]

=item *

TajimaNei [tajimanei|tajima\-nei]

=item *

JinNei [jinnei|jin\-nei] (not implemented)

=back

There are also three methods to calculate the ratio of synonymous to
non-synonymous mutations.  All are implementations of the Nei-Gojobori
evolutionary pathway method and use the Jukes-Cantor method of
nucleotide substitution. This method works well so long as the
nucleotide frequencies are roughly equal and there is no significant
transition/transversion bias.  In order to use these methods there are
several pre-requisites for the alignment.

=over 3

=item 1

DNA alignment must be based on protein alignment. Use the subroutine
L<Bio::Align::Utilities/aa_to_dna_aln> to achieve this.

=item 2

Therefore alignment gaps must be in multiples of 3 (representing an aa
deletion/insertion) and at present must be indicated by a '-' symbol.

=item 3

Alignment must be solely of coding region and be in reading frame 0 to
achieve meaningful results

=item 4

Alignment must therefore be a multiple of 3 nucleotides long.

=item 5

All sequences must be the same length (including gaps). This should be
the case anyway if the sequences have been automatically aligned using
a program like Clustal.

=item 6

Only the standard codon alphabet is supported at present.

=back

calc_KaKs_pair() calculates a number of statistics for a named pair of
sequences in the alignment.

calc_all_KaKs_pairs() calculates these statistics for all pairwise
comparisons in an MSA.  The statistics returned are:

=over 3

=item *

S_d - Number of synonymous mutations between the 2 sequences.

=item *

N_d - Number of non-synonymous mutations between the 2 sequences.

=item *

S -  Mean number of  synonymous sites in both sequences.

=item *

N -  mean number of  synonymous sites in both sequences.

=item *

P_s - proportion of synonymous differences in both sequences given by
P_s = S_d/S.

=item *

P_n - proportion of non-synonymous differences in both sequences given
by P_n = S_n/S.

=item *

D_s - estimation of synonymous mutations per synonymous site (by
Jukes-Cantor).

=item *

D_n - estimation of non-synonymous mutations per non-synonymous site (by
Jukes-Cantor).

=item *

D_n_var - estimation of variance of D_n .

=item *

D_s_var - estimation of variance of S_n.

=item *

z_value - calculation of z value.Positive value indicates D_n E<gt> D_s,
negative value indicates D_s E<gt> D_n.

=back

The statistics returned by calc_average_KaKs are:

=over 3

=item *

D_s - Average number of synonymous mutations/synonymous site.

=item *

D_n - Average number of non-synonymous mutations/non-synonymous site.

=item *

D_s_var - Estimated variance of Ds from bootstrapped alignments.

=item *

D_n_var - Estimated variance of Dn from bootstrapped alignments.

=item *

z_score - calculation of z value. Positive value indicates D_n E<gt>D_s,
negative values vice versa.

=back

The design of the code is based around the explanation of the
Nei-Gojobori algorithm in the excellent book "Molecular Evolution and
Phylogenetics" by Nei and Kumar, published by Oxford University
Press. The methods have been tested using the worked example 4.1 in
the book, and reproduce those results. If people like having this sort
of analysis in BioPerl other methods for estimating Ds and Dn can be
provided later.

Much of the DNA distance code is based on implementations in EMBOSS
(Rice et al, www.emboss.org) [distmat.c] and PHYLIP (J. Felsenstein et
al) [dnadist.c].  Insight also gained from Eddy, Durbin, Krogh, &
Mitchison.

=head1 REFERENCES

=over 3

=item *

D_JukesCantor 

"Phylogenetic Inference", Swoffrod, Olsen, Waddell and Hillis, in
Mol. Systematics, 2nd ed, 1996, Ch 11.  Derived from "Evolution of
Protein Molecules", Jukes & Cantor, in Mammalian Prot. Metab., III,
1969, pp. 21-132.

=item *

D_Tamura

K Tamura, Mol. Biol. Evol. 1992, 9, 678.

=item *

D_Kimura 

M Kimura, J. Mol. Evol., 1980, 16, 111.

=item *

JinNei 

Jin and Nei, Mol. Biol. Evol. 82, 7, 1990.

=item *

D_TajimaNei

Tajima and Nei, Mol. Biol. Evol. 1984, 1, 269.

=back

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

Email jason-AT-bioperl.org

=head1 CONTRIBUTORS

Richard Adams, richard.adams@ed.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Align::DNAStatistics;
use vars qw(%DNAChanges @Nucleotides %NucleotideIndexes
	    $GapChars $SeqCount $DefaultGapPenalty %DistanceMethods
            $CODONS %synchanges $synsites $Precision $GCChhars);
use strict;
use Bio::Align::PairwiseStatistics;
use Bio::Matrix::PhylipDist;
use Bio::Tools::IUPAC;

BEGIN {
    $GapChars = '[\.\-]';
    $GCChhars = '[GCS]';
    @Nucleotides = qw(A G T C);
    $SeqCount = 2;
    $Precision = 5;
    
    # these values come from EMBOSS distmat implementation
    %NucleotideIndexes = ( 'A' => 0,
			   'T' => 1,
			   'C' => 2,
			   'G' => 3,

			   'AT' => 0,
			   'AC' => 1,
			   'AG' => 2,
			   'CT' => 3,
			   'GT' => 4,
			   'CG' => 5,

# these are wrong now
#			   'S' => [ 1, 3],
#			   'W' => [ 0, 4],
#			   'Y' => [ 2, 3],
#			   'R' => [ 0, 1],
#			   'M' => [ 0, 3],
#			   'K' => [ 1, 2],
#			   'B' => [ 1, 2, 3],
#			   'H' => [ 0, 2, 3],
#			   'V' => [ 0, 1, 3],
#			   'D' => [ 0, 1, 2],
			   );

    $DefaultGapPenalty = 0;
    # could put ambiguities here?
    %DNAChanges = ( 'Transversions' => { 'A' => [ 'T', 'C'],
					 'T' => [ 'A', 'G'],
					 'C' => [ 'A', 'G'],
					 'G' => [ 'C', 'T'],
				     },
		    'Transitions'   => { 'A' => [ 'G' ],
					 'G' => [ 'A' ],
					 'C' => [ 'T' ],
					 'T' => [ 'C' ],
				     },
		    );
    %DistanceMethods = ( 'jc|jukes|jukescantor|jukes\-cantor' => 'JukesCantor',
			 'jcuncor|uncorrected'   => 'Uncorrected',
			 'f81|felsenstein81'     => 'F81',
			 'k2|k2p|k80|kimura'     => 'Kimura',
			 't92|tamura|tamura92'   => 'Tamura',
			 'f84|felsenstein84'     => 'F84',
			 'tajimanei|tajima\-nei' => 'TajimaNei',
			 'jinnei|jin\-nei'       => 'JinNei');

}
use base qw(Bio::Root::Root Bio::Align::StatisticsI);

## generate look up hashes for Nei_Gojobori methods##
$CODONS = get_codons();
my @t = split '', "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
#create look up hash of number of possible synonymous mutations per codon
$synsites = get_syn_sites();
#create reference look up hash of single basechanges in codons
%synchanges = get_syn_changes();



=head2 new

 Title   : new
 Usage   : my $obj = Bio::Align::DNAStatistics->new();
 Function: Builds a new Bio::Align::DNAStatistics object 
 Returns : Bio::Align::DNAStatistics
 Args    : none


=cut

sub new { 
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    
    $self->pairwise_stats( Bio::Align::PairwiseStatistics->new());

    return $self;
}


=head2 distance

 Title   : distance
 Usage   : my $distance_mat = $stats->distance(-align  => $aln, 
		 			       -method => $method);
 Function: Calculates a distance matrix for all pairwise distances of
           sequences in an alignment.
 Returns : L<Bio::Matrix::PhylipDist> object
 Args    : -align  => Bio::Align::AlignI object
           -method => String specifying specific distance method 
                      (implementing class may assume a default)
See also: L<Bio::Matrix::PhylipDist>

=cut

sub distance{
   my ($self,@args) = @_;
   my ($aln,$method) = $self->_rearrange([qw(ALIGN METHOD)],@args);
   if( ! defined $aln || ! ref ($aln) || ! $aln->isa('Bio::Align::AlignI') ) { 
       $self->throw("Must supply a valid Bio::Align::AlignI for the -align parameter in distance");
   }
   $method ||= 'JukesCantor';
   foreach my $m ( keys %DistanceMethods ) {
       if(defined $m &&  $method =~ /$m/i ) {
	   my $mtd = "D_$DistanceMethods{$m}";
	   return $self->$mtd($aln);
       }
   }
   $self->warn("Unrecognized distance method $method must be one of [".
	       join(',',$self->available_distance_methods())."]");
   return;
}

=head2 available_distance_methods

 Title   : available_distance_methods
 Usage   : my @methods = $stats->available_distance_methods();
 Function: Enumerates the possible distance methods
 Returns : Array of strings
 Args    : none


=cut

sub available_distance_methods{
   my ($self,@args) = @_;
   return values %DistanceMethods;
}

=head2 D - distance methods


=cut


=head2 D_JukesCantor

 Title   : D_JukesCantor
 Usage   : my $d = $stat->D_JukesCantor($aln)
 Function: Calculates D (pairwise distance) between 2 sequences in an 
           alignment using the Jukes-Cantor 1 parameter model. 
 Returns : L<Bio::Matrix::PhylipDist>
 Args    : L<Bio::Align::AlignI> of DNA sequences
           double - gap penalty


=cut

sub D_JukesCantor{
   my ($self,$aln,$gappenalty) = @_;
   return 0 unless $self->_check_arg($aln);
   $gappenalty = $DefaultGapPenalty unless defined $gappenalty;
   # ambiguities ignored at this point
   my (@seqs,@names,@values,%dist);
   my $seqct = 0;
   foreach my $seq ( $aln->each_seq) {
       push @names, $seq->display_id;
       push @seqs, uc $seq->seq();
       $seqct++;
   }
   my $precisionstr = "%.$Precision"."f";
   for(my $i = 0; $i < $seqct-1; $i++ ) {
       # (diagonals) distance is 0 for same sequence
       $dist{$names[$i]}->{$names[$i]} = [$i,$i];
       $values[$i][$i] = sprintf($precisionstr,0);        

       for( my $j = $i+1; $j < $seqct; $j++ ) {
	   my ($matrix,$pfreq,$gaps) = $self->_build_nt_matrix($seqs[$i],
							       $seqs[$j]);
	   # just want diagonals
	   my $m = ( $matrix->[0]->[0] + $matrix->[1]->[1] + 
		     $matrix->[2]->[2] + $matrix->[3]->[3] );
	   my $D = 1 - ( $m / ($aln->length - $gaps + ( $gaps * $gappenalty)));
	   my $d = (- 3 / 4) * log ( 1 - (4 * $D/ 3));
	   # fwd and rev lookup
	   $dist{$names[$i]}->{$names[$j]} = [$i,$j];
	   $dist{$names[$j]}->{$names[$i]} = [$i,$j];	   
	   $values[$j][$i] = $values[$i][$j] = sprintf($precisionstr,$d);
           # (diagonals) distance is 0 for same sequence
	   $dist{$names[$j]}->{$names[$j]} = [$j,$j];   
	   $values[$j][$j] = sprintf($precisionstr,0);
       }
   }
   return Bio::Matrix::PhylipDist->new(-program => 'bioperl_DNAstats',
				       -matrix  => \%dist,
				       -names   => \@names,
				       -values  => \@values);   
}

=head2 D_F81

 Title   : D_F81
 Usage   : my $d = $stat->D_F81($aln)
 Function: Calculates D (pairwise distance) between 2 sequences in an 
           alignment using the Felsenstein 1981 distance model. 
           Relaxes the assumption of equal base frequencies that is
           in JC.
 Returns : L<Bio::Matrix::PhylipDist>
 Args    : L<Bio::Align::AlignI> of DNA sequences


=cut

sub D_F81{
   my ($self,$aln,$gappenalty) = @_;
   return 0 unless $self->_check_arg($aln);
   $gappenalty = $DefaultGapPenalty unless defined $gappenalty;
   # ambiguities ignored at this point
   my (@seqs,@names,@values,%dist);
   my $seqct = 0;
   foreach my $seq ( $aln->each_seq) {
       push @names, $seq->display_id;;
       push @seqs, uc $seq->seq();
       $seqct++;
   }
   my $precisionstr = "%.$Precision"."f";
   for(my $i = 0; $i < $seqct-1; $i++ ) {
       # (diagonals) distance is 0 for same sequence
       $dist{$names[$i]}->{$names[$i]} = [$i,$i];
       $values[$i][$i] = sprintf($precisionstr,0);        

       for( my $j = $i+1; $j < $seqct; $j++ ) {
	   
	   my ($matrix,$pfreq,$gaps) = $self->_build_nt_matrix($seqs[$i],
							       $seqs[$j]);
	   # just want diagonals
	   my $m = ( $matrix->[0]->[0] + $matrix->[1]->[1] + 
		     $matrix->[2]->[2] + $matrix->[3]->[3] );
	   my $D = 1 - ( $m / ($aln->length - $gaps + ( $gaps * $gappenalty)));
	   my $d = (- 3 / 4) * log ( 1 - (4 * $D/ 3));
	   # fwd and rev lookup
	   $dist{$names[$i]}->{$names[$j]} = [$i,$j];
	   $dist{$names[$j]}->{$names[$i]} = [$i,$j];	   
	   $values[$j][$i] = $values[$i][$j] = sprintf($precisionstr,$d);
           # (diagonals) distance is 0 for same sequence
	   $dist{$names[$j]}->{$names[$j]} = [$j,$j];	   
	   $values[$j][$j] = sprintf($precisionstr,0); 
       }
   }
   return Bio::Matrix::PhylipDist->new(-program => 'bioperl_DNAstats',
				       -matrix  => \%dist,
				       -names   => \@names,
				       -values  => \@values);   
}

=head2 D_Uncorrected

 Title   : D_Uncorrected
 Usage   : my $d = $stats->D_Uncorrected($aln)
 Function: Calculate a distance D, no correction for multiple substitutions 
           is used.  In rare cases where sequences may not overlap, 'NA' is
           substituted for the distance.
 Returns : L<Bio::Matrix::PhylipDist>
 Args    : L<Bio::Align::AlignI> (DNA Alignment)
           [optional] gap penalty

=cut

sub D_Uncorrected {
   my ($self,$aln,$gappenalty) = @_;
   $gappenalty = $DefaultGapPenalty unless defined $gappenalty;
   return 0 unless $self->_check_arg($aln);
   # ambiguities ignored at this point
   my (@seqs,@names,@values,%dist);
   my $seqct = 0;
   foreach my $seq ( $aln->each_seq) {
       push @names, $seq->display_id;
       push @seqs, uc $seq->seq();
       $seqct++;
   }

   my $precisionstr = "%.$Precision"."f";
   my $len = $aln->length;
   for( my $i = 0; $i < $seqct-1; $i++ ) {
       # (diagonals) distance is 0 for same sequence
       $dist{$names[$i]}->{$names[$i]} = [$i,$i];
       $values[$i][$i] = sprintf($precisionstr,0);
       
       for( my $j = $i+1; $j < $seqct; $j++ ) {
	   my ($matrix,$pfreq,$gaps) = $self->_build_nt_matrix($seqs[$i],
							       $seqs[$j]);
	   my $m = ( $matrix->[0]->[0] + 
		     $matrix->[1]->[1] +
		     $matrix->[2]->[2] +
		     $matrix->[3]->[3] );
       my $denom = ( $len - $gaps + ( $gaps * $gappenalty));
       
       $self->warn("No distance calculated between $names[$i] and $names[$j], inserting -1")
            unless $denom;
       
	   my $D = $denom ? 1 - ( $m / $denom) : -1;
	   # fwd and rev lookup
	   $dist{$names[$i]}->{$names[$j]} = [$i,$j];
	   $dist{$names[$j]}->{$names[$i]} = [$i,$j];
	   $values[$j][$i] = $values[$i][$j] = $denom ? sprintf($precisionstr,$D)
                                                  : sprintf("%-*s", $Precision + 2, $D);
           # (diagonals) distance is 0 for same sequence
	   $dist{$names[$j]}->{$names[$j]} = [$j,$j];	   
	   $values[$j][$j] = sprintf($precisionstr,0); 
       }
   }
   return Bio::Matrix::PhylipDist->new(-program => 'bioperl_DNAstats',
				       -matrix  => \%dist,
				       -names   => \@names,
				       -values  => \@values); 
}


# M Kimura, J. Mol. Evol., 1980, 16, 111.

=head2 D_Kimura

 Title   : D_Kimura
 Usage   : my $d = $stat->D_Kimura($aln)
 Function: Calculates D (pairwise distance) between all pairs of sequences 
           in an alignment using the Kimura 2 parameter model.
 Returns : L<Bio::Matrix::PhylipDist>
 Args    : L<Bio::Align::AlignI> of DNA sequences


=cut

sub D_Kimura {
   my ($self,$aln) = @_;
   return 0 unless $self->_check_arg($aln);
   # ambiguities ignored at this point
   my (@names,@values,%dist);
   my $seqct = 0;
   foreach my $seq ( $aln->each_seq) {
       push @names, $seq->display_id;
       $seqct++;
   }

   my $precisionstr = "%.$Precision"."f";

   for( my $i = 0; $i < $seqct-1; $i++ ) {
       # (diagonals) distance is 0 for same sequence
       $dist{$names[$i]}->{$names[$i]} = [$i,$i];
       $values[$i][$i] = sprintf($precisionstr,0);

       for( my $j = $i+1; $j < $seqct; $j++ ) {
	   my $pairwise = $aln->select_noncont($i+1,$j+1);
	   my $L = $self->pairwise_stats->number_of_comparable_bases($pairwise);
	   unless( $L ) { 
	       $L = 1;
	   }
	   my $P = $self->transitions($pairwise) / $L;
	   my $Q = $self->transversions($pairwise) / $L;
	   my $K = 0;
	   my $denom = ( 1 - (2 * $P) - $Q);
	   if( $denom == 0 ) {
	       $self->throw("cannot find distance for ",$i+1,
			    ",",$j+1," $P, $Q\n");
	   }
	   my $a = 1 / ( 1 - (2 * $P) - $Q);
	   my $b = 1 / ( 1 - 2 * $Q );
	   if( $a < 0 || $b < 0 ) { 
	       $K = -1;
	   } else{ 
	       $K = (1/2) * log ( $a ) + (1/4) * log($b);
	   }
	   # fwd and rev lookup
	   $dist{$names[$i]}->{$names[$j]} = [$i,$j];
	   $dist{$names[$j]}->{$names[$i]} = [$i,$j];	   
	   $values[$j][$i] = $values[$i][$j] = sprintf($precisionstr,$K);
           # (diagonals) distance is 0 for same sequence
	   $dist{$names[$j]}->{$names[$j]} = [$j,$j];	   
	   $values[$j][$j] = sprintf($precisionstr,0); 
       }
   }
   return Bio::Matrix::PhylipDist->new(-program => 'bioperl_DNAstats',
				       -matrix  => \%dist,
				       -names   => \@names,
				       -values  => \@values); 
}


=head2 D_Kimura_variance

 Title   : D_Kimura
 Usage   : my $d = $stat->D_Kimura_variance($aln)
 Function: Calculates D (pairwise distance) between all pairs of sequences 
           in an alignment using the Kimura 2 parameter model.
 Returns : array of 2 L<Bio::Matrix::PhylipDist>,
           the first is the Kimura distance and the second is
           a matrix of variance V(K)
 Args    : L<Bio::Align::AlignI> of DNA sequences


=cut

sub D_Kimura_variance {
   my ($self,$aln) = @_;
   return 0 unless $self->_check_arg($aln);
   # ambiguities ignored at this point
   my (@names,@values,%dist,@var);
   my $seqct = 0;
   foreach my $seq ( $aln->each_seq) {
       push @names, $seq->display_id;
       $seqct++;
   }

   my $precisionstr = "%.$Precision"."f";

   for( my $i = 0; $i < $seqct-1; $i++ ) {
       # (diagonals) distance is 0 for same sequence
       $dist{$names[$i]}->{$names[$i]} = [$i,$i];
       $values[$i][$i] = sprintf($precisionstr,0);

       for( my $j = $i+1; $j < $seqct; $j++ ) {
	   my $pairwise = $aln->select_noncont($i+1,$j+1);
	   my $L = $self->pairwise_stats->number_of_comparable_bases($pairwise);
	   unless( $L ) { 
	       $L = 1;
	   }
	   my $P = $self->transitions($pairwise) / $L;
	   my $Q = $self->transversions($pairwise) / $L;
	   my ($a,$b,$K,$var_k);
	   my $a_denom = ( 1 - (2 * $P) - $Q);
	   my $b_denom = 1 - 2 * $Q;
	   unless( $a_denom > 0 && $b_denom > 0 ) {
	       $a = 1;
	       $b = 1;
	       $K = -1;
	       $var_k = -1;
	   } else { 
	       $a = 1 / $a_denom;
	       $b = 1 / $b_denom;
	       $K = (1/2) * log ( $a ) + (1/4) * log($b);
	       # from Wu and Li 1985 which in turn is from Kimura 1980
	       my $c = ( $a - $b ) / 2;
	       my $d = ( $a + $b ) / 2;
	       $var_k = ( $a**2 * $P + $d**2 * $Q - ( $a * $P + $d * $Q)**2 ) / $L;
	   }

	   # fwd and rev lookup
	   $dist{$names[$i]}->{$names[$j]} = [$i,$j];
	   $dist{$names[$j]}->{$names[$i]} = [$i,$j];	   
	   $values[$j][$i] = $values[$i][$j] = sprintf($precisionstr,$K);
           # (diagonals) distance is 0 for same sequence
	   $dist{$names[$j]}->{$names[$j]} = [$j,$j];   
	   $values[$j]->[$j] = sprintf($precisionstr,0); 
	   
	   $var[$j]->[$i] = $var[$i]->[$j] = sprintf($precisionstr,$var_k);
	   $var[$j]->[$j] = $values[$j]->[$j];
       }
   }
   return ( Bio::Matrix::PhylipDist->new(-program => 'bioperl_DNAstats',
					 -matrix  => \%dist,
					 -names   => \@names,
					 -values  => \@values),
	    Bio::Matrix::PhylipDist->new(-program => 'bioperl_DNAstats',
					 -matrix  => \%dist,
					 -names   => \@names,
					 -values  => \@var)
	    );
}


#  K Tamura, Mol. Biol. Evol. 1992, 9, 678.

=head2 D_Tamura

 Title   : D_Tamura
 Usage   : Calculates D (pairwise distance) between 2 sequences in an 
           alignment using Tamura 1992 distance model. 
 Returns : L<Bio::Matrix::PhylipDist>
 Args    : L<Bio::Align::AlignI> of DNA sequences


=cut

sub D_Tamura {
   my ($self,$aln) = @_;
   return 0 unless $self->_check_arg($aln);
   # ambiguities ignored at this point
   my (@seqs,@names,@values,%dist,$i,$j);
   my $seqct = 0;
   my $length = $aln->length;
   foreach my $seq ( $aln->each_seq) {
       push @names, $seq->display_id;;
       push @seqs, uc $seq->seq();
       $seqct++;
   }

   my $precisionstr = "%.$Precision"."f";
   my (@gap,@gc,@trans,@tranv,@score);
   $i = 0;
   for my $t1 ( @seqs ) {
       $j = 0;
       for my $t2 ( @seqs ) {
	   $gap[$i][$j] = 0;
	   for( my $k = 0; $k < $length; $k++ ) {
	       my ($c1,$c2) = ( substr($seqs[$i],$k,1),
				substr($seqs[$j],$k,1) );
	       if( $c1 =~ /^$GapChars$/ ||
		   $c2 =~ /^$GapChars$/ ) {
		   $gap[$i][$j]++;	
	       } elsif( $c2 =~ /^$GCChhars$/i ) {
		   $gc[$i][$j]++;
	       } 
	   }
	   $gc[$i][$j] = ( $gc[$i][$j] / 
			   ($length - $gap[$i][$j]) );
	   $j++;
       }
       $i++;
   }
   
   for( $i = 0; $i < $seqct-1; $i++ ) {
       # (diagonals) distance is 0 for same sequence
       $dist{$names[$i]}->{$names[$i]} = [$i,$i];
       $values[$i][$i] = sprintf($precisionstr,0);
       
       for( $j = $i+1; $j < $seqct; $j++ ) {
	   
	   my $pairwise = $aln->select_noncont($i+1,$j+1);
	   my $L = $self->pairwise_stats->number_of_comparable_bases($pairwise);
	   my $P = $self->transitions($pairwise) / $L;
	   my $Q = $self->transversions($pairwise) / $L;
	   my $C = $gc[$i][$j] + $gc[$j][$i]- 
	       ( 2 * $gc[$i][$j] * $gc[$j][$i] );
	   if( $P ) {
	       $P = $P / $C;
	   }
	   my $d = -($C * log(1- $P - $Q)) -(0.5* ( 1 - $C) * log(1 - 2 * $Q));
           # fwd and rev lookup
	   $dist{$names[$i]}->{$names[$j]} = [$i,$j];
	   $dist{$names[$j]}->{$names[$i]} = [$i,$j];	   
	   $values[$j][$i] = $values[$i][$j] = sprintf($precisionstr,$d);
           # (diagonals) distance is 0 for same sequence
	   $dist{$names[$j]}->{$names[$j]} = [$j,$j];
	   $values[$j][$j] = sprintf($precisionstr,0); 
       }
   }
   return Bio::Matrix::PhylipDist->new(-program => 'bioperl_DNAstats',
				       -matrix  => \%dist,
				       -names   => \@names,
				       -values  => \@values); 

}

=head2 D_F84

 Title   : D_F84
 Usage   : my $d = $stat->D_F84($aln)
 Function: Calculates D (pairwise distance) between 2 sequences in an 
           alignment using the Felsenstein 1984 distance model. 
 Returns : L<Bio::Matrix::PhylipDist>
 Args    : L<Bio::Align::AlignI> of DNA sequences
           [optional] double - gap penalty

=cut

sub D_F84 {
   my ($self,$aln,$gappenalty) = @_;
   return 0 unless $self->_check_arg($aln);
   $self->throw_not_implemented();
   # ambiguities ignored at this point
   my (@seqs,@names,@values,%dist);
   my $seqct = 0;
   foreach my $seq ( $aln->each_seq) {
       # if there is no name, 
       my $id = $seq->display_id;
       if( ! length($id) ||       # deal with empty names
	   $id =~ /^\s+$/ ) {
	   $id = $seqct+1;
       }
       push @names, $id;
       push @seqs, uc $seq->seq();
       $seqct++;
   }

   my $precisionstr = "%.$Precision"."f";

   for( my $i = 0; $i < $seqct-1; $i++ ) {
       # (diagonals) distance is 0 for same sequence
       $dist{$names[$i]}->{$names[$i]} = [$i,$i];
       $values[$i][$i] = sprintf($precisionstr,0);

       for( my $j = $i+1; $j < $seqct; $j++ ) {
       }
   }   
}

# Tajima and Nei, Mol. Biol. Evol. 1984, 1, 269.
#  Tajima-Nei correction used for multiple substitutions in the calc
# of the distance matrix. Nucleic acids only.
#
#  D = p-distance = 1 - (matches/(posns_scored + gaps)
#
#  distance = -b * ln(1-D/b)
#

=head2 D_TajimaNei

 Title   : D_TajimaNei
 Usage   : my $d = $stat->D_TajimaNei($aln)
 Function: Calculates D (pairwise distance) between 2 sequences in an 
           alignment using the TajimaNei 1984 distance model. 
 Returns : L<Bio::Matrix::PhylipDist>
 Args    : Bio::Align::AlignI of DNA sequences


=cut

sub D_TajimaNei{
   my ($self,$aln) = @_;
   return 0 unless $self->_check_arg($aln);
   # ambiguities ignored at this point
   my (@seqs,@names,@values,%dist);
   my $seqct = 0;
   foreach my $seq ( $aln->each_seq) {
       # if there is no name, 
       push @names, $seq->display_id;
       push @seqs, uc $seq->seq();
       $seqct++;
   }
   my $precisionstr = "%.$Precision"."f";
   my ($i,$j,$bs);
   # pairwise
   for( $i =0; $i < $seqct -1; $i++ ) {
       $dist{$names[$i]}->{$names[$i]} = [$i,$i];
       $values[$i][$i] = sprintf($precisionstr,0);

       for ( $j = $i+1; $j <$seqct;$j++ ) {
	   my ($matrix,$pfreq,$gaps) = $self->_build_nt_matrix($seqs[$i],
							       $seqs[$j]);
	   my $pairwise = $aln->select_noncont($i+1,$j+1);
	   my $slen = $self->pairwise_stats->number_of_comparable_bases($pairwise);	    
	   my $fij2 = 0;
	   for( $bs = 0; $bs < 4; $bs++ ) {
	       my $fi = 0;
	       map {$fi += $matrix->[$bs]->[$_] } 0..3;
	       my $fj = 0;
	       # summation 
	       map { $fj += $matrix->[$_]->[$bs] } 0..3;
	       my $fij = ( $fi && $fj ) ? ($fi + $fj) /( 2 * $slen) : 0;
	       $fij2 += $fij**2;
	   }
	   
	   my ($pair,$h) = (0,0);
	   for( $bs = 0; $bs < 3; $bs++ ) {
	       for(my $bs1 = $bs+1; $bs1 <= 3; $bs1++ ) {
		   my $fij = $pfreq->[$pair++] / $slen;
		   if( $fij ) {
		       
		       my ($ci1,$ci2,$cj1,$cj2) = (0,0,0,0);

		       map { $ci1 += $matrix->[$_]->[$bs] } 0..3;
		       map { $cj1 += $matrix->[$bs]->[$_] } 0..3;
		       map { $ci2 += $matrix->[$_]->[$bs1] } 0..3;
		       map { $cj2 += $matrix->[$bs1]->[$_] } 0..3;
		       
		       if( $fij ) {
			   $h += ( ($fij**2) / 2 ) / 
			       (  ( ( $ci1 + $cj1 ) / (2 * $slen) ) *
				  ( ( $ci2 + $cj2 ) / (2 * $slen) ) 
				  );
		       }
		       $self->debug( "slen is $slen h is $h fij = $fij ci1 =$ci1 cj1=$cj1 ci2=$ci2 cj2=$cj2\n");
		   }
	       }
	   }
	   # just want diagonals which are matches (A matched A, C -> C)

	   my $m = ( $matrix->[0]->[0] + $matrix->[1]->[1] + 
		     $matrix->[2]->[2] + $matrix->[3]->[3] );
	   my $D = 1 - ( $m / $slen);
	   my $d;
	   if( $h == 0 ) {
	       $d = -1;
	   } else {
	       my $b = (1 - $fij2 + (($D**2)/$h)) / 2;
	       my $c = 1- $D/ $b;

	       if( $c < 0 ) {
		   $d = -1;
	       } else { 
		   $d = (-1 * $b) * log ( $c);
	       }
	   }
	   # fwd and rev lookup
	   $dist{$names[$i]}->{$names[$j]} = [$i,$j];
	   $dist{$names[$j]}->{$names[$i]} = [$i,$j];	   
	   $values[$j][$i] = $values[$i][$j] = sprintf($precisionstr,$d);

           # (diagonals) distance is 0 for same sequence
	   $dist{$names[$j]}->{$names[$j]} = [$j,$j];	   
	   $values[$j][$j] = sprintf($precisionstr,0); 
       }
   }
   return Bio::Matrix::PhylipDist->new(-program => 'bioperl_DNAstats',
				       -matrix  => \%dist,
				       -names   => \@names,
				       -values  => \@values); 

}

# Jin and Nei, Mol. Biol. Evol. 82, 7, 1990.

=head2 D_JinNei

 Title   : D_JinNei
 Usage   : my $d = $stat->D_JinNei($aln)
 Function: Calculates D (pairwise distance) between 2 sequences in an 
           alignment using the Jin-Nei 1990 distance model. 
 Returns : L<Bio::Matrix::PhylipDist>
 Args    : L<Bio::Align::AlignI> of DNA sequences


=cut

sub D_JinNei{
   my ($self,@args) = @_;
   $self->warn("JinNei implementation not completed");
   return;
}

=head2 transversions

 Title   : transversions
 Usage   : my $transversions = $stats->transversion($aln);
 Function: Calculates the number of transversions between two sequences in 
           an alignment
 Returns : integer
 Args    : Bio::Align::AlignI


=cut

sub transversions{
   my ($self,$aln) = @_;
   return $self->_trans_count_helper($aln, $DNAChanges{'Transversions'});
}

=head2 transitions

 Title   : transitions
 Usage   : my $transitions = Bio::Align::DNAStatistics->transitions($aln);
 Function: Calculates the number of transitions in a given DNA alignment
 Returns : integer representing the number of transitions
 Args    : Bio::Align::AlignI object


=cut

sub transitions{
   my ($self,$aln) = @_;
   return $self->_trans_count_helper($aln, $DNAChanges{'Transitions'});
}


sub _trans_count_helper {
    my ($self,$aln,$type) = @_;
    return 0 unless( $self->_check_arg($aln) );
    if( ! $aln->is_flush ) { $self->throw("must be flush") }
    my (@tcount);
    my ($first,$second) = ( uc $aln->get_seq_by_pos(1)->seq(),
			    uc $aln->get_seq_by_pos(2)->seq() );
    my $alen = $aln->length; 
    for (my $i = 0;$i<$alen; $i++ ) { 
	my ($c1,$c2) = ( substr($first,$i,1),
			 substr($second,$i,1) );
	if( $c1 ne $c2 ) { 
	    foreach my $nt ( @{$type->{$c1}} ) {
		if( $nt eq $c2) {
		   $tcount[$i]++;
	       }
	    }
	}
    }
    my $sum = 0;
    map { if( $_) { $sum += $_} } @tcount;
    return $sum;
}

# this will generate a matrix which records across the row, the number
# of DNA subst 
# 
sub _build_nt_matrix {
    my ($self,$seqa,$seqb) = @_;
    

    my $basect_matrix = [ [ qw(0 0 0 0) ],  # number of bases that match
			  [ qw(0 0 0 0) ],
			  [ qw(0 0 0 0) ],
			  [ qw(0 0 0 0) ] ];
    my $gaps = 0;                           # number of gaps
    my $pfreq = [ qw( 0 0 0 0 0 0)];        # matrix for pair frequency
    my $len_a = length($seqa);
    for( my $i = 0; $i < $len_a; $i++) {
	my ($ti,$tj) = (substr($seqa,$i,1),substr($seqb,$i,1));
	$ti =~ tr/U/T/;
	$tj =~ tr/U/T/;

	if( $ti =~ /^$GapChars$/) { $gaps++; next; }
	if( $tj =~ /^$GapChars$/) { $gaps++; next }

	my $ti_index = $NucleotideIndexes{$ti};		
	my $tj_index = $NucleotideIndexes{$tj};	    

	if( ! defined $ti_index ) {
	    $self->warn("ti_index not defined for $ti\n");
	    next;
	}
	
	$basect_matrix->[$ti_index]->[$tj_index]++;
	
	if( $ti ne $tj ) {
	    $pfreq->[$NucleotideIndexes{join('',sort ($ti,$tj))}]++;
	}
    }
    return ($basect_matrix,$pfreq,$gaps);
}

sub _check_ambiguity_nucleotide {
    my ($base1,$base2) = @_;
    my %iub = Bio::Tools::IUPAC->iupac_iub();
    my @amb1 = @{ $iub{uc($base1)} };
    my @amb2 = @{ $iub{uc($base2)} };    
    my ($pmatch) = (0);
    for my $amb ( @amb1 ) {
	if( grep { $amb eq $_ } @amb2 ) {
	    $pmatch = 1;
	    last;
	}
    }
    if( $pmatch ) { 
	return (1 / scalar @amb1) * (1 / scalar @amb2);
    } else { 
	return 0;
    }
}


sub _check_arg {
    my($self,$aln ) = @_;
    if( ! defined $aln || ! $aln->isa('Bio::Align::AlignI') ) {
	$self->warn("Must provide a Bio::Align::AlignI compliant object to Bio::Align::DNAStatistics");
	return 0;
    } elsif( $aln->get_seq_by_pos(1)->alphabet ne 'dna' ) { 
	$self->warn("Must provide a DNA alignment to Bio::Align::DNAStatistics, you provided a " . $aln->get_seq_by_pos(1)->alphabet);
	return 0;
    }
    return 1;
}

=head2 Data Methods

=cut

=head2 pairwise_stats

 Title   : pairwise_stats
 Usage   : $obj->pairwise_stats($newval)
 Function: 
 Returns : value of pairwise_stats
 Args    : newvalue (optional)


=cut

sub pairwise_stats{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_pairwise_stats'} = $value;
    }
    return $self->{'_pairwise_stats'};

}

=head2 calc_KaKs_pair

 Title    : calc_KaKs_pair
 Useage   : my $results = $stats->calc_KaKs_pair($alnobj,
            $name1, $name2).
 Function : calculates Nei-Gojobori statistics for pairwise 
            comparison.
 Args     : A Bio::Align::AlignI compliant object such as a 
            Bio::SimpleAlign object, and 2 sequence name strings.
 Returns  : a reference to a hash of statistics with keys as 
            listed in Description.

=cut

sub calc_KaKs_pair {
    my ( $self, $aln, $seq1_id, $seq2_id) = @_;
    $self->throw("Needs 3 arguments - an alignment object, and 2 sequence ids") 
	if @_!= 4;
    $self->throw ("This calculation needs a Bio::Align::AlignI compatible object, not a [ " . ref($aln) . " ]object") unless $aln->isa('Bio::Align::AlignI');
    my @seqs = (
		#{id => $seq1_id, seq =>($aln->each_seq_with_id($seq1_id))[0]->seq},
		#{id => $seq2_id, seq =>($aln->each_seq_with_id($seq2_id))[0]->seq}
		{id => $seq1_id, seq => uc(($aln->each_seq_with_id($seq1_id))[0]->seq)},
                {id => $seq2_id, seq => uc(($aln->each_seq_with_id($seq2_id))[0]->seq)}
	       ) ;
    if (length($seqs[0]{'seq'}) != length($seqs[1]{'seq'})) {
	$self->throw(" aligned sequences must be of equal length!");
    }
    my $results = [];
    $self->_get_av_ds_dn(\@seqs, $results);
    return $results;

}

=head2 calc_all_KaKs_pairs

 Title    : calc_all_KaKs_pairs
 Useage   : my $results2 = $stats->calc_KaKs_pair($alnobj).
 Function : Calculates Nei_gojobori statistics for all pairwise
            combinations in sequence.
 Arguments: A Bio::Align::ALignI compliant object such as
            a Bio::SimpleAlign object.
 Returns  : A reference to an array of hashes of statistics of
            all pairwise comparisons in the alignment.

=cut



sub calc_all_KaKs_pairs {
#returns a multi_element_array with all pairwise comparisons
	my ($self,$aln) = @_;
	$self->throw ("This calculation needs a Bio::Align::AlignI compatible object, not a [ " . ref($aln) . " ]object") unless $aln->isa('Bio::Align::AlignI');
	my @seqs;
	for my $seq ($aln->each_seq) {
		push @seqs, {id => $seq->display_id, seq=>$seq->seq};
		}
	my $results ;
	$results = $self->_get_av_ds_dn(\@seqs, $results);
	return $results;
}

=head2 calc_average_KaKs

 Title    : calc_average_KaKs.  
 Useage   : my $res= $stats->calc_average_KaKs($alnobj, 1000).
 Function : calculates Nei_Gojobori stats for average of all 
            sequences in the alignment.
 Args     : A Bio::Align::AlignI compliant object such as a
            Bio::SimpleAlign object, number of bootstrap iterations
            (default 1000).
 Returns  : A reference to a hash of statistics as listed in Description.

=cut

sub calc_average_KaKs {
#calculates global value for sequences in alignment using bootstrapping
#this is quite slow (~10 seconds per  3 X 200nt seqs); 
    my ($self, $aln, $bootstrap_rpt) = @_;
    $bootstrap_rpt ||= 1000;
    $self->throw ("This calculation needs a Bio::Align::AlignI compatible object, not a [ " . ref($aln) . " ]object") unless $aln->isa('Bio::Align::AlignI');
    my @seqs;
    for my $seq ($aln->each_seq) {
	push @seqs, {id => $seq->display_id, seq=>$seq->seq};
    }
    my $results ;
    my ($ds_orig, $dn_orig) = $self->_get_av_ds_dn(\@seqs);
    #print "ds = $ds_orig, dn = $dn_orig\n";
    $results = {D_s => $ds_orig, D_n => $dn_orig};
    $self->_run_bootstrap(\@seqs, $results, $bootstrap_rpt);
    return $results;
}

############## primary internal subs for alignment comparisons ########################

sub _run_bootstrap {
    ### generates sampled sequences, calculates Ds and Dn values,
    ### then calculates variance of sampled sequences and add results to results hash
    ### 
    my ($self,$seq_ref, $results, $bootstrap_rpt) = @_;	
    my @seqs = @$seq_ref;
    my @btstrp_aoa; # to hold array of array of nucleotides for resampling
    my %bootstrap_values = (ds => [], dn =>[]);	# to hold list of av values 

    #1st make alternative array of codons;
    my $c = 0;
    while ($c < length $seqs[0]{'seq'}) {
	for (0..$#seqs) {
	    push @{$btstrp_aoa[$_]}, substr ($seqs[$_]{'seq'}, $c, 3);
	}
	$c+=3;
    }

    for (1..$bootstrap_rpt) {
	my $sampled = _resample (\@btstrp_aoa);
	my ($ds, $dn) = $self->_get_av_ds_dn ($sampled) ; # is array ref
	push @{$bootstrap_values{'ds'}}, $ds;
	push @{$bootstrap_values{'dn'}}, $dn;
    }	

    $results->{'D_s_var'} = sampling_variance($bootstrap_values{'ds'});
    $results->{'D_n_var'} = sampling_variance($bootstrap_values{'dn'});
    $results->{'z_score'} = 	($results->{'D_n'} - $results->{'D_s'}) / 
	sqrt($results->{'D_s_var'} + $results->{'D_n_var'} ); 
    #print "bootstrapped var_syn = 	$results->{'D_s_var'} \n" ;
    #print "bootstrapped var_nc = 	$results->{'D_n_var'} \n"; 
    #print "z is $results->{'z_score'}\n";	### end of global set up of/perm look up data
}

sub _resample {
    my $ref = shift;
    my $codon_num = scalar (@{$ref->[0]});
    my @altered;
    for (0..$codon_num -1) {	#for each codon
	my $rand = int (rand ($codon_num));
	for (0..$#$ref) {
	    push @{$altered[$_]}, $ref->[$_][$rand];
	}
    }
    my @stringed = map {join '', @$_}@altered;
    my @return;
    #now out in random name to keep other subs happy
    for (@stringed) {
	push @return, {id=>'1', seq=> $_};
    }
    return \@return;
}

sub _get_av_ds_dn {
    # takes array of hashes of sequence strings and ids   #
    my $self = shift;
    my $seq_ref = shift;
    my $result = shift if @_;
    my @caller = caller(1);
    my @seqarray = @$seq_ref;
    my $bootstrap_score_list;
    #for a multiple alignment considers all pairwise combinations#
    my %dsfor_average = (ds => [], dn => []); 
    for (my $i = 0; $i < scalar @seqarray; $i++) {
	for (my $j = $i +1; $j<scalar @seqarray; $j++ ){
#			print "comparing $i and $j\n";
	    if (length($seqarray[$i]{'seq'}) != length($seqarray[$j]{'seq'})) {
		$self->warn(" aligned sequences must be of equal length!");
		next;
	    }

	    my $syn_site_count = count_syn_sites($seqarray[$i]{'seq'}, $synsites);
	    my $syn_site_count2 = count_syn_sites($seqarray[$j]{'seq'}, $synsites);
#			print "syn 1 is $syn_site_count , syn2 is $syn_site_count2\n";
	    my ($syn_count, $non_syn_count, $gap_cnt) = analyse_mutations($seqarray[$i]{'seq'}, $seqarray[$j]{'seq'});	
	    #get averages
	    my $av_s_site = ($syn_site_count + $syn_site_count2)/2;
	    my $av_ns_syn_site = length($seqarray[$i]{'seq'}) - $gap_cnt- $av_s_site ;

	    #calculate ps and pn  (p54)
	    my $syn_prop = $syn_count / $av_s_site;
	    my $nc_prop = $non_syn_count / $av_ns_syn_site	;

	    #now use jukes/cantor to calculate D_s and D_n, would alter here if needed a different method
	    my $d_syn = $self->jk($syn_prop);
	    my $d_nc = $self->jk($nc_prop);

	    #JK calculation must succeed for continuation of calculation
	    #ret_value = -1 if error
	    next unless $d_nc >=0 && $d_syn >=0;


	    push @{$dsfor_average{'ds'}}, $d_syn;
	    push @{$dsfor_average{'dn'}}, $d_nc;

	    #if not doing bootstrap, calculate the pairwise comparisin stats
	    if ($caller[3] =~ /calc_KaKs_pair/ || $caller[3] =~ /calc_all_KaKs_pairs/) {
				#now calculate variances assuming large sample
		my $d_syn_var =  jk_var($syn_prop, length($seqarray[$i]{'seq'})  - $gap_cnt );
		my $d_nc_var =  jk_var($nc_prop, length ($seqarray[$i]{'seq'}) - $gap_cnt);
		#now calculate z_value
		#print "d_syn_var is  $d_syn_var,and d_nc_var is $d_nc_var\n";
		#my $z = ($d_nc - $d_syn) / sqrt($d_syn_var + $d_nc_var);
		my $z = ($d_syn_var + $d_nc_var) ? 
		  ($d_nc - $d_syn) / sqrt($d_syn_var + $d_nc_var) : 0;
		#	print "z is $z\n";
		push @$result , {S => $av_s_site, N=>$av_ns_syn_site,
				 S_d => $syn_count, N_d =>$non_syn_count,
				 P_s => $syn_prop, P_n=>$nc_prop,
				 D_s => @{$dsfor_average{'ds'}}[-1],
				 D_n => @{$dsfor_average{'dn'}}[-1],
				 D_n_var =>$d_nc_var, D_s_var => $d_syn_var,
				 Seq1 => $seqarray[$i]{'id'},
				 Seq2 => $seqarray[$j]{'id'},
				 z_score => $z,
			     };
		$self->warn (" number of mutations too small to justify normal test for  $seqarray[$i]{'id'} and $seqarray[$j]{'id'}\n- use Fisher's exact, or bootstrap a MSA")
		    if ($syn_count < 10 || $non_syn_count < 10 ) && $self->verbose > -1 ;
	    }#endif
	    }
    }

    #warn of failure if no results hashes are present
    #will fail if Jukes Cantor has failed for all pairwise combinations
    #$self->warn("calculation failed!") if scalar @$result ==0;

    #return results unless bootstrapping
    return $result if $caller[3]=~ /calc_all_KaKs/ || $caller[3] =~ /calc_KaKs_pair/; 
    #else if getting average for bootstrap
    return( mean ($dsfor_average{'ds'}),mean ($dsfor_average{'dn'})) ;
}


sub jk {
    my ($self, $p) = @_;
    if ($p > 0.75) {
	$self->warn( " Jukes Cantor won't  work -too divergent!");
	return -1;
    }
    return -1 * (3/4) * (log(1 - (4/3) * $p));
}

#works for large value of n (50?100?)
sub jk_var {
    my ($p, $n) = @_;
    return (9 * $p * (1 -$p))/(((3 - 4 *$p) **2) * $n);
}


# compares 2 sequences to find the number of synonymous/non
# synonymous mutations between them

sub analyse_mutations {
    my ($seq1, $seq2) = @_;
    my %mutator = ( 2=> {0=>[[1,2],  # codon positions to be altered 
			     [2,1]], # depend on which is the same
			 1=>[[0,2],
			     [2,0]],
			 2=>[[0,1],
			     [1,0]],	
		     },
		    3=> [ [0,1,2],  # all need to be altered 
			  [1,0,2],
			  [0,2,1],
			  [1,2,0],
			  [2,0,1],
			  [2,1,0] ],
		    );
    my $TOTAL   = 0;    # total synonymous changes
    my $TOTAL_n = 0;	# total non-synonymous changes
    my $gap_cnt = 0;

    my %input;
    my $seqlen = length($seq1);
    for (my $j=0; $j< $seqlen; $j+=3) {
	$input{'cod1'} = substr($seq1, $j,3);
	$input{'cod2'} = substr($seq2, $j,3);

	#ignore codon if beeing compared with gaps! 
	if ($input{'cod1'} =~ /\-/ || $input{'cod2'} =~ /\-/){
	    $gap_cnt += 3; #just increments once if there is a pair of gaps
	    next;
	}

	my ($diff_cnt, $same) = count_diffs(\%input);

	#ignore if codons are identical
	next if $diff_cnt == 0 ;
	if ($diff_cnt == 1) {
	    $TOTAL += $synchanges{$input{'cod1'}}{$input{'cod2'}};
	    $TOTAL_n += 1 - $synchanges{$input{'cod1'}}{$input{'cod2'}};
	    #print " \nfordiff is 1 , total now $TOTAL, total n now $TOTAL_n\n\n"
	}
	elsif ($diff_cnt ==2) {
	    my $s_cnt = 0;
	    my $n_cnt = 0;
	    my $tot_muts = 4;
	    #will stay 4 unless there are stop codons at intervening point
	  OUTER:for my $perm (@{$mutator{'2'}{$same}}) {
	      my $altered = $input{'cod1'};
	      my $prev= $altered;
	      #		print "$prev -> (", $t[$CODONS->{$altered}], ")";
	      for 	my $mut_i (@$perm) { #index of codon mutated
		  substr($altered, $mut_i,1) = substr($input{'cod2'}, $mut_i, 1);
		  if ($t[$CODONS->{$altered}] eq '*') {
		      $tot_muts -=2;
		      #print "changes to stop codon!!\n";
		      next OUTER;
		  }
		  else {
		      $s_cnt += $synchanges{$prev}{$altered};
		      #					print "$altered ->(", $t[$CODONS->{$altered}], ") ";
		  }
		  $prev = $altered;
	      }
	      #		print "\n";
	  }
	    if ($tot_muts != 0) {
		$TOTAL += ($s_cnt/($tot_muts/2));
		$TOTAL_n += ($tot_muts - $s_cnt)/ ($tot_muts / 2);
	    }

	}
	elsif ($diff_cnt ==3 ) {
	    my $s_cnt = 0;
	    my $n_cnt = 0;
	    my $tot_muts = 18;	#potential number  of mutations
	  OUTER: for my $perm (@{$mutator{'3'}}) {
	      my $altered = $input{'cod1'};
	      my $prev= $altered;
	      #	print "$prev -> (", $t[$CODONS->{$altered}], ")";
	      for my $mut_i (@$perm) { #index of codon mutated
		  substr($altered, $mut_i,1) = substr($input{'cod2'}, $mut_i, 1);
		  if ($t[$CODONS->{$altered}] eq '*') {
		      $tot_muts -=3;
		      #	print "changes to stop codon!!\n";
		      next OUTER;

		  }
		  else {
		      $s_cnt += $synchanges{$prev}{$altered};
		      #			print "$altered ->(", $t[$CODONS->{$altered}], ") ";
		  }
		  $prev = $altered;
	      }
	      #	print "\n";

	  }#end OUTER loop
	      #calculate number of synonymous/non synonymous mutations for that codon
	      # and add to total
	      if ($tot_muts != 0) {
		  $TOTAL += ($s_cnt / ($tot_muts /3));
		  $TOTAL_n += 3 - ($s_cnt / ($tot_muts /3));
	      }
	}			#endif $diffcnt = 3
    }				#end of sequencetraversal
    return ($TOTAL, $TOTAL_n, $gap_cnt);
}


sub count_diffs {
    #counts the number of nucleotide differences between 2 codons
    # returns this value plus the codon index of which nucleotide is the same when 2
    #nucleotides are different. This is so analyse_mutations() knows which nucleotides
    # to change.
    my $ref = shift;
    my $cnt = 0;
    my $same= undef;
    #just for 2 differences
    for (0..2) {
	if (substr($ref->{'cod1'}, $_,1) ne substr($ref->{'cod2'}, $_, 1)){
	    $cnt++;
	} else {
	    $same = $_;
	}
    }
    return ($cnt, $same);
}

=head2 get_syn_changes

 Title   : get_syn_changes
 Usage   : Bio::Align::DNAStatitics->get_syn_changes
 Function: Generate a hashref of all pairwise combinations of codns
           differing by 1
 Returns : Symetic matrix using hashes
           First key is codon
           and each codon points to a hashref of codons
           the values of which describe type of change.
           my $type = $hash{$codon1}->{$codon2};
           values are :
             1   synonymous
             0   non-syn
            -1   either codon is a stop codon
 Args    : none

=cut

sub get_syn_changes {
#hash of all pairwise combinations of codons differing by 1
# 1 = syn, 0 = non-syn, -1 = stop
    my %results;
    my @codons = _make_codons ();
    my $arr_len = scalar @codons;
    for (my $i = 0; $i < $arr_len -1; $i++) {
	my $cod1 = $codons[$i];
	for (my $j = $i +1; $j < $arr_len; $j++) {
	    my $diff_cnt = 0;
	    for my $pos(0..2) {
		$diff_cnt++ if substr($cod1, $pos, 1) ne substr($codons[$j], $pos, 1);
	    }
	    next if $diff_cnt !=1;

	    #synon change
	    if($t[$CODONS->{$cod1}] eq $t[$CODONS->{$codons[$j]}]) {
		$results{$cod1}{$codons[$j]} =1;
		$results{$codons[$j]}{$cod1} = 1;
	    }
	    #stop codon
	    elsif ($t[$CODONS->{$cod1}] eq '*' or $t[$CODONS->{$codons[$j]}] eq '*') {
		$results{$cod1}{$codons[$j]} = -1;
		$results{$codons[$j]}{$cod1} = -1;
	    }
	    # nc change
	    else {
		$results{$cod1}{$codons[$j]} = 0;
		$results{$codons[$j]}{$cod1} = 0;
	    }
	}
    }
    return %results;
}

=head2 dnds_pattern_number

 Title   : dnds_pattern_number
 Usage   : my $patterns = $stats->dnds_pattern_number($alnobj);
 Function: Counts the number of codons with no gaps in the MSA
 Returns : Number of codons with no gaps ('patterns' in PAML notation)
 Args    : A Bio::Align::AlignI compliant object such as a
            Bio::SimpleAlign object.

=cut

sub dnds_pattern_number{
    my ($self, $aln) = @_;
    return ($aln->remove_gaps->length)/3;
}

sub count_syn_sites {
    #counts the number of possible synonymous changes for sequence
    my ($seq, $synsite) = @_;
    __PACKAGE__->throw("not integral number of codons") if length($seq) % 3 != 0;
    my $S = 0;
    for (my $i = 0; $i< length($seq); $i+=3) {
	my $cod = substr($seq, $i, 3);
	next if $cod =~ /\-/;	#deal with alignment gaps
	$S +=  $synsite->{$cod}{'s'};
    }
    #print "S is $S\n";
    return $S;
}

	

sub get_syn_sites {
    #sub to generate lookup hash for the number of synonymous changes per codon
    my @nucs = qw(T C A G);
    my %raw_results;
    for my $i (@nucs) {
	for my $j (@nucs) {
	    for my $k (@nucs) {
		# for each possible codon
          	my $cod = "$i$j$k";
           	my $aa = $t[$CODONS->{$cod}];
		#calculate number of synonymous mutations vs non syn mutations
            	for my $i (qw(0 1 2)){
		    my $s = 0;
		    my $n = 3;
		    for my $nuc (qw(A T C G)) {
			next if substr ($cod, $i,1) eq $nuc;
			my $test = $cod;
			substr($test, $i, 1) = $nuc ;
			if ($t[$CODONS->{$test}] eq $aa) {
			    $s++;
			}
			if ($t[$CODONS->{$test}] eq '*') {
			    $n--;
			}	
		    }
		    $raw_results{$cod}[$i] = {'s' => $s ,
					      'n' => $n };
		}
		
	    } #end analysis of single codon
	}
    } #end analysis of all codons
    my %final_results;
    
    for my $cod (sort keys %raw_results) {
    	my $t = 0;
    	map{$t += ($_->{'s'} /$_->{'n'})} @{$raw_results{$cod}};
    	$final_results{$cod} = { 's'=>$t, 'n' => 3 -$t};
    }
    return \%final_results;
}

sub _make_codons {
#makes all codon combinations, returns array of them
    my @nucs = qw(T C A G);
    my @codons;
    for my $i (@nucs) {
        for my $j (@nucs) {
            for my $k (@nucs) {
            	push @codons, "$i$j$k";
	    }
	}
    }    
    return @codons;
}

sub get_codons {
 #generates codon translation look up table#
 my $x = 0;
 my  $CODONS = {};
 for my $codon (_make_codons) {
     $CODONS->{$codon} = $x;
     $x++;
 } 
 return $CODONS;
}

#########stats subs, can go in another module? Here for speed. ###
sub mean {
    my $ref = shift;
    my $el_num = scalar @$ref;
    my $tot = 0;
    map{$tot += $_}@$ref;
    return ($tot/$el_num);
}

sub variance {
    my $ref = shift;
    my $mean = mean($ref);
    my $sum_of_squares = 0;
    map{$sum_of_squares += ($_ - $mean) **2}@$ref;
    return $sum_of_squares;
}

sub sampling_variance {
    my $ref = shift;
    return variance($ref) / (scalar @$ref -1);
}

1;
