# $Id$
#
# BioPerl module for Bio::Align::DNAStatistics
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Align::DNAStatistics - Calculate some statistics for a DNA alignment

=head1 SYNOPSIS

    use Bio::Align::DNAStatistics;
    use Bio::AlignIO;

    my $stats = new Bio::Align::PairwiseStatistics;
    my $alignin = new Bio::AlignIO(-format => 'emboss',
				   -file   => 't/data/insulin.water');
    my $jc = $stats->distance($aln, 'Jukes-Cantor');
    foreach my $r ( @$jc )  {
	print "\t";
	foreach my $r ( @$d ) {
	    print "$r\t";
	} 
	print "\n";
    }

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
  http://bioperl.org/bioperl-bugs/

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


package Bio::Align::DNAStatistics;
use vars qw(@ISA %DNAChanges @Nucleotides %NucleotideIndexes 
	    $GapChars $SeqCount $DefaultGapPenalty %DistanceMethods);
use strict;
use Bio::Align::PairwiseStatistics;
use Bio::Root::Root;

BEGIN {
    $GapChars = '(\.|\-)';
    @Nucleotides = qw(A G T C);
    $SeqCount = 2;
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
    %DistanceMethods = ( 'jc|jukes|jukes-cantor' => 'JukesCantor',
			 'f81'                   => 'F81',
			 'k2|k2p|k80|kimura'        => 'Kimura',
			 't92|tamura|tamura92'   => 'Tamura',
			 'f84'                   => 'F84',
			 'tajimanei|tajima-nei'  => 'TajimaNei' );
}

@ISA = qw( Bio::Root::Root Bio::Align::StatisticsI );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Align::DNAStatistics();
 Function: Builds a new Bio::Align::DNAStatistics object 
 Returns : Bio::Align::DNAStatistics
 Args    : none


=cut

sub new { 
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    
    $self->pairwise_stats( new Bio::Align::PairwiseStatistics());

    return $self;
}


=head2 distance

 Title   : distance
 Usage   : my $distance_mat = $stats->distance(-align  => $aln, 
		 			       -method => $method);
 Function: Calculates a distance matrix for all pairwise distances of
           sequences in an alignment.
 Returns : Array ref
 Args    : -align  => Bio::Align::AlignI object
           -method => String specifying specific distance method 
                      (implementing class may assume a default)

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
   return undef;
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
 Returns : ArrayRef of all pairwise distances of all sequence pairs in the alignment
 Args    : Bio::Align::AlignI of DNA sequences
           double - gap penalty


=cut

sub D_JukesCantor{
   my ($self,$aln,$gappenalty) = @_;
   return 0 unless $self->_check_arg($aln);
   $gappenalty = $DefaultGapPenalty unless defined $gappenalty;
   # ambiguities ignored at this point
   
   my (@seqs);
   foreach my $seq ( $aln->each_seq) {
       push @seqs, [ split(//,uc $seq->seq())];
   }
   my $seqct = scalar @seqs;
   my @DVals; 
   for(my $i = 1; $i <= $seqct; $i++ ) {
       for( my $j = $i+1; $j <= $seqct; $j++ ) {
	   my ($matrix,$pfreq,$gaps) = $self->_build_nt_matrix($seqs[$i-1],
							       $seqs[$j-1]);
	   # just want diagonals
	   my $m = ( $matrix->[0]->[0] + $matrix->[1]->[1] + 
		     $matrix->[2]->[2] + $matrix->[3]->[3] );
	   my $D = 1 - ( $m / ($aln->length - $gaps + ( $gaps * $gappenalty)));
	   my $d = (- 3 / 4) * log ( 1 - (4 * $D/ 3));
	   $DVals[$i]->[$j] = $DVals[$j]->[$i] = $d;
       }
   }
   return \@DVals;
}

=head2 D_F81

 Title   : D_F81
 Usage   : my $d = $stat->D_F81($aln)
 Function: Calculates D (pairwise distance) between 2 sequences in an 
           alignment using the Felsenstein 1981 distance model. 
 Returns : ArrayRef of a 2d array of all pairwise distances in the alignment
 Args    : Bio::Align::AlignI of DNA sequences


=cut

sub D_F81{
   my ($self,$aln) = @_;
   return 0 unless $self->_check_arg($aln);
   $self->throw("This isn't implemented yet - sorry");
}


# M Kimura, J. Mol. Evol., 1980, 16, 111.

=head2 D_Kimura

 Title   : D_Kimura
 Usage   : my $d = $stat->D_Kimura($aln)
 Function: Calculates D (pairwise distance) between 2 sequences in an 
           alignment using the Kimura 2 parameter model.
 Returns : ArrayRef of pairwise distances between all sequences in alignment
 Args    : Bio::Align::AlignI of DNA sequences


=cut

sub D_Kimura{
   my ($self,$aln) = @_;
   return 0 unless $self->_check_arg($aln);
   my $seqct = $aln->no_sequences;
   my @KVals;
   for( my $i = 1; $i <= $seqct; $i++ ) {
       for( my $j = $i+1; $j <= $seqct; $j++ ) {
	   my $pairwise = $aln->select_noncont($i,$j);
	   my $L = $self->pairwise_stats->number_of_comparable_bases($pairwise);
	   my $P = $self->transitions($pairwise) / $L;
	   my $Q = $self->transversions($pairwise) / $L;
	   
	   my $a = 1 / ( 1 - (2 * $P) - $Q);
	   my $b = 1 / ( 1 - 2 * $Q );
	   my $K = (1/2) * log ( $a ) + (1/4) * log($b);
	   $KVals[$i]->[$j] = $K;
	   $KVals[$j]->[$i] = $K;
       }
   }
   return \@KVals;
}

#  K Tamura, Mol. Biol. Evol. 1992, 9, 678.

=head2 D_Tamura

 Title   : D_Tamura
 Usage   :
 Function:
 Returns : 
 Args    :


=cut

sub D_Tamura{
   my ($self,$aln) = @_;
   my $seqct = $aln->no_sequences;
   my @KVals;
   for( my $i = 1; $i <= $seqct; $i++ ) {
       for( my $j = $i+1; $j <= $seqct; $j++ ) {
       }
   }
	   my $O = 0.25;
   my $t = 0;
   my $a = 0;
   my $b = 0;
   

   my $d = 4 * $O * ( 1 - $O ) * $a * $t  + 2 * $b * $t;
   return $d;
}

=head2 D_F84

 Title   : D_F84
 Usage   : my $d = $stat->D_F84($aln)
 Function: Calculates D (pairwise distance) between 2 sequences in an 
           alignment using the Felsenstein 1984 distance model. 
 Returns : Distance value
 Args    : Bio::Align::AlignI of DNA sequences
           double - gap penalty

=cut

sub D_F84{
   my ($self,$aln) = @_;
   return 0 unless $self->_check_arg($aln);
}

# Tajima and Nei, Mol. Biol. Evol. 1984, 1, 269.

=head2 D_TajimaNei

 Title   : D_TajimaNei
 Usage   : my $d = $stat->D_TajimaNei($aln)
 Function: Calculates D (pairwise distance) between 2 sequences in an 
           alignment using the TajimaNei 1984 distance model. 
 Returns : Distance value
 Args    : Bio::Align::AlignI of DNA sequences


=cut

sub D_TajimaNei{
   my ($self,$aln) = @_;
   $self->warn("The result from this method is not correct right now");
   my (@seqs);
   foreach my $seq ( $aln->each_seq) {
       push @seqs, [ split(//,uc $seq->seq())];
   }
   my $seqct = scalar @seqs;
   my @DVals; 
   for(my $i = 1; $i <= $seqct; $i++ ) {
       for( my $j = $i+1; $j <= $seqct; $j++ ) {
	   my ($matrix,$pfreq,$gaps) = $self->_build_nt_matrix($seqs[$i-1],
							       $seqs[$j-1]);
	   my $fij2;
	   my $slen = $aln->length - $gaps;
	   for( my $bs = 0; $bs < 4; $bs++ ) {
	       my $fi = 0;
	       map {$fi += $matrix->[$bs]->[$_] } 0..3;
	       my $fj = 0;
	       map { $fj += $matrix->[$_]->[$bs] } 0..3;
	       my $fij = ( $fi && $fj ) ? ($fi + $fj) /( 2 * $slen) : 0;
	       $fij2 += $fij**2;
	   }
	   my ($pair,$h) = (0,0);
	   for( my $bs = 0; $bs < 3; $bs++ ) {
	       for( my $bs1 = $bs+1; $bs1 <= 3; $bs1++ ) {
		   my $fij = $pfreq->[$pair++] / $slen;
		   if( $fij ) {
		       
		       my ($ci1,$ci2,$cj1,$cj2) = (0,0,0,0);

		       map { $ci1 += $matrix->[$_]->[$bs] } 0..3;
		       map { $cj1 += $matrix->[$bs]->[$_] } 0..3;
		       map { $ci2 += $matrix->[$_]->[$bs1] } 0..3;
		       map { $cj2 += $matrix->[$bs1]->[$_] } 0..3;
		       
		       $h += ( $fij*$fij / 2 ) / 
			   (  ( ( $ci1 + $cj1 ) / 2 * $slen ) *
			      ( ( $ci2 + $cj2 ) /2 * $slen ) 
			      );
		       $self->debug( "h is $h fij = $fij ci1 =$ci1 cj1=$cj1 ci2=$ci2 cj2=$cj2\n");
		   }
	       }
	   }
	   # just want diagonals first

	   my $m = ( $matrix->[0]->[0] + $matrix->[1]->[1] + 
		     $matrix->[2]->[2] + $matrix->[3]->[3] );
	   my $D = 1 - ( $m / $slen);

	   my $b = (1-$fij2+(($D**2)/$h)) / 2;
	   $self->debug("h is $h fij2 is $fij2 b is $b\n");

	   my $d = (-1 * $b) * log ( 1 - $D/ $b);
	   $DVals[$i]->[$j] = $DVals[$j]->[$i] = $d;
       }
   }
   return \@DVals;


}

# HKY -- HASEGAWA, M., H. KISHINO, and T. YANO. 1985
# Tamura and Nei 1993?
# GTR? 

=head2 K - sequence substitution methods

=cut

=head2 K_JukesCantor

 Title   : K_JukesCantor
 Usage   : my $k = $stats->K_JukesCantor($aln)
 Function: Calculates K - the number of nucleotide substitutions between 
           2 seqs - according to the Jukes-Cantor 1 parameter model
           This only involves the number of changes between two sequences.
 Returns : double
 Args    : Bio::Align::AlignI


=cut

sub K_JukesCantor{
   my ($self,$aln) = @_;
   return 0 unless $self->_check_arg($aln);
   my $seqct = $aln->no_sequences;
   my @KVals;
   for( my $i = 1; $i <= $seqct; $i++ ) {
       for( my $j = $i+1; $j <= $seqct; $j++ ) {
	   my $pairwise = $aln->select_noncont($i,$j);
	   my $L = $self->pairwise_stats->number_of_comparable_bases($pairwise);
	   my $N = $self->pairwise_stats->number_of_differences($pairwise);
	   my $p = $N / $L;   
	   my $K = - ( 3 / 4) * log ( 1 - (( 4 * $p) / 3 ));
	   $KVals[$i]->[$j] = $KVals[$j]->[$i] = $K;
       }
   }
   return \@KVals;
}

=head2 K_TajimaNei

 Title   : K_TajimaNei
 Usage   : my $k = $stats->K_TajimaNei($aln)
 Function: Calculates K - the number of nucleotide substitutions between 
           2 seqs - according to the Kimura 2 parameter model.
           This does not assume equal frequencies among all the nucleotides.
 Returns : ArrayRef of 2d matrix which contains pairwise K values for 
           all sequences in the alignment
 Args    : Bio::Align::AlignI

=cut

sub K_TajimaNei {
    my ($self,$aln) = @_;
    return 0 unless $self->_check_arg($aln);

    my @seqs;
    foreach my $seq ( $aln->each_seq) {
	push @seqs, [ split(//,uc $seq->seq())];
    }
    my @KVals;
    my $L = $self->pairwise_stats->number_of_comparable_bases($aln);	    
    my $seqct = scalar @seqs;
    for( my $i = 1; $i <= $seqct; $i++ ) {
	for( my $j = $i+1; $j <= $seqct; $j++ ) {
	    my (%q,%y);
	    my ($first,$second) = ($seqs[$i-1],$seqs[$j-1]);
	    
	    for (my $k = 0;$k<$aln->length; $k++ ) {	
		next if( $first->[$k] =~ /^$GapChars$/ ||
			 $second->[$k]  =~ /^$GapChars$/);
		
		$q{$second->[$k]}++;
		$q{$first->[$k]}++;
		if( $first->[$k] ne $second->[$k] ) {		
		    $y{$first->[$k]}->{$second->[$k]}++;
		}
	    }
	    
	    my $q_sum = 0;
	    foreach my $let ( @Nucleotides ) {              
		# ct is the number of sequences compared (2)
		# L is the length of the alignment without gaps
		# $ct * $L = total number of nt compared
		my $avg = $q{$let} / ( $SeqCount * $L );
		$q_sum += $avg**2;
	    }
	    my $b1 = 1 - $q_sum;
	    my $h = 0;
	    for( my $i = 0; $i <= 2; $i++ ) {
		for( my $j = $i+1; $j <= 3; $j++) {	    
		    $y{$Nucleotides[$i]}->{$Nucleotides[$j]} ||= 0;
		    $y{$Nucleotides[$j]}->{$Nucleotides[$i]} ||= 0;
		    my $x = ($y{$Nucleotides[$i]}->{$Nucleotides[$j]} + 
			     $y{$Nucleotides[$j]}->{$Nucleotides[$i]}) / $L;
		    $h += ($x ** 2) / ( 2 * $q{$Nucleotides[$i]} * 
					$q{$Nucleotides[$j]} );
		}
	    }
	    my $N = $self->pairwise_stats->number_of_differences($aln);
	    my $p = $N / $L;
	    my $b = ( $b1 + $p ** 2 / $h ) / 2;	    
	    my $K = - $b * log ( 1 - $p / $b );
	    $KVals[$i]->[$j] = $KVals[$j]->[$i] = $K;
	}
    }
    return \@KVals;
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
    my (@seqs,@tcount);
    foreach my $seq ( $aln->get_seq_by_pos(1), $aln->get_seq_by_pos(2) ) {
	push @seqs, [ split(//,$seq->seq())];
    }
    my ($first,$second) = @seqs;

    for (my $i = 0;$i<$aln->length; $i++ ) { 
	next if( $first->[$i]  =~ /^$GapChars$/ ||
		 $second->[$i]  =~ /^$GapChars$/);	
	if( $first->[$i] ne $second->[$i] ) {
	    foreach my $nt ( @{$type->{$first->[$i]}} ) {
		if( $nt eq $second->[$i]) {
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
    
    for( my $i = 0; $i < scalar @$seqa; $i++) {
	
	my ($ti,$tj) = ($seqa->[$i],$seqb->[$i]);
	$ti =~ tr/U/T/;
	$tj =~ tr/U/T/;

	if( $ti =~ /^$GapChars$/) { $gaps++; next; }
	if( $tj =~ /^$GapChars$/) { $gaps++; next }

	my $ti_index = $NucleotideIndexes{$ti};		
	my $tj_index = $NucleotideIndexes{$tj};	    

	if( ! defined $ti_index ) {
	    print "ti_index not defined for $ti\n";
	    next;
	}
	
	$basect_matrix->[$ti_index]->[$tj_index]++;
	
	if( $ti ne $tj ) {
	    $pfreq->[$NucleotideIndexes{join('',sort ($ti,$tj))}]++;
	}
    }
    return ($basect_matrix,$pfreq,$gaps);
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

1;
