#
# BioPerl module for Bio::Align::ProteinStatistics
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Align::ProteinStatistics - Calculate Protein Alignment statistics (mostly distances)

=head1 SYNOPSIS

  use Bio::Align::ProteinStatistics;
  use Bio::AlignIO;
  my $in = Bio::AlignIO->new(-format => 'fasta',
			    -file   => 'pep-104.fasaln');
  my $aln = $in->next_aln;

  my $pepstats = Bio::Align::ProteinStatistics->new();
  $kimura = $protstats->distance(-align => $aln,
			         -method => 'Kimura');
  print $kimura->print_matrix;


=head1 DESCRIPTION

This object is for generating various statistics from a protein
alignment.  Mostly it is where pairwise protein distances can be
calculated.

=head1 REFERENCES 

D_Kimura - Kimura, M. 1983. The Neutral Theory of Molecular Evolution. CUP, 
           Cambridge.

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

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Align::ProteinStatistics;
use vars qw(%DistanceMethods $Precision $DefaultGapPenalty);
use strict;

use Bio::Align::PairwiseStatistics;
use Bio::Matrix::PhylipDist;

%DistanceMethods = ('kimura|k' => 'Kimura',
		    );
$Precision = 5;
$DefaultGapPenalty = 0;

use base qw(Bio::Root::Root Bio::Align::StatisticsI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Align::ProteinStatistics->new();
 Function: Builds a new Bio::Align::ProteinStatistics object 
 Returns : an instance of Bio::Align::ProteinStatistics
 Args    :


=cut

sub new {
  my($class,@args) = @_;

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

=cut

sub distance{
   my ($self,@args) = @_;
   my ($aln,$method) = $self->_rearrange([qw(ALIGN METHOD)],@args);
   if( ! defined $aln || ! ref ($aln) || ! $aln->isa('Bio::Align::AlignI') ) { 
       $self->throw("Must supply a valid Bio::Align::AlignI for the -align parameter in distance");
   }
   $method ||= 'Kimura';
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


=head2 D_Kimura

 Title   : D_Kimura
 Usage   : my $matrix = $pepstats->D_Kimura($aln);
 Function: Calculate Kimura protein distance (Kimura 1983) which 
           approximates PAM distance
           D = -ln ( 1 - p - 0.2 * p^2 )
 Returns : L<Bio::Matrix::PhylipDist>
 Args    : L<Bio::Align::AlignI>


=cut

# Kimura, M. 1983. The Neutral Theory of Molecular Evolution. CUP, Cambridge.

sub D_Kimura{
   my ($self,$aln) = @_;
   return 0 unless $self->_check_arg($aln);
   # ambiguities ignored at this point
   my (@seqs,@names,@values,%dist);
   my $seqct = 0;
   foreach my $seq ( $aln->each_seq) {
       push @names, $seq->display_id;
       push @seqs, uc $seq->seq();
       $seqct++;
   }
   my $len = $aln->length;
   my $precisionstr = "%.$Precision"."f";

   for( my $i = 0; $i < $seqct-1; $i++ ) {
       # (diagonals) distance is 0 for same sequence
       $dist{$names[$i]}->{$names[$i]} = [$i,$i];
       $values[$i][$i] = sprintf($precisionstr,0);
       for( my $j = $i+1; $j < $seqct; $j++ ) {
	   my ($scored,$match) = (0,0);
	   for( my $k=0; $k < $len; $k++ ) {
	       my $m1 = substr($seqs[$i],$k,1);
	       my $m2 = substr($seqs[$j],$k,1);
	       if( $m1 ne '-' && $m2 ne '-' ) {
		   # score is number of scored bases (alignable bases)
		   # it could have also come from 
		   # my $L = $self->pairwise_stats->number_of_comparable_bases($pairwise);
		   # match is number of matches weighting ambiguity bases
		   # as well
		   $match += _check_ambiguity_protein($m1,$m2);
		   $scored++;
	       }
	   }
	   # From Felsenstein's PHYLIP documentation:
	   # This is very quick to do but has some obvious
	   # limitations. It does not take into account which amino
	   # acids differ or to what amino acids they change, so some
	   # information is lost. The units of the distance measure
	   # are fraction of amino acids differing, as also in the
	   # case of the PAM distance. If the fraction of amino acids
	   # differing gets larger than 0.8541 the distance becomes
	   # infinite.

	   my $D = 1 - ( $match / $scored );
	   if( $D < 0.8541 ) {
	       $D = - log ( 1 - $D - (0.2 * ($D ** 2)));
	       $values[$j][$i] = $values[$i][$j] = sprintf($precisionstr,$D);
	   } else { 
	       $values[$j][$i] = $values[$i][$j] = '    NaN';
	   }
	   # fwd and rev lookup
	   $dist{$names[$i]}->{$names[$j]} = [$i,$j];
	   $dist{$names[$j]}->{$names[$i]} = [$i,$j];	   

           # (diagonals) distance is 0 for same sequence
	   $dist{$names[$j]}->{$names[$j]} = [$j,$j];	   
	   $values[$j][$j] = sprintf($precisionstr,0); 

       }
   }
   return Bio::Matrix::PhylipDist->new(-program => 'bioperl_PEPstats',
				       -matrix  => \%dist,
				       -names   => \@names,
				       -values  => \@values); 
   
}

# some methods from EMBOSS distmat
sub _check_ambiguity_protein
{
    my ($t1,$t2) = @_;
    my $n = 0;

    if( $t1 ne 'X' && $t1 eq $t2 ) { 
        $n = 1.0;
    } elsif(  (($t1 eq 'B' && $t2 =~ /[DN]/ ) ||
	       ($t2 eq 'B' && $t1 =~ /[DN]/ )) ||
	      
	      (($t1 eq 'Z' && $t2 =~ /[EQ]/) ||
	       ($t2 eq 'Z' && $t1 =~ /[EQ]/ ))) {
        $n = 0.5;
    } elsif ( $t1 eq 'X' && $t2 eq 'X' ) {
        $n = 0.0025;
    } elsif(  $t1 eq 'X' || $t2 eq 'X' ) {
        $n = 0.05;
    }
    return $n;
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

sub _check_arg {
    my($self,$aln ) = @_;
    if( ! defined $aln || ! $aln->isa('Bio::Align::AlignI') ) {
	$self->warn("Must provide a Bio::Align::AlignI compliant object to Bio::Align::DNAStatistics");
	return 0;
    } elsif( $aln->get_seq_by_pos(1)->alphabet ne 'protein' ) { 
	$self->warn("Must provide a protein alignment to Bio::Align::ProteinStatistics, you provided a " . $aln->get_seq_by_pos(1)->alphabet);
	return 0;
    }
    return 1;
}

1;
