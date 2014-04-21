#
# BioPerl module for Bio::Tools::SeqStats
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by
#
# Copyright Peter Schattner
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::SeqStats - Object holding statistics for one 
particular sequence

=head1 SYNOPSIS

  # build a primary nucleic acid or protein sequence object somehow
  # then build a statistics object from the sequence object

  $seqobj = Bio::PrimarySeq->new(-seq      => 'ACTGTGGCGTCAACTG',
                                 -alphabet => 'dna',
                                 -id       => 'test');
  $seq_stats  =  Bio::Tools::SeqStats->new(-seq => $seqobj);

  # obtain a hash of counts of each type of monomer
  # (i.e. amino or nucleic acid)
  print "\nMonomer counts using statistics object\n";
  $seq_stats  =  Bio::Tools::SeqStats->new(-seq=>$seqobj);
  $hash_ref = $seq_stats->count_monomers();  # e.g. for DNA sequence
  foreach $base (sort keys %$hash_ref) {
      print "Number of bases of type ", $base, "= ", 
         %$hash_ref->{$base},"\n";
  }

  # obtain the count directly without creating a new statistics object
  print "\nMonomer counts without statistics object\n";
  $hash_ref = Bio::Tools::SeqStats->count_monomers($seqobj);
  foreach $base (sort keys %$hash_ref) {
      print "Number of bases of type ", $base, "= ", 
         %$hash_ref->{$base},"\n";
  }


  # obtain hash of counts of each type of codon in a nucleic acid sequence
  print "\nCodon counts using statistics object\n";
  $hash_ref = $seq_stats-> count_codons();  # for nucleic acid sequence
  foreach $base (sort keys %$hash_ref) {
      print "Number of codons of type ", $base, "= ", 
         %$hash_ref->{$base},"\n";
  }

  #  or
  print "\nCodon counts without statistics object\n";
  $hash_ref = Bio::Tools::SeqStats->count_codons($seqobj);
  foreach $base (sort keys %$hash_ref) {
      print "Number of codons of type ", $base, "= ", 
         %$hash_ref->{$base},"\n";
  }

  # Obtain the molecular weight of a sequence. Since the sequence 
  # may contain ambiguous monomers, the molecular weight is returned 
  # as a (reference to) a two element array containing greatest lower 
  # bound (GLB) and least upper bound (LUB) of the molecular weight
  $weight = $seq_stats->get_mol_wt();
  print "\nMolecular weight (using statistics object) of sequence ", 
          $seqobj->id(), " is between ", $$weight[0], " and " ,
          $$weight[1], "\n";

  #  or
  $weight = Bio::Tools::SeqStats->get_mol_wt($seqobj);
  print "\nMolecular weight (without statistics object) of sequence ", 
        $seqobj->id(), " is between ", $$weight[0], " and " ,
        $$weight[1], "\n";

  # Calculate mean Kyte-Doolittle hydropathicity (aka "gravy" score)
  my $prot = Bio::PrimarySeq->new(-seq=>'MSFVLVAPDMLATAAADVVQIGSAVSAGS',
                                  -alphabet=>'protein');
  my $gravy = Bio::Tools::SeqStats->hydropathicity($seqobj);
  print "might be hydropathic" if $gravy > 1;  

=head1 DESCRIPTION

Bio::Tools::SeqStats is a lightweight object for the calculation of
simple statistical and numerical properties of a sequence. By
"lightweight" I mean that only "primary" sequences are handled by the
object.  The calling script needs to create the appropriate primary
sequence to be passed to SeqStats if statistics on a sequence feature
are required.  Similarly if a codon count is desired for a
frame-shifted sequence and/or a negative strand sequence, the calling
script needs to create that sequence and pass it to the SeqStats
object.

Nota that nucleotide sequences in bioperl do not strictly separate RNA
and DNA sequences. By convention, sequences from RNA molecules are
shown as is they were DNA. Objects are supposed to make the
distinction when needed. This class is one of the few where this
distinctions needs to be made. Internally, it changes all Ts into Us
before weight and monomer count.

SeqStats can be called in two distinct manners.  If only a single
computation is required on a given sequence object, the method can be
called easily using the SeqStats object directly:

  $weight = Bio::Tools::SeqStats->get_mol_wt($seqobj);

Alternately, if several computations will be required on a given
sequence object, an "instance" statistics object can be constructed
and used for the method calls:

  $seq_stats = Bio::Tools::SeqStats->new($seqobj);
  $monomers = $seq_stats->count_monomers();
  $codons = $seq_stats->count_codons();
  $weight = $seq_stats->get_mol_wt();
  $gravy = $seq_stats->hydropathicity();

As currently implemented the object can return the following values
from a sequence:

=over

=item *

The molecular weight of the sequence: get_mol_wt()

=item *

The number of each type of monomer present: count_monomers()

=item *

The number of each codon present in a nucleic acid sequence:
count_codons()

=item *

The mean hydropathicity ("gravy" score) of a protein:
hydropathicity()

=back

For DNA and RNA sequences single-stranded weights are returned. The
molecular weights are calculated for neutral, or not ionized,
nucleic acids. The returned weight is the sum of the
base-sugar-phosphate residues of the chain plus one weight of water to
to account for the additional OH on the phosphate of the 5' residue
and the additional H on the sugar ring of the 3' residue.  Note that
this leads to a difference of 18 in calculated molecular weights
compared to some other available programs (e.g. Informax VectorNTI).

Note that since sequences may contain ambiguous monomers (e.g. "M",
meaning "A" or "C" in a nucleic acid sequence), the method get_mol_wt
returns a two-element array containing the greatest lower bound and
least upper bound of the molecule. For a sequence with no ambiguous
monomers, the two elements of the returned array will be equal. The
method count_codons() handles ambiguous bases by simply counting all
ambiguous codons together and issuing a warning to that effect.


=head1 DEVELOPERS NOTES

Ewan moved it from Bio::SeqStats to Bio::Tools::SeqStats

Heikki made tiny adjustments (+/- 0.01 daltons) to amino acid
molecular weights to have the output match values in SWISS-PROT.

Torsten added hydropathicity calculation.

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
the bugs and their resolution.  Bug reports can be submitted the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Peter Schattner

Email schattner AT alum.mit.edu

=head1 CONTRIBUTOR - Torsten Seemann 

Email torsten.seemann AT infotech.monash.edu.au

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


package Bio::Tools::SeqStats;
use strict;
use vars qw(%Alphabets %Alphabets_strict $amino_weights
	    $rna_weights $dna_weights %Weights $amino_hydropathicity);
use Bio::Seq;
use base qw(Bio::Root::Root);

BEGIN {
	%Alphabets =   (
			 'dna'     => [ qw(A C G T R Y M K S W H B V D X N) ],
		    'rna'     => [ qw(A C G U R Y M K S W H B V D X N) ],
		    'protein' => [ qw(A R N D C Q E G H I L K M F U
									 P S T W X Y V B Z J O *) ], # sac: added B, Z
						);

# SAC: new strict alphabet: doesn't allow any ambiguity characters.
    %Alphabets_strict = (
			 'dna'     => [ qw( A C G T ) ],
			 'rna'     => [ qw( A C G U ) ],
			 'protein'    => [ qw(A R N D C Q E G H I L K M F U
					      P S T W Y V O) ],
			 );


#  IUPAC-IUB SYMBOLS FOR NUCLEOTIDE NOMENCLATURE:
#   Cornish-Bowden (1985) Nucl. Acids Res. 13: 3021-3030.

#  Amino Acid alphabet

# ------------------------------------------
# Symbol           Meaning
# ------------------------------------------

    my $amino_A_wt = 89.09;
    my $amino_C_wt = 121.15;
    my $amino_D_wt = 133.1;
    my $amino_E_wt = 147.13;
    my $amino_F_wt = 165.19;
    my $amino_G_wt = 75.07;
    my $amino_H_wt = 155.16;
    my $amino_I_wt = 131.17;
    my $amino_K_wt = 146.19;
    my $amino_L_wt = 131.17;
    my $amino_M_wt = 149.21;
    my $amino_N_wt = 132.12;
    my $amino_O_wt = 255.31;
    my $amino_P_wt = 115.13;
    my $amino_Q_wt = 146.15;
    my $amino_R_wt = 174.20;
    my $amino_S_wt = 105.09;
    my $amino_T_wt = 119.12;
    my $amino_U_wt = 168.06;
    my $amino_V_wt = 117.15;
    my $amino_W_wt = 204.23;
    my $amino_Y_wt = 181.19;


    $amino_weights = {
	'A'     => [$amino_A_wt, $amino_A_wt], # Alanine
	'B'     => [$amino_N_wt, $amino_D_wt], # Aspartic Acid, Asparagine
	'C'     => [$amino_C_wt, $amino_C_wt], # Cysteine
	'D'     => [$amino_D_wt, $amino_D_wt], # Aspartic Acid
	'E'     => [$amino_E_wt, $amino_E_wt], # Glutamic Acid
	'F'     => [$amino_F_wt, $amino_F_wt], # Phenylalanine
	'G'     => [$amino_G_wt, $amino_G_wt], # Glycine
	'H'     => [$amino_H_wt, $amino_H_wt], # Histidine
	'I'     => [$amino_I_wt, $amino_I_wt], # Isoleucine
	'J'     => [$amino_L_wt, $amino_I_wt], # Leucine, Isoleucine
	'K'     => [$amino_K_wt, $amino_K_wt], # Lysine
	'L'     => [$amino_L_wt, $amino_L_wt], # Leucine
	'M'     => [$amino_M_wt, $amino_M_wt], # Methionine
	'N'     => [$amino_N_wt, $amino_N_wt], # Asparagine
	'O'     => [$amino_O_wt, $amino_O_wt], # Pyrrolysine
	'P'     => [$amino_P_wt, $amino_P_wt], # Proline
	'Q'     => [$amino_Q_wt, $amino_Q_wt], # Glutamine
	'R'     => [$amino_R_wt, $amino_R_wt], # Arginine
	'S'     => [$amino_S_wt, $amino_S_wt], # Serine
	'T'     => [$amino_T_wt, $amino_T_wt], # Threonine
	'U'     => [$amino_U_wt, $amino_U_wt], # SelenoCysteine
	'V'     => [$amino_V_wt, $amino_V_wt], # Valine
	'W'     => [$amino_W_wt, $amino_W_wt], # Tryptophan
	'X'     => [$amino_G_wt, $amino_W_wt], # Unknown
	'Y'     => [$amino_Y_wt, $amino_Y_wt], # Tyrosine
	'Z'     => [$amino_Q_wt, $amino_E_wt], # Glutamic Acid, Glutamine
    };

    # Extended Dna / Rna alphabet
    use vars ( qw($C $O $N $H $P $water) );
    use vars ( qw($adenine   $guanine   $cytosine   $thymine   $uracil));
    use vars ( qw($ribose_phosphate   $deoxyribose_phosphate   $ppi));
    use vars ( qw($dna_A_wt   $dna_C_wt   $dna_G_wt  $dna_T_wt
		  $rna_A_wt   $rna_C_wt   $rna_G_wt   $rna_U_wt));
    use vars ( qw($dna_weights   $rna_weights   %Weights));

    $C = 12.01;
    $O = 16.00;
    $N = 14.01;
    $H = 1.01;
    $P = 30.97;
    $water = 18.015;

    $adenine = 5 * $C + 5 * $N + 5 * $H;
    $guanine = 5 * $C + 5 * $N + 1 * $O + 5 * $H;
    $cytosine = 4 * $C + 3 * $N + 1 * $O + 5 * $H;
    $thymine = 5 * $C + 2 * $N + 2 * $O + 6 * $H;
    $uracil = 4 * $C + 2 * $N + 2 * $O + 4 * $H;

    $ribose_phosphate = 5 * $C + 7 * $O + 9 * $H + 1 * $P;
    # neutral (unionized) form
    $deoxyribose_phosphate = 5 * $C + 6 * $O + 9 * $H + 1 * $P;

    # the following are single strand molecular weights / base
    $dna_A_wt = $adenine + $deoxyribose_phosphate - $water;
    $dna_C_wt = $cytosine + $deoxyribose_phosphate - $water;
    $dna_G_wt = $guanine + $deoxyribose_phosphate - $water;
    $dna_T_wt = $thymine + $deoxyribose_phosphate - $water;

    $rna_A_wt = $adenine + $ribose_phosphate - $water;
    $rna_C_wt = $cytosine + $ribose_phosphate - $water;
    $rna_G_wt = $guanine + $ribose_phosphate - $water;
    $rna_U_wt = $uracil + $ribose_phosphate - $water;

    $dna_weights = {
	'A'             => [$dna_A_wt,$dna_A_wt],            # Adenine
	'C'             => [$dna_C_wt,$dna_C_wt],            # Cytosine
	'G'             => [$dna_G_wt,$dna_G_wt],            # Guanine
	'T'             => [$dna_T_wt,$dna_T_wt],            # Thymine
	'M'             => [$dna_C_wt,$dna_A_wt],            # A or C
	'R'             => [$dna_A_wt,$dna_G_wt],            # A or G
	'W'             => [$dna_T_wt,$dna_A_wt],            # A or T
	'S'             => [$dna_C_wt,$dna_G_wt],            # C or G
	'Y'             => [$dna_C_wt,$dna_T_wt],            # C or T
	'K'             => [$dna_T_wt,$dna_G_wt],            # G or T
	'V'             => [$dna_C_wt,$dna_G_wt],            # A or C or G
	'H'             => [$dna_C_wt,$dna_A_wt],            # A or C or T
	'D'             => [$dna_T_wt,$dna_G_wt],            # A or G or T
	'B'             => [$dna_C_wt,$dna_G_wt],            # C or G or T
	'X'             => [$dna_C_wt,$dna_G_wt],            # G or A or T or C
	'N'             => [$dna_C_wt,$dna_G_wt],            # G or A or T or C
    };

    $rna_weights =  {
	'A'             => [$rna_A_wt,$rna_A_wt],            # Adenine
	'C'             => [$rna_C_wt,$rna_C_wt],            # Cytosine
	'G'             => [$rna_G_wt,$rna_G_wt],            # Guanine
	'U'             => [$rna_U_wt,$rna_U_wt],            # Uracil
	'M'             => [$rna_C_wt,$rna_A_wt],            # A or C
	'R'             => [$rna_A_wt,$rna_G_wt],            # A or G
	'W'             => [$rna_U_wt,$rna_A_wt],            # A or U
	'S'             => [$rna_C_wt,$rna_G_wt],            # C or G
	'Y'             => [$rna_C_wt,$rna_U_wt],            # C or U
	'K'             => [$rna_U_wt,$rna_G_wt],            # G or U
	'V'             => [$rna_C_wt,$rna_G_wt],            # A or C or G
	'H'             => [$rna_C_wt,$rna_A_wt],            # A or C or U
	'D'             => [$rna_U_wt,$rna_G_wt],            # A or G or U
	'B'             => [$rna_C_wt,$rna_G_wt],            # C or G or U
	'X'             => [$rna_C_wt,$rna_G_wt],            # G or A or U or C
	'N'             => [$rna_C_wt,$rna_G_wt],            # G or A or U or C
    };

    %Weights =   (
		  'dna'     =>  $dna_weights,
		  'rna'     =>  $rna_weights,
		  'protein' =>  $amino_weights,
		  );
	
	# Amino acid scale: Hydropathicity.
	# Ref: Kyte J., Doolittle R.F. J. Mol. Biol. 157:105-132(1982).
	# http://au.expasy.org/tools/pscale/Hphob.Doolittle.html
	
	$amino_hydropathicity = {
    A =>  1.800,  
    R => -4.500,  
    N => -3.500,  
    D => -3.500,  
    C =>  2.500,  
    Q => -3.500,  
    E => -3.500,  
    G => -0.400,  
    H => -3.200,  
    I =>  4.500,  
    L =>  3.800,  
    K => -3.900,  
    M =>  1.900,  
    F =>  2.800,  
    P => -1.600,  
    S => -0.800,  
    T => -0.700,  
    W => -0.900,  
    Y => -1.300,  
    V =>  4.200,  
	};

}

sub new {
	my($class,@args) = @_;
	my $self = $class->SUPER::new(@args);

	my ($seqobj) = $self->_rearrange([qw(SEQ)],@args);
	unless  ($seqobj->isa("Bio::PrimarySeqI")) {
		$self->throw("SeqStats works only on PrimarySeqI objects");
	}
	if ( !defined $seqobj->alphabet || 
		  !defined $Alphabets{$seqobj->alphabet}) {
		$self->throw("Must have a valid alphabet defined for seq (".
						 join(",",keys %Alphabets));
	}
	$self->{'_seqref'} = $seqobj;
	# check the letters in the sequence
	$self->{'_is_strict'} = _is_alphabet_strict($seqobj); 
	return $self;
}

=head2 count_monomers

 Title   : count_monomers
 Usage   : $rcount = $seq_stats->count_monomers();
           or $rcount = $seq_stats->Bio::Tools::SeqStats->($seqobj);
 Function: Counts the number of each type of monomer (amino acid or
	        base) in the sequence.
           Ts are counted as Us in RNA sequences.
 Example :
 Returns : Reference to a hash in which keys are letters of the
           genetic alphabet used and values are number of occurrences
           of the letter in the sequence.
 Args    : None or reference to sequence object
 Throws  : Throws an exception if type of sequence is unknown (ie amino
           or nucleic)or if unknown letter in alphabet. Ambiguous
           elements are allowed.

=cut

sub count_monomers{
	my %count  = ();
	my $seqobj;
	my $_is_strict;
	my $element = '';
	my $_is_instance = 1 ;
	my $self = shift @_;
	my $object_argument = shift @_;

	# First we need to determine if the present object is an instance
	# object or if the sequence object has been passed as an argument

	if (defined $object_argument) {
		$_is_instance = 0;
	}

	# If we are using an instance object...
	if ($_is_instance) {
		if ($self->{'_monomer_count'}) {
			return $self->{'_monomer_count'}; # return count if previously calculated
		}
		$_is_strict =  $self->{'_is_strict'}; # retrieve "strictness"
		$seqobj =  $self->{'_seqref'};
	} else {
		#  otherwise...
		$seqobj =  $object_argument;

		#  Following two lines lead to error in "throw" routine
		$seqobj->isa("Bio::PrimarySeqI") ||
		  $self->throw("SeqStats works only on PrimarySeqI objects");
		# is alphabet OK? Is it strict?
		$_is_strict =  _is_alphabet_strict($seqobj);
	}

	my $alphabet =  $_is_strict ? $Alphabets_strict{$seqobj->alphabet} :
	  $Alphabets{$seqobj->alphabet}  ; # get array of allowed letters

	# convert everything to upper case to be safe
	my $seqstring = uc $seqobj->seq();

	# Since T is used in RichSeq RNA sequences, do conversion locally
	$seqstring =~ s/T/U/g if $seqobj->alphabet eq 'rna';

	#  For each letter, count the number of times it appears in
	#  the sequence
 LETTER:
	foreach $element (@$alphabet) {
		# skip terminator symbol which may confuse regex
		next LETTER if $element eq '*';
		$count{$element} = ( $seqstring =~ s/$element/$element/g);
	}

	if ($_is_instance) {
		$self->{'_monomer_count'} = \%count;  # Save in case called again later
	}

	return \%count;
}

=head2  get_mol_wt

 Title   : get_mol_wt
 Usage   : $wt = $seqobj->get_mol_wt() or
           $wt = Bio::Tools::SeqStats ->get_mol_wt($seqobj);
 Function: Calculate molecular weight of sequence
           Ts are counted as Us in RNA sequences.
 Example :

 Returns : Reference to two element array containing lower and upper
           bounds of molecule molecular weight. For DNA and RNA
           sequences single-stranded weights are returned. If
           sequence contains no ambiguous elements, both entries in
           array are equal to molecular weight of molecule.
 Args    : None or reference to sequence object
 Throws  : Exception if type of sequence is unknown (ie not amino or
           nucleic) or if unknown letter in alphabet. Ambiguous
           elements are allowed.

=cut

sub get_mol_wt {
	my $seqobj;
	my $_is_strict;
	my $element = '';
	my $_is_instance = 1 ;
	my $self = shift @_;
	my $object_argument = shift @_;
	my ($weight_array, $rcount);

	if (defined $object_argument) {
		$_is_instance = 0;
	}

	if ($_is_instance) {
		if ($weight_array = $self->{'_mol_wt'}) {
			# return mol. weight if previously calculated
			return $weight_array;
		}
		$seqobj =  $self->{'_seqref'};
		$rcount = $self->count_monomers();
	} else {
		$seqobj =  $object_argument;
		$seqobj->isa("Bio::PrimarySeqI") ||
		  $self->throw("Error: SeqStats works only on PrimarySeqI objects");
		$_is_strict =  _is_alphabet_strict($seqobj); # is alphabet OK?
		$rcount =  $self->count_monomers($seqobj);
	}

	# We will also need to know what type of monomer we are dealing with
	my $moltype = $seqobj->alphabet();

	# In general,the molecular weight is bounded below by the sum of the
	# weights of lower bounds of each alphabet symbol times the number of
	# occurrences of the symbol in the sequence. A similar upper bound on
	# the weight is also calculated.

	# Note that for "strict" (i.e. unambiguous) sequences there is an
	# inefficiency since the upper bound = the lower bound and there are
	# two calculations.  However, this decrease in performance will be
	# minor and leads to significantly more readable code.

	my $weight_lower_bound = 0;
	my $weight_upper_bound = 0;
	my $weight_table =  $Weights{$moltype};
    my $total_res;
    
	# compute weight of all the residues
	foreach $element (keys %$rcount) {
		$weight_lower_bound += $$rcount{$element} * $$weight_table{$element}->[0];
		$weight_upper_bound += $$rcount{$element} * $$weight_table{$element}->[1];
        
        # this tracks only the residues used for counting MW
        $total_res += $$rcount{$element};
	}
	if ($moltype =~ /protein/) {
        # remove H2O during peptide bond formation.
    	$weight_lower_bound -= $water * ($total_res - 1);
    	$weight_upper_bound -= $water * ($total_res - 1);
	} else {
    	# Correction because phosphate of 5' residue has additional OH and
    	# sugar ring of 3' residue has additional H
    	$weight_lower_bound += $water;
    	$weight_upper_bound += $water;
	}

	$weight_lower_bound = sprintf("%.1f", $weight_lower_bound);
	$weight_upper_bound = sprintf("%.1f", $weight_upper_bound);

	$weight_array = [$weight_lower_bound, $weight_upper_bound];

	if ($_is_instance) {
		$self->{'_mol_wt'} = $weight_array;  # Save in case called again later
	}
	return $weight_array;
}


=head2  count_codons

 Title   : count_codons
 Usage   : $rcount = $seqstats->count_codons() or
           $rcount = Bio::Tools::SeqStats->count_codons($seqobj)
 Function: Counts the number of each type of codons for a dna or rna 
           sequence, starting at the 1st triple of the input sequence.
 Example :
 Returns : Reference to a hash in which keys are codons of the genetic
           alphabet used and values are number of occurrences of the
           codons in the sequence. All codons with "ambiguous" bases
           are counted together.
 Args    : None or sequence object
 Throws  : an exception if type of sequence is unknown or protein.

=cut

sub count_codons {
	my $rcount = {};
	my $codon ;
	my $seqobj;
	my $_is_strict;
	my $element = '';
	my $_is_instance = 1 ;
	my $self = shift @_;
	my $object_argument = shift @_;

	if (defined $object_argument) {
		$_is_instance = 0;
	}

	if ($_is_instance) {
		if ($rcount = $self->{'_codon_count'}) {
			return $rcount;        # return count if previously calculated
		}
		$_is_strict =  $self->{'_is_strict'}; # retrieve "strictness"
		$seqobj =  $self->{'_seqref'};
	} else {
		$seqobj =  $object_argument;
		$seqobj->isa("Bio::PrimarySeqI") ||
		  $self->throw("Error: SeqStats works only on PrimarySeqI objects");
		$_is_strict =  _is_alphabet_strict($seqobj);
	}

	# Codon counts only make sense for nucleic acid sequences
	my $alphabet = $seqobj->alphabet();

	unless ($alphabet =~ /[dr]na/i) {
		$seqobj->throw("Codon counts only meaningful for dna or rna, ".
							"not for $alphabet sequences.");
	}

	# If sequence contains ambiguous bases, warn that codons
	# containing them will all be lumped together in the count.

	if (!$_is_strict ) {
		$seqobj->warn("Sequence $seqobj contains ambiguous bases.".
		" All codons with ambiguous bases will be added together in count.")
                    if $self->verbose >= 0 ;
	}

	my $seq = $seqobj->seq();

	# Now step through the string by threes and count the codons

 CODON:
	while (length($seq) > 2) {
		$codon = uc substr($seq,0,3);
		$seq = substr($seq,3);
		if ($codon =~ /[^ACTGU]/i) {
			$$rcount{'ambiguous'}++; #lump together ambiguous codons
			next CODON;
		}
		if (!defined $$rcount{$codon}) {
			$$rcount{$codon}= 1 ;
			next CODON;
		}
		$$rcount{$codon}++;  # default
	}

	if ($_is_instance) {
		$self->{'_codon_count'} = $rcount;  # Save in case called again later
	}

	return $rcount;
}


=head2  hydropathicity

 Title   : hydropathicity
 Usage   : $gravy = $seqstats->hydropathicity(); or
           $gravy = Bio::Tools::SeqStats->hydropathicity($seqobj);

 Function: Calculates the mean Kyte-Doolittle hydropathicity for a
           protein sequence. Also known as the "gravy" score. Refer to 
           Kyte J., Doolittle R.F., J. Mol. Biol. 157:105-132(1982). 
 Example :
 Returns : float 
 Args    : None or reference to sequence object

 Throws  : an exception if type of sequence is not protein.

=cut

sub hydropathicity {
	my $seqobj;
	my $_is_strict;
	my $element = '';
	my $_is_instance = 1 ;
	my $self = shift @_;
	my $object_argument = shift @_;

	if (defined $object_argument) {
		$_is_instance = 0;
	}

	if ($_is_instance) {
		if (my $gravy = $self->{'_hydropathicity'}) {
			return $gravy;        # return value if previously calculated
		}
		$_is_strict =  $self->{'_is_strict'}; # retrieve "strictness"
		$seqobj =  $self->{'_seqref'};
	} else {
		$seqobj =  $object_argument;
		$seqobj->isa("Bio::PrimarySeqI") ||
		  $self->throw("Error: SeqStats works only on PrimarySeqI objects");
		$_is_strict =  _is_alphabet_strict($seqobj);
	}
	
	# hydropathicity not menaingful for empty sequences
	unless ($seqobj->length() > 0) {
	  $seqobj->throw("hydropathicity not defined for zero-length sequences");
        }

	# hydropathicity only make sense for protein sequences
	my $alphabet = $seqobj->alphabet();

	unless ($alphabet =~ /protein/i) {
		$seqobj->throw("hydropathicity only meaningful for protein, ".
							"not for $alphabet sequences.");
	}

	# If sequence contains ambiguous bases, warn that codons
	# containing them will all be lumped together in the count.

	unless ($_is_strict ) {
		$seqobj->throw("Sequence $seqobj contains ambiguous amino acids. ".
		"Hydropathicity can not be caculated.")
	}

	my $seq = $seqobj->seq();

	# Now step through the string and add up the hydropathicity values

    my $gravy = 0;
    for my $i ( 0 .. length($seq) ) {
       my $codon = uc(substr($seq,$i,1));
       $gravy += $amino_hydropathicity->{$codon}||0; # table look-up
    }
    $gravy /= length($seq);


	if ($_is_instance) {
		$self->{'_hydropathicity'} = $gravy;  # Save in case called again later
	}

	return $gravy;
}


=head2  _is_alphabet_strict

 Title   :  _is_alphabet_strict
 Usage   :
 Function: internal function to determine whether there are
           any ambiguous elements in the current sequence
 Example :
 Returns : 1 if strict alphabet is being used,
           0 if ambiguous elements are present
 Args    :

 Throws  : an exception if type of sequence is unknown (ie amino or
           nucleic) or if unknown letter in alphabet. Ambiguous
           monomers are allowed.

=cut

sub _is_alphabet_strict {

	my ($seqobj) = @_;
	my $moltype = $seqobj->alphabet();

	# convert everything to upper case to be safe
	my $seqstring = uc $seqobj->seq();

	# Since T is used in RichSeq RNA sequences, do conversion locally
	$seqstring =~ s/T/U/g if $seqobj->alphabet eq 'rna';

	# First we check if only the 'strict' letters are present in the
	# sequence string If not, we check whether the remaining letters
	# are ambiguous monomers or whether there are illegal letters in
	# the string

	# $alpha_array is a ref to an array of the 'strictly' allowed letters
	my $alpha_array =   $Alphabets_strict{$moltype} ;

	# $alphabet contains the allowed letters in string form
	my $alphabet = join ('', @$alpha_array) ;
	unless ($seqstring =~ /[^$alphabet]/)  {
		return 1 ;
	}

	# Next try to match with the alphabet's ambiguous letters
	$alpha_array =   $Alphabets{$moltype} ;
	$alphabet = join ('', @$alpha_array) ;

	unless ($seqstring =~ /[^$alphabet]/)  {
		return 0 ;
	}

	# If we got here there is an illegal letter in the sequence
	$seqobj->throw("Alphabet not OK for $seqobj");
}

=head2   _print_data

 Title   : _print_data
 Usage   : $seqobj->_print_data() or Bio::Tools::SeqStats->_print_data();
 Function: Displays dna / rna parameters (used for debugging)
 Returns : 1
 Args    : None

Used for debugging.

=cut

sub _print_data {

    print "\n adenine  = :  $adenine \n";
    print "\n guanine  = :  $guanine \n";
    print "\n cytosine = :  $cytosine \n";
    print "\n thymine  = :  $thymine \n";
    print "\n uracil   = :  $uracil \n";

    print "\n dna_A_wt = :  $dna_A_wt \n";
    print "\n dna_C_wt = :  $dna_C_wt \n";
    print "\n dna_G_wt = :  $dna_G_wt \n";
    print "\n dna_T_wt = :  $dna_T_wt \n";

    print "\n rna_A_wt = :  $rna_A_wt \n";
    print "\n rna_C_wt = :  $rna_C_wt \n";
    print "\n rna_G_wt = :  $rna_G_wt \n";
    print "\n rna_U_wt = :  $rna_U_wt \n";

    return 1;
}

1;
