
#
# BioPerl module for Bio::Tools::SeqStats
#
# Cared for by
#
# Copyright Peter Schattner
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::SeqStats - Object holding statistics for one particular sequence

=head1 SYNOPSIS

    # build a primary nucleic acid or protein sequence object somehow
    # then build a statistics object from the sequence object

	$seqobj = Bio::PrimarySeq->new(-seq=>'ACTGTGGCGTCAACTG', 
				       -moltype = 'dna', -id = 'test');
	$seq_stats  =  Bio::Tools::SeqStats->new($seqobj);

    # obtain a hash of counts of each type of monomer (ie amino or
    # nucleic acid)

	$hash_ref = $seq_stats->count_monomers();  # eg for DNA sequence
	foreach $base ( sort keys $$hash_ref) {
	    print "Number of bases of type ",$base "= ",%$hash_ref{$base}"\n";
	}
    # or obtain the count directly without creating a new statistics object
	$hash_ref = Bio::Tools::SeqStats->count_monomers($seqobj);
	foreach $base ( sort keys $$hash_ref) {
	    print "Number of bases of type ",$base "= ",%$hash_ref{$base}"\n";
	}


    # obtain hash of counts of each type of codon in a nucleic acid sequence
	$hash_ref = $seq_stats-> count_codons();  # for nucleic acid sequence
    #  or
	$hash_ref = Bio::Tools::SeqStats->count_codons($seqobj);


    # Obtain the molecular weight of a sequence. Since the sequence
    # may contain ambiguous monomers, the molecular weight is returned
    # as a (reference to) a two element array containing greatest
    # lower bound (GLB) and least upper bound (LUB) of the molecular
    # weight

	$weight = $seq_stats->get_mol_wt();
    #  or
	$weight = Bio::Tools::SeqStats->get_mol_wt($seqobj);
	print "Molecular weight of sequence ", $seqobj->id(), " is greater than ", $$weight[0], " and less than " , $$weight[1], "\n";



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

SeqStats can be called in two distinct manners.  If only a single
computation is required on a given sequence object, the method can be
called easily using the SeqStats object directly:

	$weight = Bio::SeqStats->get_mol_wt($seqobj);

Alternately, if several computations will be required on a given
sequence object, an "instance" statistics object can be constructed
and used for the method calls:

	$seq_stats  =  Bio::SeqStats->new($seqobj);
	$monomers = $seq_stats->count_monomers();
	$codons = $seq_stats->count_codons();
	$weight = $seq_stats->get_mol_wt();

As currently implemented the object can return the following values
from a sequence:

      * The molecular weight of the sequence: get_mol_wt()
      * The number of each type of monomer present: count_monomers()
      * The number of each codon present in a nucleic acid sequence:
        count_codons()

Note that since sequences may contain ambiguous monomers (eg "M"
meaning "A" or "C" in a nucleic acid sequence), the method get_mol_wt
returns a two-element array containing the greatest lower bound and
least upper bound of the molecule. (For a sequence with no ambiguous
monomers, the two elements of the returned array will be equal.) The
method count_codons() handles ambiguous bases by simply counting all
ambiguous codons together and issuing a warning to that effect.


=head1 DEVELOPERS NOTES

Ewan moved it from Bio::SeqStats to Bio::Tools::SeqStats

=head1 STILL TO WRITE


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org               - General discussion
  http://bio.perl.org/MailList.html   - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR -  Peter Schattner

Email schattner@alum.mit.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Tools::SeqStats;
use vars qw(@ISA %Alphabets %Alphabets_strict %Weights 
	    $dna_weights $rna_weights $amino_weights
	    $dna_A_wt $dna_C_wt $dna_G_wt $dna_T_wt
	    $oxy_wt $rna_A_wt $rna_C_wt $rna_G_wt $rna_U_wt
	    $amino_A_wt $amino_C_wt $amino_D_wt
	    $amino_E_wt $amino_F_wt $amino_G_wt
	    $amino_H_wt $amino_I_wt $amino_K_wt
	    $amino_L_wt $amino_M_wt $amino_N_wt $amino_P_wt
	    $amino_Q_wt $amino_R_wt $amino_S_wt $amino_T_wt
	    $amino_V_wt $amino_W_wt $amino_Y_wt 
	    );
use strict;

use Bio::Root::RootI;
use Bio::Seq;
@ISA = qw(Bio::Root::RootI);
BEGIN { 
    %Alphabets =   
	('dna'     => [ "A","C","G","T","R","Y","M","K",
			"S","W","H","B","V","D","X", "N" ],
	 'rna'     => [ "A","C","G","U","R","Y","M","K","S",
			"W","H","B","V","D","X", "N" ],
	 'protein'    => [ "A","R","N","D","C","Q","E","G","H",
			   "I","L","K","M","F",
			   "P","S","T","W","X","Y","V",
			   "B","Z","*" ], # sac: added B, Z
	 );

# SAC: new strict alphabet: doesn't allow any ambiguity characters.
    %Alphabets_strict = (
			 'dna'     => [ "A","C","G","T" ],
			 'rna'     => [ "A","C","G","U" ],
			 'protein'    => [ "A","R","N","D","C",
					   "Q","E","G","H","I","L",
					   "K","M","F", "P","S","T","W",
					   "Y","V"],
			 );

# Extended Dna / Rna alphabet

# (includes symbols for nucleotide ambiguity)
# ------------------------------------------
# Symbol       Meaning      Nucleic Acid
# ------------------------------------------

    $dna_A_wt = 347;
    $dna_C_wt = 323;
    $dna_G_wt = 363;
    $dna_T_wt = 322;


    $dna_weights = {
	'A'              => [$dna_A_wt,$dna_A_wt], #  Adenine
	'C'              => [$dna_C_wt,$dna_C_wt], #  Cytosine
	'G'              => [$dna_G_wt,$dna_G_wt], #   Guanine
	'T'              => [$dna_T_wt,$dna_T_wt], #  Thymine
	'M'             => [$dna_C_wt,$dna_A_wt], # A or C
	'R'             => [$dna_A_wt,$dna_G_wt], # A or G
	'W'             => [$dna_T_wt,$dna_A_wt], # A or T
	'S'             => [$dna_C_wt,$dna_G_wt], # C or G
	'Y'             => [$dna_T_wt,$dna_C_wt], # C or T
	'K'             => [$dna_T_wt,$dna_G_wt], # G or T
	'V'             => [$dna_C_wt,$dna_G_wt], # A or C or G
	'H'             => [$dna_T_wt,$dna_A_wt], # A or C or T
	'D'             => [$dna_T_wt,$dna_G_wt], # A or G or T
	'B'             => [$dna_T_wt,$dna_G_wt], # C or G or T
	'X'             => [$dna_T_wt,$dna_G_wt], # G or A or T or C
	'N'             => [$dna_T_wt,$dna_G_wt], # G or A or T or C
    };

    $oxy_wt = 16;
    $rna_A_wt = $dna_A_wt + $oxy_wt;
    $rna_C_wt = $dna_C_wt + $oxy_wt;
    $rna_G_wt = $dna_G_wt + $oxy_wt;
    $rna_U_wt = 340;



    $rna_weights =  {
	'A'              => [$rna_A_wt,$rna_A_wt], #  Adenine
	'C'              => [$rna_C_wt,$rna_C_wt], #  Cytosine
	'G'              => [$rna_G_wt,$rna_G_wt], #   Guanine
	'U'              => [$rna_U_wt,$rna_U_wt], #   Uracil
	'M'             => [$rna_C_wt,$rna_A_wt], # A or C
	'R'             => [$rna_A_wt,$rna_G_wt], # A or G
	'W'             => [$rna_U_wt,$rna_A_wt], # A or U
	'S'             => [$rna_C_wt,$rna_G_wt], # C or G
	'Y'             => [$rna_C_wt,$rna_U_wt], # C or U
	'K'             => [$rna_U_wt,$rna_G_wt], # G or U
	'V'             => [$rna_C_wt,$rna_G_wt], # A or C or G
	'H'             => [$rna_C_wt,$rna_A_wt], # A or C or U
	'D'             => [$rna_U_wt,$rna_G_wt], # A or G or U
	'B'             => [$rna_C_wt,$rna_G_wt], # C or G or U
	'X'             => [$rna_C_wt,$rna_G_wt], # G or A or U or C
	'N'             => [$rna_C_wt,$rna_G_wt], # G or A or U or C
    };

#  IUPAC-IUB SYMBOLS FOR NUCLEOTIDE NOMENCLATURE:
#   Cornish-Bowden (1985) Nucl. Acids Res. 13: 3021-3030.


#  Amino Acid alphabet

# ------------------------------------------
# Symbol           Meaning   
# ------------------------------------------


    $amino_A_wt = 89.09;
    $amino_C_wt = 121.15;
    $amino_D_wt = 133.1;
    $amino_E_wt = 147.13;
    $amino_F_wt = 165.19;
    $amino_G_wt = 75.07;
    $amino_H_wt = 155.16;
    $amino_I_wt = 131.18;
    $amino_K_wt = 146.19;
    $amino_L_wt = 131.18;
    $amino_M_wt = 149.22;
    $amino_N_wt = 132.12;
    $amino_P_wt = 115.13;
    $amino_Q_wt = 146.15;
    $amino_R_wt = 174.21;
    $amino_S_wt = 105.09;
    $amino_T_wt = 119.12;
    $amino_V_wt = 117.15;
    $amino_W_wt = 204.22;
    $amino_Y_wt = 181.19;


    $amino_weights = {
	'A'     => [$amino_A_wt, $amino_A_wt], #    Alanine
	'B'      => [$amino_N_wt, $amino_D_wt],	#   Aspartic Acid, Asparagine
	'C'      => [$amino_C_wt, $amino_C_wt],	#   Cystine
	'D'         => [$amino_D_wt, $amino_D_wt], # Aspartic Acid
	'E'        => [$amino_E_wt, $amino_E_wt], # Glutamic Acid
	'F'        => [$amino_F_wt, $amino_F_wt], # Phenylalanine
	'G'        => [$amino_G_wt, $amino_G_wt], # Glycine
	'H'        => [$amino_H_wt, $amino_H_wt], # Histidine
	'I'        => [$amino_I_wt, $amino_I_wt], # Isoleucine
	'K'        => [$amino_K_wt, $amino_K_wt], # Lysine
	'L'        => [$amino_L_wt, $amino_L_wt], # Leucine
	'M'        => [$amino_M_wt, $amino_M_wt], # Methionine
	'N'        => [$amino_N_wt, $amino_N_wt], # Asparagine
	'P'        => [$amino_P_wt, $amino_P_wt], # Proline
	'Q'        => [$amino_Q_wt, $amino_Q_wt], # Glutamine
	'R'        => [$amino_R_wt, $amino_R_wt], # Arginine
	'S'        => [$amino_S_wt, $amino_S_wt], # Serine
	'T'        => [$amino_T_wt, $amino_T_wt], # Threonine
	'V'        => [$amino_V_wt, $amino_V_wt], # Valine
	'W'        => [$amino_W_wt, $amino_W_wt], # Tryptophan
	'X'        => [$amino_G_wt, $amino_W_wt], # Unknown
	'Y'        => [$amino_Y_wt, $amino_Y_wt], # Tyrosine
	'Z'        => [$amino_Q_wt, $amino_E_wt], # Glutamic Acid, Glutamine
    };


    %Weights =   (
		  'dna'     =>  $dna_weights,
		  'rna'     =>  $rna_weights,
		  'protein'    => $amino_weights,
		  );
}

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my $seqobj = shift (@args);
    unless  ($seqobj->isa("Bio::PrimarySeqI")) {
	die(" SeqStats works only on PrimarySeqI objects  \n");
    }
    $self->{'_seqref'} = $seqobj;
    $self->{'_is_strict'} = &_is_alphabet_strict($seqobj); # check the letters in the sequence
    return $self; 
}

=head2 count_monomers

 Title   : count_monomers
 Usage   : $rcount = $seq_stats->count_monomers(); 
           or $rcount = $seq_stats->Bio::SeqStats->($seqobj);
 Function: Counts the number of each type of monomer (amino acid or
           base) in the sequence.
 Example :

 Returns : Reference to a hash in which keys are letters of the
           genetic alphabet used and values are number of occurrences
           of the letter in the sequence.
 Args    : None or reference to sequence object

  Throws an exception if type of sequence is unknown (ie amino or
  nucleic)or if unknown letter in alphabet. Ambiguous elements are
  allowed.

=cut

sub count_monomers{
    my $rcount;
    my %count  = ();
    my $seqobj;
    my $_is_strict;
    my $element = '';
    my $_is_instance = 1 ;
    my $self = shift @_;
    my $object_argument = shift @_;

    # First we need to determine if the present object is an instance object or if
    # the sequence object has been passed as an argument
    if (defined $object_argument) {
	$_is_instance = 0;
    }

    # If we are using an instance object...
    if ($_is_instance) {
	if ($rcount = $self->{'_monomer_count'}) {
	    return $rcount;        # return count if previously calculated
	}
	$_is_strict =  $self->{'_is_strict'}; # retrieve "strictness"
        $seqobj =  $self->{'_seqref'};
    } else {
#  otherwise...
	$seqobj =  $object_argument;
	
	#  Following two lines lead to error in "throw" routine
	$seqobj->isa("Bio::PrimarySeqI") ||
	    die(" SeqStats works only on PrimarySeqI objects  \n");
        # is alphabet OK? Is it strict?
	$_is_strict =  &_is_alphabet_strict($seqobj); 
    }
	
    my $alphabet =  $_is_strict ? $Alphabets_strict{$seqobj->moltype} :
	$Alphabets{$seqobj->moltype}  ; # get array of allowed letters
	
 # convert everything to upper case to be safe
    my $seqstring = uc $seqobj->seq();  

   #  For each letter, count the number of times it appears in
   #  the sequence
  LETTER:
    foreach $element (@$alphabet) {
# skip terminator symbol which may confuse regex
	next LETTER if ($element eq '*'); 
	$count{$element} = ( $seqstring =~ s/$element/$element/g);
    }
    
    $rcount = \%count;
    
    if ($_is_instance) {
        # Save in case called again later
	$self->{'_monomer_count'} = $rcount;  
    }    
    return  $rcount;    
}


=head2  get_mol_wt

 Title   : get_mol_wt
 Usage   : $wt = $seqobj->get_mol_wt() or 
           $wt = Bio::SeqStats->get_mol_wt($seqobj);
 Function: Calculate molecular weight of sequence
 Example :

 Returns : Reference to two element array containing lower and upper
           bounds of molecule's molecular weight. If sequence contains
           no ambiguous elements, both entries in array are equal to
           molecular weight of molecule.  Args : None or reference to
           sequence object

 Throws an exception if type of sequence is unknown (ie not amino or
 nucleic)or if unknown letter in alphabet. Ambiguous elements are
 allowed.


=cut
#'

sub get_mol_wt {
    my ($self,$object_argument) = @_;

    my $seqobj;
    my $_is_strict;
    my $element = '';
    my $_is_instance = 1 ;
    my ($weight_array, $rcount);

    if (defined $object_argument) {
	$_is_instance = 0;
    }

    if ($_is_instance) {	
	if ($weight_array = $self->{'_mol_wt'}) {
	    return $weight_array; # return mol. weight if previously calculated
	}
        $seqobj =  $self->{'_seqref'};
        $rcount = $self->count_monomers();
    } else {
	$seqobj =  $object_argument;
	$seqobj->isa("Bio::PrimarySeqI") ||
	    die("Error: SeqStats works only on PrimarySeqI objects  \n");
	$_is_strict =  &_is_alphabet_strict($seqobj); # is alphabet OK?
        $rcount =  $self->count_monomers($seqobj);
    }

# We will also need to know what type of monomer we are dealing with

    my $moltype = $seqobj->moltype();

# In general,the molecular weight is bounded below by the sum of the
# weights of lower bounds of each alphabet symbol times the number of
# occurrences of the symbol in the sequence. A similar upper bound on
# the weight is also calculated.

#  Note that for "strict" (ie unambiguous) sequences there is an
# inefficiency since the upper bound = the lower bound (and is
# calculated twice).  However, this decrease in performance will be
# minor and leads to (IMO) significantly more readable code.

    my $weight_lower_bound = 0;
    my $weight_upper_bound = 0;
    my $weight_table =  $Weights{$moltype};
    my $water = 18.015;

    foreach $element (keys %$rcount) {
	$weight_lower_bound += $$rcount{$element} * $$weight_table{$element}->[0];
	$weight_upper_bound += $$rcount{$element} * $$weight_table{$element}->[1];
    }

    # Added by kdj: removal of H2O during peptide bond formation!
    $weight_lower_bound -= $water * ($seqobj->length - 1);
    $weight_upper_bound -= $water * ($seqobj->length - 1);
    $weight_lower_bound = sprintf("%.0f", $weight_lower_bound);
    $weight_upper_bound = sprintf("%.0f", $weight_upper_bound);

    $weight_array = [$weight_lower_bound, $weight_upper_bound];

    if ($_is_instance) {
	$self->{'_mol_wt'} = $weight_array; # Save in case called again later
    }
    return $weight_array;
}


=head2  count_codons

 Title   : count_codons
 Usage   : $rcount = $seqstats->count_codons (); or 
           $rcount = Bio::SeqStats->count_codons($seqobj);
 Function: Counts the number of each type of codons in a given frame
           for a dna or rna sequence.
 Example :
 Returns : Reference to a hash in which keys are codons of the genetic
           alphabet used and values are number of occurrences of the
           codons in the sequence. All codons with "ambiguous" bases
           are counted together.
 Args    : None or reference to sequence object

 Throws an exception if type of sequence is unknown or protein.


=cut

sub count_codons {
    my ($self, $object_argument) = @_;
    my $rcount = {};
    my $codon ;
    my $seqobj;
    my $_is_strict;
    my $element = '';
    my $_is_instance = 1 ;

    if (defined $object_argument) {
	$_is_instance = 0;
    }

    if ($_is_instance) {
	if ($rcount = $self->{'_codon_count'}) {
	    return $rcount;	# return count if previously calculated
	}
 	$_is_strict =  $self->{'_is_strict'}; # retrieve "strictness"
        $seqobj =  $self->{'_seqref'};
    } else {
	$seqobj =  $object_argument;
	$seqobj->isa("Bio::PrimarySeqI") ||
	    die(" Error: SeqStats works only on PrimarySeqI objects  \n");
	$_is_strict =  &_is_alphabet_strict($seqobj);
    }

# Codon counts only make sense for nucleic acid sequences
    my $moltype = $seqobj->moltype();

    unless ($moltype =~ /[dr]na/) {
	$seqobj->throw(" Codon counts only meaningful for dna or rna, not for $moltype sequences. \n");
    }

# If sequence contains ambiguous bases, warn that codons containing
# them will all be lumped together in the count.

    if (!$_is_strict ) {
	$seqobj->warn(" Sequence $seqobj contains ambiguous bases.  \n All codons with ambiguous bases will be added together in count.  \n");
    }

    my $seq = $seqobj->seq();

# Now step through the string by threes and count the codons

  CODON:
    while (length($seq) > 2) {
	$codon = substr($seq,0,3);
	$seq = substr($seq,3);
	if ($codon =~ /[^ACTGU]/) {
	    $$rcount{'ambiguous'}++; #lump together ambiguous codons
	    next CODON;
	}
	if (!defined $$rcount{$codon}) {
	    $$rcount{$codon}= 1 ;
	    next CODON;
	}
	$$rcount{$codon}++;	# default
    }    

    if ($_is_instance) {
	$self->{'_codon_count'} = $rcount; # Save in case called again later
    }    
    return $rcount;
}




=head2  _is_alphabet_strict

 Title   :   _is_alphabet_strict
 Usage   :
 Function: internal function to determine whether there are any
           ambiguous elements in the current sequence
 Example :
 Returns : 1 if strict alphabet is being used, 0 if ambiguous elements
           are present
 Args    :

  Throws an exception if type of sequence is unknown (ie amino or
  nucleic) or if unknown letter in alphabet. Ambiguous monomers are
  allowed.

=cut

sub  _is_alphabet_strict {
    my ($seqobj) = @_;
    my $moltype = $seqobj->moltype();
    # convert everything to upper case to be safe    
    my $seqstring = uc $seqobj->seq();   

# First we check if only the 'strict' letters are present in the
# sequence string If not, we check whether the remaining letters are
# ambiguous monomers or whether there are illegal letters in the
# string

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
    $seqobj->throw(" Alphabet not OK for $seqobj \n");

}
