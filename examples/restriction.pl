#!/usr/bin/perl -w

#-------------------------------------------------------------------------
# PROGRAM  : restriction.pl
# PURPOSE  : Simple driver for Bio::Tools::RestrictionEnzyme.pm
# AUTHOR   : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED  : 2 Jan 1997
#
# INSTALLATION
#    Edit the use lib "...." line to point the directory
#    containing your Bioperl modules.
#
# MODIFIED :
#   sac --- Tue Oct  6 09:22:16 1998
#     Assorted updates associated with the new Bio:: modules.
#----------------------------------------------------------

use lib "/home/sac";
use Bio::Seq;
use Bio::Tools::RestrictionEnzyme qw(@RE_available);

## Get string representing a DNA sequence and "purify" the DNA.
$sequence = get_sequence();  
$sequence =~ /[AGCT]{6,}/i;     
$& || die "\n*** CAN'T PROCESS SEQUENCE: Not DNA or not long enough\n$sequence\n";
$sequence = $&;

## Build sequence and restriction enzyme objects.
eval{
    $seq = new Bio::Seq(-ID  =>'test_seq', 
			-SEQ =>$sequence); 
    
    $re  = new Bio::Tools::RestrictionEnzyme(-NAME=>'EcoRI');
};
$@ && die "\n*** TROUBLE BUILDING OBJECTS:\n$@\n"; 

## Print some basic information about the new objects.
printf "Cutting sequence %s, length = %d bp\n", 
    $seq->id, $seq->seq_len;

printf "with enzyme %s, recognition sequence = %s", 
    $re->seq->id, $re->site;

printf " (%spalindromic)\n\n", 
    ($re->palindromic ? '' : 'not ');

## Cut the sequence with the restriction enzyme object.
@fragments = eval{ $re->cut_seq($seq); };

$@ && die "\n*** TROUBLE CUTTING:\n$@\n"; 

## Display the results and do some simply verification.
printf "%d fragments produced: ", scalar @fragments;
$total = 0;
foreach(@fragments) {
    printf "%d ",length $_; 
    $total += length $_; 
}

if($total == $seq->seq_len) {
    print "\nTotal fragment length equals size of original sequence.\n";  
} else {
    printf "\nRestriction Error: Total fragment length (%d) does not equal size of original sequence(%d)", $total, $seq->seq_len;  
}    


printf "\n%d enzymes are currently available:\n", scalar @RE_available;
for($i=1; $i<=@RE_available; $i++) {
    printf "%10s", $RE_available[$i-1];
    print "\n" if $i % 6 == 0;
}
print "\n";

@sixcutters = $re->available_list(6);
printf "\n%d enzymes have 6-base recognition sequences:\n", scalar @sixcutters;
for($i=1; $i<=@sixcutters; $i++) {
    printf "%10s", $sixcutters[$i-1];
    print "\n" if $i % 6 == 0;
}
print "\n\n";

printf "%s is%s available\n", 'HindIII', $re->is_available('HindIII') ? '' : ' not';
printf "%s is%s available\n", 'HindIV', $re->is_available('HindIV') ? '' : ' not';


exit(0);

#-----------------
sub get_sequence {
#-----------------
    ## This sequence has two EcoRI recognition sites (GAATTC).
    ## Cutting with an EcoRI restriction enzyme should produce 3 fragments.
    my $seq = <<DNA;
TTGCTGTGTGAGGCAGAACCTGCGGGGGCAGGGGCGGGCTGGTTCCCTGGCCAGCCATTGGCAGAGTCCG
CAGGCTAGGGCTGTCAATCATGCTGGCCGGCGTGGCCCCGCCTCCGCCGGCGCGGCCCCGCCTCCGCCGG
CGCACGTCTGGGACGCAAGGCGCCGTGGGGGCTGCCGGGACGGGTCCAAGATGGACGGCCGCTCAGGTTC
TGCTTTTACCTGCGGCCCAGAGCCCCATTCATTGCCCCGGTGCTGAGCGGCGCCGCGAGTCGGCCCGAGG
CCTCCGGGGACTGCCGTGCCGGGCGGGAATTCGCCATGGCGACCCTGGAAAAGCTGATGAAGGCCTTCGA
GTCCCTCAAGTCCTTCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG
CAGCAGCAGCAACAGCCGCCACCGCCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAGCCGCCGCCGC
AGGCACAGCCGCTGCTGCCTCAGCCGCAGCCGCCCCCGCCGCCGCCCCCGCCGCCACCCGGCCCGGCTGT
GGCTGAGGAGCCGCTGCACCGACCAAAGAAAGAACTTTCAGCTACCAAGAAAGACCGTGTGAATCATTGT
CTGACAATATGTGAAAACATAGTGGCACAGTCTGTCAGAAATTCTCCAGAATTCCAGAAACTTCTGGGCA
TCGCTATGGAACTTTTTCTGCTGTGCAGTGATGACGCAGAGTCAGATGTCAGGATGGTGGCTGACGAATG
CCTCAACAAAGTTATCAAAGCTTTGATGGATTCTAATCTTCCAAGGTTACAGCTCGAGCTCTATAAGGAA
DNA
    $seq =~ s/\s//g;
    $seq;
}
