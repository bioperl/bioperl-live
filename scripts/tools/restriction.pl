#!/usr/bin/perl -w

#-------------------------------------------------------------------------
# PROGRAM  : restriction.pl
# PURPOSE  : Demonstrate how to use Bio::Tools::RestrictionEnzyme
# AUTHOR   : Steve Chervitz (sac@bioperl.org>
# CREATED  : 2 Jan 1997
#
# INSTALLATION
#    Edit the use lib "...." line to point the directory
#    containing your Bioperl modules.
#
# MODIFIED :
#   sac --- Sat Feb 22 00:14:13 2002
#     Added example showing how to cut a seq at runs of Ns
#     Added example trying to create a bogus enzyme.
#     Using PrimarySeq instead of Seq.
#   sac --- Tue Oct  6 09:22:16 1998
#     Assorted updates associated with the new Bio:: modules.
#----------------------------------------------------------

use lib ".";
use Bio::PrimarySeq;
use Bio::Tools::RestrictionEnzyme qw(@RE_available);

## Get string representing a DNA sequence and "purify" the DNA.
$sequence = get_sequence();  
$sequence =~ /[AGCT]{6,}/i;     
$& || die "\n*** CAN'T PROCESS SEQUENCE: Not DNA or not long enough\n$sequence\n";
$sequence = $&;

## Build sequence and restriction enzyme objects.
eval{
    $seq = new Bio::PrimarySeq(-ID  =>'test_seq', 
			       -SEQ =>$sequence); 
    
    $re  = new Bio::Tools::RestrictionEnzyme(-name=>'EcoRI');
};
$@ && die "\n*** TROUBLE BUILDING OBJECTS:\n$@\n"; 

## Print some basic information about the new objects.
printf "Cutting sequence %s, length = %d bp\n", 
    $seq->id, $seq->length;

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

if($total == $seq->length) {
    print "\nTotal fragment length equals size of original sequence.\n";  
} else {
    printf "\nRestriction Error: Total fragment length (%d) does not equal size of original sequence(%d)", $total, $seq->length;  
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


# Special restriction enzyme named 'N' will match actual N's in the
# sequence to be "digested". Normally, N's represent [AGCT].
# Specify the number of N's you want to match. This cannot be a regexp.
# It must be a specific number of N's.
$re_n  = new Bio::Tools::RestrictionEnzyme(-name=>'N--NNNNN', -make=>'custom');

printf "\nSpecial enzyme %s, recognition sequence = %s\n",
    $re_n->seq->id, $re_n->site;

print " This 'enzyme' recognizes runs of 1 or more N's in the sequence being 'digested'\n";


$seq_n = new Bio::PrimarySeq(-ID  =>'test_seq', 
			     -SEQ =>get_sequence_with_n() ); 

## Cut the sequence with the restriction enzyme object.
@locs= @{$re_n->cut_locations($seq_n)};

## Display the results and do some simply verification.
printf "Cut locations: @locs\n";

## Cut the sequence with the restriction enzyme object.
@fragments = eval{ $re_n->cut_seq($seq_n); };

$@ && die "\n*** TROUBLE CUTTING:\n$@\n"; 

## Display the results and do some simply verification.
printf "%d fragments produced: ", scalar @fragments;
$total = 0;
foreach(@fragments) {
    printf "%d ",length $_; 
    print "\n$_\n";
    $total += length $_; 
}

print "\n\n";

print "Attempt to build a non-custom enzyme that is not available.\n";
print "Should throw an exception.\n\n";

$re  = new Bio::Tools::RestrictionEnzyme(-name=>'FooBAR');

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

#----------------------
sub get_sequence_with_n {
#---------------------
    ## This sequence has internal runs of N's
    ## Cutting with a special 'N' restriction enzyme with 5 N's
    ## should produce 3 fragments.
    my $seq = <<DNA;
TTGCTGTGTGAGGCAGAACCTGCGGGGGCAGGGGCGGGCTGGTTCCCTGGCCAGCCATTGGCAGAGTCCG
CAGGCTAGGGCTGTCAATCATGCTGGCCGGCGTGGCCCCGCCTCCGCCGGCGCGGCCCCGCCTCCGCCGG
CGCACGTCTGGGACGCAAGGCGCCGTNNNNNCTGCCGGGACGGGTCCAAGATGGACGGCCGCTCAGGTTC
TGCTTTTACCTGCGGCCCAGAGCCCCATTCATTGCCCCGGTGCTGAGCGGCGCCGCGAGTCGGCCCGAGG
CCTCCGGGGACTGCCGTGCCGGGCGGGAATTCGCCATGGCGACCCTGGAAAAGCTGATGAAGGCCTTCGA
GTCCCTCAAGTCCTTCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG
CAGCAGCAGCAACAGCCGCCACCGCCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAGCCGCCGCCGC
AGGCACAGCCGCTGCTGCCTCAGCCGCAGCCGCCCCCGCCGCCGCCCCCGCCGCCACCCGGCCCGGCTGT
GGCTGAGGAGCCGCTGCACCGACCAAAGAAAGAACTTTCAGCTACCAAGAAAGACCGTGTGAATCATTGT
CTGACAATATGTGAAAACATAGTGGCACAGTCTGTCAGAAATTCTCCAGAATTCCAGAAACTTCTGGGCA
TCGCTATGGAACTTTTTCTGCTGTGCAGTGATGACGCAGAGTCAGATGNNNNNNNNNNNNNNNACGAATG
CCTCAACAAAGTTATCAAAGCTTTGATGGATTCTAATCTTCCAAGGTTACAGCTCGAGCTCTATAAGGAA
DNA
    $seq =~ s/\s//g;
    $seq;
}
