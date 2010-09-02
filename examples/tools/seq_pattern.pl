#!/usr/bin/env perl

#-----------------------------------------------------------------------------
# PROGRAM : seq_pattern.pl 
# PURPOSE : This is a simple driver used to test the Bio::Tools::SeqPattern.pm
#           module for working with sequence patterns (regexps that recognize
#           nucleotide or peptide sequences).
# AUTHOR  : Steve Chervitz (sac@bioperl.org)
# CREATED : 28 Aug 1997
# USAGE   : seq_pattern.pl -h
# COMMENTS: 
#    This is a driver script for the Bio::Tools::SeqPattern.pm Bioperl module
#    that can be used for working with both nucleotide and peptide sequence and 
#    offers features such as:
#    
#      -- generate reverse complement of sequence pattern
#      -- ensuring pattern has no invalid characters
#      -- untainting pattern
#      -- expanding ambiguity codes. 
#    
#     Functionality is not yet complete but it may be of use as-is. 
# 
# INSTALLATION
#    Edit the use lib "...." line to point the directory
#    containing your Bioperl modules.
#
# DOCUMENTATION:
#  http://genome-www.stanford.edu/perlOOP/bioperl/lib/Bio/Tools/SeqPattern.pm.html
#
#-----------------------------------------------------------------------------

use lib "/Users/steve/lib/perl";
use Bio::Tools::SeqPattern ();
use Getopt::Std;

$opt_h = 0;
$opt_n = 0;
$opt_p = 0;
$opt_r = 0;

getopts('hnprv:');
$pat = $ARGV[0] || '';

$opt_h and die <<"QQ_USAGE_QQ";

Usage: seq_pattern.pl [-n|p|r|h] 'REGEXP' 

 regexp : full-regular expression for a nucleotide or peptide sequence.
          Must be listed *after* one of the following options:
 -n     : interpret regexp as a nucleotide pattern.
 -p     : interpret regexp as a peptide pattern.
 -r     : output only the reverse complement of the nucleotide pattern.
 -h     : print usage.
 
QQ_USAGE_QQ


## Nucleotide test patterns (most are based on actual patterns submitted by users):

%nucpat = (1 =>'YR...CG(CCG){5,7}CG[^G]TN{10,}[SA]{4}NN(ACA){2,}GCGTTT.{20,40}GT>',
	   2 =>'cggnnn[ta][ta][ta]n{3,5}[ta][ta][ta]nnnccg',
	   3 =>'<ATGX{6,10}RTTRTT',
	   4 =>'cggnnnwwwn{3,5}wwwnnnccg',
	   5 =>'(CCCCT)N{1,200}(agggg)N{1,200}(agggg)', 
	   6 =>'cccct{2,}', 
	   7 =>'(a){10,40}', 
	   8 =>'(cag){36,}', 
	   9 =>'rgaatgx{2,}ygtttca(cag){5,}',
	   10 =>'yattgtt(n){20,80}yattgtt', 
	   11 =>'yattgtt(aca){20,80}yattgtt', 
	   12 =>'TATAAAN{30,100}[AT][CAT][AT]YCAAR[CAT][AT][CAT]', 
	   13 =>'TGACTC[N]{1,300}TGACTC',
	   14 =>'TGACTCN*GAGTCAN*GAGTCAN*TGACTC',
	   15 =>'TGACTC(TCA)*GAGTCA',
	   16 =>'TGACTCN*GAG(TCA)*GAGTCA',
	   17 =>'[at][at]ttcacatgy',
	   );

%peppat = (1 =>'<X{10,}[WFY]XXXDN[BK][ST]Z{5,}>',
	   2 =>'<x{10,40}[gas]x[gasct]x*[gascdn]x[gas]x{0,10}[bst]{8,}x{0,8}>', 
	   );

#----------------------
# Main

if($opt_r) {
    print Bio::Tools::SeqPattern->new(-SEQ =>$pat, -TYPE =>'Dna')->revcom->str,"\n";
    
} else {
    test_nuc($pat) if ($opt_n and !$opt_p);
    test_pep($pat) if ($opt_p and !$opt_n);
    (test_nuc($pat), test_pep($pat)) if !($opt_p or $opt_n);
}
    
exit 0;

#----------------------

sub test_nuc {
# Create nucleotide pattern object:
    my $pat = shift;
    $pat ||= $nucpat{9};

    $npat = new Bio::Tools::SeqPattern(-seq =>$pat, -type =>'Dna');

    print "\nNucleotide Pattern:\n";
    print "-----------------------\n";
    printf "%18s: %s\n", 'Type',    $npat->type;
    printf "%18s: %s\n", 'Original',$npat->str;
    printf "%18s: %s\n", 'Expanded', $npat->expand;
    printf "%18s: %s\n", 'Reverse-Comp', $npat->revcom->str;
    printf "%18s: %s\n", 'Rev-Comp+Expanded', $npat->revcom(1)->str; # Hate this syntax. May change.
    print "\n";
}


sub test_pep {
# Create peptide pattern object:
    my $pat = shift;
    $pat ||= $peppat{1};

    $ppat = new Bio::Tools::SeqPattern(-seq =>$pat, -type =>'Amino');

    print "\nPeptide Pattern:\n";
    print "-----------------------\n";
    printf "%18s: %s\n", 'Type',    $ppat->type;
    printf "%18s: %s\n", 'Original',$ppat->str;
    printf "%18s: %s\n", 'Expanded', $ppat->expand;
    print "\n";
}    






